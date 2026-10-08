"""Parallel jobs with time limits and a per-job cache (documentation: Stage 6).

A job is a dict with a readable 'key' and whatever the job function needs. Each job runs
in its own forked process (the child inherits the loaded Model), is stopped after the
time limit, and when it finishes is saved straight away to
    <cache_dir>/<safe key>_<hash>.pkl
The hash covers the job itself plus the model context (mechanism file contents, phase
names, solver settings), so a saved result is never reused after any of those change.
"""
import hashlib
import multiprocessing as mp
import os
import pickle
import re
import time
import traceback

import numpy as np


def canonical(obj):
    """Deterministic representation for hashing: dicts sorted, floats to 12 significant digits."""
    if isinstance(obj, dict):
        return tuple(sorted((str(k), canonical(v)) for k, v in obj.items()))
    if isinstance(obj, (list, tuple)):
        return tuple(canonical(v) for v in obj)
    if isinstance(obj, (float, np.floating)):
        return float(f'{float(obj):.12g}')
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.ndarray):
        return canonical(obj.tolist())
    return obj


def safe_name(key):
    return re.sub(r'[^A-Za-z0-9.+-]+', '_', key.replace(' | ', '__')).strip('_')


def cache_path(cache_dir, job, context):
    spec = {k: v for k, v in job.items() if k != 'key'}
    digest = hashlib.md5(repr((canonical(spec), canonical(context))).encode()).hexdigest()[:12]
    return os.path.join(cache_dir, f"{safe_name(job['key'])}_{digest}.pkl")


def _worker(func, job, conn):
    try:
        conn.send({'status': 'done', **func(job)})
    except Exception:
        conn.send({'status': 'error', 'error': traceback.format_exc()})
    finally:
        conn.close()


def default_workers():
    return max(1, (os.cpu_count() or 2) // 2)


def run_jobs(func, jobs, cache_dir, context=(), workers=None, timeout=1800.0, recompute=False,
             describe=None, log=print):
    """Run func(job) -> dict for every job; returns {key: result}. Each result also has
    'status' ('done' / 'error' / 'timeout'), 'job' and 'wall_seconds'.
    describe(result) -> str is printed for each finished job (optional)."""
    keys = [j['key'] for j in jobs]
    if len(set(keys)) != len(keys):
        raise ValueError('job keys must be unique')
    workers = workers or default_workers()
    os.makedirs(cache_dir, exist_ok=True)
    ctx = mp.get_context('fork')
    results, pending = {}, []
    for job in jobs:
        path = cache_path(cache_dir, job, context)
        if os.path.exists(path) and not recompute:
            with open(path, 'rb') as f:
                results[job['key']] = pickle.load(f)
        else:
            pending.append(job)
    log(f'{len(results)} loaded from cache, {len(pending)} to solve with {workers} workers')

    n_total, running, t_batch = len(jobs), {}, time.time()
    while pending or running:
        while pending and len(running) < workers:
            job = pending.pop(0)
            parent_conn, child_conn = ctx.Pipe(duplex=False)
            proc = ctx.Process(target=_worker, args=(func, job, child_conn), daemon=True)
            proc.start()
            child_conn.close()
            running[job['key']] = (proc, parent_conn, time.time(), job)
        for key in list(running):
            proc, conn, t_start, job = running[key]
            elapsed = time.time() - t_start
            if conn.poll():
                try:
                    res = conn.recv()
                except EOFError:
                    res = {'status': 'error', 'error': 'worker closed the pipe without a result'}
                proc.join()
            elif not proc.is_alive():
                res = {'status': 'error', 'error': f'worker exited with code {proc.exitcode}'}
            elif elapsed > timeout:
                proc.terminate()
                proc.join()
                res = {'status': 'timeout', 'error': f'stopped after {timeout:g} s'}
            else:
                continue
            res.update({'job': job, 'wall_seconds': elapsed})
            results[key] = res
            del running[key]
            if res['status'] == 'done':
                with open(cache_path(cache_dir, job, context), 'wb') as f:
                    pickle.dump(res, f)
                note = describe(res) if describe else ''
                log(f'  [{len(results):3d}/{n_total}] {key}: {elapsed:5.0f} s  {note}')
            else:
                log(f"  [{len(results):3d}/{n_total}] {key}: {res['status'].upper()} after {elapsed:.0f} s - "
                    f"{res['error'].strip().splitlines()[-1][:150]}")
        time.sleep(0.2)
    log(f'finished in {time.time() - t_batch:.0f} s')
    return {k: results[k] for k in keys}
