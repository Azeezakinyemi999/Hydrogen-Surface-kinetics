"""Fixed-gas steady-state solver (documentation: Stages 3-5).

The gas composition is held fixed and only the surface coverages are solved for.
Method 1 ('direct') is Cantera's Newton solve; Method 3 ('reactor') integrates the
surface in time inside a reactor whose gas volume is huge compared with the catalyst
area. Every attempt starts from a fresh surface object (see Model.new_surface), and
every result is accepted only if it passes check_convergence().

Every function takes the Model first; solver settings come from model.cfg.solver.
"""
import os
import time

import cantera as ct
import numpy as np
import pandas as pd

from .conditions import condition_key
from .runner import run_jobs


def _cantera_error(e):
    lines = [l for l in str(e).splitlines() if l.strip() and '***' not in l]
    return ' | '.join(lines[:3])


def clean_surface(model):
    return {model.site: 1.0}


def set_state(model, T, P, X, theta0, multipliers=None, h_shifts=None):
    """Fresh surface object; gas set to (T, P, X); coverages set to theta0.
    The gas is only SET here, never integrated. Returns a copy of the gas mole
    fractions, to check later that the gas did not change."""
    gas = model.gas
    surf = model.new_surface(multipliers, h_shifts)
    gas.TPX = T, P, X
    surf.TP = T, P
    unknown = [name for name in theta0 if name not in surf.species_names]
    if unknown:
        raise ValueError(f'unknown surface species in theta0: {unknown}')
    total = sum(theta0.values())
    if not np.isclose(total, 1.0, atol=1e-8):
        raise ValueError(f'starting coverages sum to {total:.10f}, not 1')
    surf.coverages = theta0
    assert gas.T == surf.T == T, 'gas and surface temperatures differ'
    assert np.isclose(gas.P, P), 'gas pressure not set'
    X_total = sum(X.values())
    gas_X = gas.mole_fraction_dict()
    for name, x in X.items():
        assert np.isclose(gas_X.get(name, 0.0), x / X_total), f'gas mole fraction of {name} not set'
    return gas.X.copy()


def check_convergence(model, gas_X_start):
    """Is the current surface state a converged steady state? Each species' net rate must
    satisfy |dtheta/dt| <= abs_tol + rel_gross_tol * (its gross rate)."""
    s, surf, gas = model.cfg.solver, model.surf, model.gas
    to_theta = model.site_size / model.gamma
    dtheta_dt = surf.net_production_rates[model.surf_kin_idx] * to_theta
    gross = (surf.creation_rates[model.surf_kin_idx] + surf.destruction_rates[model.surf_kin_idx]) * to_theta
    allowed = s['abs_tol'] + s['rel_gross_tol'] * gross
    k_worst = int(np.argmax(np.abs(dtheta_dt) / allowed))
    residual = abs(dtheta_dt[k_worst])
    tof_fwd = model.probe_rates()[0] / model.gamma
    sum_error = abs(surf.coverages.sum() - 1.0)
    min_coverage = surf.coverages.min()
    gas_drift = np.abs(gas.X - gas_X_start).max()
    checks = {
        'residual_ok': bool(np.all(np.abs(dtheta_dt) <= allowed)),
        'site_balance_ok': bool(sum_error <= s['sum_tol']),
        'non_negative_ok': bool(min_coverage >= s['neg_tol']),
        'gas_fixed_ok': bool(gas_drift <= s['gas_tol']),
    }
    return {
        'converged': all(checks.values()), **checks,
        'residual_1_per_s': residual,
        'residual_species': surf.species_names[k_worst],
        'residual_allowed_1_per_s': allowed[k_worst],
        'residual_over_tof': residual / tof_fwd if tof_fwd > 0 else np.inf,
        'tof_fwd_1_per_s': tof_fwd,
        'sum_error': sum_error,
        'min_coverage': min_coverage,
        'gas_drift': gas_drift,
    }


def method_direct(model, gas_X_start):
    """Method 1: Cantera's direct steady-state solve. Returns (ok, error, seconds, t_end=None)."""
    t0 = time.time()
    try:
        model.surf.advance_coverages_to_steady_state()
        return True, None, time.time() - t0, None
    except ct.CanteraError as e:
        return False, _cantera_error(e), time.time() - t0, None


def fixed_gas_reactor(model):
    """A ReactorNet holding the current surface in a huge isothermal gas volume.
    Call set_state() first. Gas-phase reactions must be switched off by the caller."""
    s = model.cfg.solver
    reactor = ct.IdealGasReactor(model.gas, energy='off')
    reactor.volume = 1.0
    ct.ReactorSurface(model.surf, reactor, A=reactor.volume * s['reactor_area_per_volume'])
    net = ct.ReactorNet([reactor])
    net.rtol, net.atol = s['reactor_rtol'], s['reactor_atol']
    return net


def method_reactor(model, gas_X_start):
    """Method 3: integrate in time, half-decade checkpoints from 1e-8 s to reactor_t_max.
    Stops at the first checkpoint where the residual test passes AND the TOF and the
    coverages have stopped changing since the previous checkpoint.
    Returns (ok, error, seconds, t_end)."""
    s = model.cfg.solver
    t0 = time.time()
    net = fixed_gas_reactor(model)
    model.gas.set_multiplier(0.0)
    previous = None
    try:
        n_checkpoints = int(round(2 * (np.log10(s['reactor_t_max']) + 8))) + 1
        for t in np.logspace(-8, np.log10(s['reactor_t_max']), n_checkpoints):
            net.advance(t)
            diag = check_convergence(model, gas_X_start)
            current = (diag['tof_fwd_1_per_s'], model.surf.coverages.copy())
            if previous is not None and diag['residual_ok']:
                tof_change = abs(current[0] - previous[0]) / max(abs(previous[0]), 1e-300)
                theta_change = np.abs(current[1] - previous[1]).max()
                if tof_change < s['plateau_tof_rtol'] and theta_change < s['plateau_theta_atol']:
                    return True, None, time.time() - t0, t
            previous = current
        return False, f"no plateau by t = {s['reactor_t_max']:g} s", time.time() - t0, s['reactor_t_max']
    except ct.CanteraError as e:
        return False, _cantera_error(e), time.time() - t0, net.time
    finally:
        model.gas.set_multiplier(1.0)


METHODS = {'direct': method_direct, 'reactor': method_reactor}


def default_methods(model, T):
    return ('reactor',) if T < model.cfg.solver['direct_min_T'] else ('direct', 'reactor')


def solve_steady_state(model, T, P, X, theta0, methods=None, multipliers=None, h_shifts=None):
    """Steady-state coverages for a fixed gas (T, P, X), starting from theta0.
    Tries each method in turn (fresh state each time) and returns the first result that
    passes check_convergence(). Returns (coverages array or None, info dict)."""
    methods = methods or default_methods(model, T)
    t_start = time.time()
    attempts = []
    for name in methods:
        gas_X_start = set_state(model, T, P, X, theta0, multipliers, h_shifts)
        ok, err, secs, t_end = METHODS[name](model, gas_X_start)
        diag = check_convergence(model, gas_X_start) if ok else None
        if ok and not diag['converged']:
            failed = [k for k, v in diag.items() if k.endswith('_ok') and not v]
            ok, err = False, f'finished but failed checks: {failed}'
        attempts.append({'method': name, 'ok': ok, 'seconds': secs, 't_end': t_end, 'error': err})
        if ok:
            info = {'method': name, **diag, 'total_seconds': time.time() - t_start, 'attempts': attempts}
            return model.surf.coverages.copy(), info
    return None, {'method': 'failed', 'converged': False,
                  'total_seconds': time.time() - t_start, 'attempts': attempts}


def describe(info):
    """One-line summary of a solve_steady_state() result."""
    tries = ', '.join(f"{a['method']} {'OK' if a['ok'] else 'failed'} ({a['seconds']:.0f} s)"
                      for a in info['attempts'])
    if not info['converged']:
        return f"FAILED after {info['total_seconds']:.0f} s [{tries}]"
    return (f"{info['method']} in {info['total_seconds']:.0f} s [{tries}] | "
            f"TOF fwd {info['tof_fwd_1_per_s']:.4e} 1/s | residual {info['residual_1_per_s']:.1e} 1/s")


# --- Two starting surfaces (Stage 5) ------------------------------------------------
def coverages_to_dict(model, theta):
    """Coverage array -> {species: theta}, clipped at 0 and renormalised to sum to 1."""
    theta = np.clip(theta, 0.0, None)
    theta = theta / theta.sum()
    return {name: c for name, c in zip(model.species_names, theta) if c > 0}


def probe_covered_start(model):
    f = model.cfg.solver['probe_covered_fraction']
    return {model.probe_adsorbate: f, model.site: 1.0 - f}


_impurity_start_cache = {}


def impurity_covered_start(model, T, P, X):
    """Start B for a case with impurities: a clean surface exposed to the impurity gas
    ALONE (probe gas removed, same impurity partial pressures) for pre_exposure_time.
    Returns (theta0 dict or None, info dict)."""
    s = model.cfg.solver
    x_total = sum(X.values())
    imp_X = {sp: x for sp, x in X.items() if sp != model.probe_gas}
    imp_P = P * sum(imp_X.values()) / x_total
    key = (model.mechanism_sha256, T, tuple(sorted((sp, round(P * x / x_total, 6)) for sp, x in imp_X.items())))
    if key in _impurity_start_cache:
        return _impurity_start_cache[key]
    t0 = time.time()
    set_state(model, T, imp_P, imp_X, clean_surface(model))
    net = fixed_gas_reactor(model)
    model.gas.set_multiplier(0.0)
    try:
        net.advance(s['pre_exposure_time'])
        theta0, err = coverages_to_dict(model, model.surf.coverages), None
    except ct.CanteraError as e:
        theta0, err = None, _cantera_error(e)
    finally:
        model.gas.set_multiplier(1.0)
    info = {'pre_exposure_s': s['pre_exposure_time'], 'seconds': time.time() - t0, 'error': err}
    _impurity_start_cache[key] = (theta0, info)
    return theta0, info


def compare_starts(model, theta_a, info_a, theta_b, info_b):
    """Do the two starting surfaces reach the same steady state?"""
    s = model.cfg.solver
    if theta_a is None or theta_b is None:
        return {'starts_agree': None, 'max_theta_diff': np.nan, 'tof_rel_diff': np.nan}
    max_theta_diff = np.abs(theta_a - theta_b).max()
    tof_a, tof_b = info_a['tof_fwd_1_per_s'], info_b['tof_fwd_1_per_s']
    tof_rel_diff = abs(tof_a - tof_b) / max(abs(tof_a), abs(tof_b), 1e-300)
    agree = bool(max_theta_diff <= s['multi_ss_theta_atol'] and tof_rel_diff <= s['multi_ss_tof_rtol'])
    return {'starts_agree': agree, 'max_theta_diff': max_theta_diff, 'tof_rel_diff': tof_rel_diff}


def solve_condition(model, row, multipliers=None, h_shifts=None, both_starts=True):
    """Solve one condition row from start A (clean) and, if both_starts, start B
    (probe-covered when there are no impurities, otherwise impurity-pre-exposed)."""
    T, P, X = row['T_K'], row['P_total_bar'] * 1e5, row['X']
    theta_a, info_a = solve_steady_state(model, T, P, X, clean_surface(model),
                                         multipliers=multipliers, h_shifts=h_shifts)
    out = {'theta_a': theta_a, 'info_a': info_a, 'species': list(model.species_names)}
    if not both_starts:
        return out
    if set(X) == {model.probe_gas}:
        start_b, start_b_info = probe_covered_start(model), None
    else:
        start_b, start_b_info = impurity_covered_start(model, T, P, X)
    if start_b is None:
        theta_b, info_b = None, {'method': 'failed', 'converged': False, 'attempts': [],
                                 'total_seconds': 0.0, 'note': 'impurity pre-exposure failed'}
    else:
        theta_b, info_b = solve_steady_state(model, T, P, X, start_b,
                                             multipliers=multipliers, h_shifts=h_shifts)
    out.update({'theta_b': theta_b, 'info_b': info_b, 'start_b': start_b, 'start_b_info': start_b_info,
                **compare_starts(model, theta_a, info_a, theta_b, info_b)})
    return out


# --- Many conditions in parallel (Stage 6) ------------------------------------------
def condition_jobs(conditions, prefix=''):
    """Conditions table -> runner jobs (one per row)."""
    jobs = []
    for _, row in conditions.iterrows():
        row = row.to_dict()
        jobs.append({'key': prefix + condition_key(row), 'row': row})
    return jobs


def run_conditions(model, conditions, cache_subdir='cache', prefix='', both_starts=True, **run_overrides):
    """Solve every row of a conditions table from both starts, in parallel, with caching.
    run_overrides (workers, job_timeout, recompute) override cfg.run. Returns {key: result};
    each result has 'row' as well as 'job'."""
    run = model.cfg.section('run', **run_overrides)

    def job_func(job):
        return solve_condition(model, job['row'], both_starts=both_starts)

    def note(res):
        text = f"A {res['info_a']['method']}"
        if both_starts:
            text += f", B {res['info_b']['method']} | starts agree: {res['starts_agree']}"
        return text

    results = run_jobs(job_func, condition_jobs(conditions, prefix),
                       cache_dir=os.path.join(model.cfg.output_dir, cache_subdir),
                       context=(model.fingerprint_context(), both_starts),
                       workers=run['workers'], timeout=run['job_timeout'], recompute=run['recompute'],
                       describe=note)
    for res in results.values():
        res['row'] = res['job']['row']
    return results


def status_table(results):
    """One line per condition: did everything converge, and do the starts agree?"""
    rows = []
    for key, res in results.items():
        r = res['row']
        done = res.get('status') == 'done'
        has_b = done and 'info_b' in res
        rows.append({
            'key': key, 'basis': r['basis'], 'case': r['case'], 'T_K': r['T_K'],
            'status': res.get('status'),
            'method_A': res['info_a']['method'] if done else None,
            'method_B': res['info_b']['method'] if has_b else None,
            'converged_A': bool(res['info_a']['converged']) if done else False,
            'converged_B': bool(res['info_b']['converged']) if has_b else None,
            'starts_agree': res.get('starts_agree') if has_b else None,
            'tof_fwd_A': res['info_a'].get('tof_fwd_1_per_s') if done else None,
            'wall_s': round(res['wall_seconds'], 1),
        })
    return pd.DataFrame(rows)
