"""Command line:

    python -m hydrogen_surface_kinetics check  runs/<name>/input.py
    python -m hydrogen_surface_kinetics run    runs/<name>/input.py [--only a,b] [--workers N]
                                               [--timeout S] [--recompute]
    python -m hydrogen_surface_kinetics status runs/<name>/input.py
"""
import argparse
import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')                                     # no windows from the command line

from .config import ANALYSIS_ORDER, ConfigError, load_config   # noqa: E402
from .model import Model, ModelError                           # noqa: E402


def cmd_check(args):
    cfg = load_config(args.input)
    model = Model(cfg)
    print(model.summary())
    from .conditions import main_conditions
    conds = main_conditions(model)
    print(f'\n{len(conds)} main conditions:')
    print(conds.drop(columns='X').to_string(index=False))
    print(f'\nanalyses: {cfg.analyses}\noutput: {cfg.output_dir}\ninput OK')


def cmd_run(args):
    from . import analyses
    only = args.only.split(',') if args.only else None
    bad = [n for n in only or [] if n not in ANALYSIS_ORDER]
    if bad:
        raise ConfigError(f'unknown analyses {bad}; choose from {ANALYSIS_ORDER}')
    analyses.run(args.input, only=only, workers=args.workers, job_timeout=args.timeout,
                 recompute=True if args.recompute else None)


def cmd_status(args):
    import pickle
    cfg = load_config(args.input)
    out = cfg.output_dir
    if not os.path.isdir(out):
        print(f'no output yet ({out})')
        return
    print(f'output: {out}')
    for d in sorted(glob.glob(os.path.join(out, 'cache*'))):
        files = glob.glob(os.path.join(d, '*.pkl'))
        bad = 0
        for p in files:
            with open(p, 'rb') as f:
                res = pickle.load(f)
            info = res.get('info_a') or res.get('info') or {}
            if info and not info.get('converged', True):
                bad += 1
        print(f'  {os.path.basename(d):20s} {len(files):5d} cached results, {bad} not converged')
    tables = sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, '*.csv')))
    figures = sorted(os.path.basename(p) for p in glob.glob(os.path.join(out, '*.png')))
    print(f'  tables:  {", ".join(tables) or "-"}\n  figures: {", ".join(figures) or "-"}')


def main(argv=None):
    p = argparse.ArgumentParser(prog='hydrogen_surface_kinetics', description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest='command', required=True)
    for name, fn, text in (('check', cmd_check, 'load the input and mechanism, list the conditions'),
                           ('run', cmd_run, 'run the analyses'),
                           ('status', cmd_status, 'what is cached and saved')):
        sp = sub.add_parser(name, help=text)
        sp.add_argument('input', help='path to the run input.py')
        sp.set_defaults(func=fn)
        if name == 'run':
            sp.add_argument('--only', help=f'comma-separated subset of {",".join(ANALYSIS_ORDER)}')
            sp.add_argument('--workers', type=int, help='parallel processes')
            sp.add_argument('--timeout', type=float, help='time limit per job (s)')
            sp.add_argument('--recompute', action='store_true', help='ignore cached results')
    args = p.parse_args(argv)
    try:
        args.func(args)
    except (ConfigError, ModelError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 2
    return 0
