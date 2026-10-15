"""Command line:

    python -m hydrogen_surface_kinetics check  runs/<name>/settings.py
    python -m hydrogen_surface_kinetics run    runs/<name>/settings.py [--only a,b] [--workers N]
                                               [--timeout S] [--recompute] [--no-report]
    python -m hydrogen_surface_kinetics status runs/<name>/settings.py
    python -m hydrogen_surface_kinetics report runs/<name>/settings.py [--tex-only]

A run of every analysis in settings.py ends with the PDF report (output/report/); `report`
rebuilds it from the saved results without solving anything.
"""
import argparse
import glob
import os
import sys

import matplotlib

matplotlib.use('Agg')                                     # no windows from the command line

from .config import ANALYSIS_ORDER, ConfigError, load_config   # noqa: E402
from .model import Model, ModelError                           # noqa: E402
from .report import ReportError, latexmk_path                  # noqa: E402


def cmd_check(args):
    cfg = load_config(args.settings)
    model = Model(cfg)
    print(model.summary())
    from .conditions import main_conditions, series_points
    from .analyses.common import Context
    from .analyses.plan import counts
    pts = series_points(cfg)
    t_pts = [p for p in pts if p['series'] in ('temperature', 'both')]
    p_pts = [p for p in pts if p['series'] in ('pressure', 'both')]
    temps = ', '.join(f"{p['T']:g}" for p in t_pts)
    print(f'\ntemperature series: {temps} K at {cfg.reference_pressure / 1e5:g} bar')
    if p_pts:
        for T in sorted({p['T'] for p in p_pts}):
            Ps = ', '.join(f"{P / 1e5:g}" for P in sorted(p['P'] for p in p_pts if p['T'] == T))
            print(f'pressure series:    {Ps} bar at {T:g} K')
    else:
        print('pressure series:    none (pressures = [])')
    conds = main_conditions(model)
    print(f'\n{len(conds)} main conditions ({len(pts)} points x cases):')
    print(conds.drop(columns='X').to_string(index=False))
    print(f'\nanalyses: {cfg.analyses}')
    table = counts(Context(model, log=lambda *a: None), cfg.analyses)
    if len(table):
        print('\nsolves per analysis (nothing is solved by check):')
        print(table.to_string(index=False))
        print(f"total to solve: {int(table['to_solve'].sum())}")
    if not cfg.report['enabled']:
        print('\nreport: off (report = {"enabled": False})')
    else:
        exe = latexmk_path(cfg.report)
        print('\nreport: built at the end of a full run, '
              + (f'with {exe}' if exe else 'but LaTeX (latexmk) was not found: only the .tex will be written'))
    print(f'\noutput: {cfg.output_dir}\nsettings OK')


def cmd_run(args):
    from . import analyses
    only = args.only.split(',') if args.only else None
    bad = [n for n in only or [] if n not in ANALYSIS_ORDER]
    if bad:
        raise ConfigError(f'unknown analyses {bad}; choose from {ANALYSIS_ORDER}')
    analyses.run(args.settings, only=only, workers=args.workers, job_timeout=args.timeout,
                 recompute=True if args.recompute else None, report=False if args.no_report else None)


def cmd_report(args):
    from . import report
    report.build(args.settings, latex=not args.tex_only)


def cmd_status(args):
    import pickle
    cfg = load_config(args.settings)
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
    for name, fn, text in (('check', cmd_check, 'load the settings and mechanism, list the conditions'),
                           ('run', cmd_run, 'run the analyses'),
                           ('status', cmd_status, 'what is cached and saved'),
                           ('report', cmd_report, 'build the PDF report from the saved results (no solving)')):
        sp = sub.add_parser(name, help=text)
        sp.add_argument('settings', help="path to the run's settings.py")
        sp.set_defaults(func=fn)
        if name == 'run':
            sp.add_argument('--only', help=f'comma-separated subset of {",".join(ANALYSIS_ORDER)}')
            sp.add_argument('--workers', type=int, help='parallel processes')
            sp.add_argument('--timeout', type=float, help='time limit per job (s)')
            sp.add_argument('--recompute', action='store_true', help='ignore cached results')
            sp.add_argument('--no-report', action='store_true', help='do not build the PDF report at the end')
        if name == 'report':
            sp.add_argument('--tex-only', action='store_true', help='write the .tex but do not run LaTeX')
    args = p.parse_args(argv)
    try:
        args.func(args)
    except (ConfigError, ModelError, ReportError) as e:
        print(f'error: {e}', file=sys.stderr)
        return 2
    return 0
