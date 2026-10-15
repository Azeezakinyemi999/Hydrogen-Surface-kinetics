"""How many steady states (or transients) each chosen analysis needs, and how many are
already cached, without solving anything. Shown by `check`.

Exact for the analyses built on condition tables (main, curves, sweeps, orders). For
rate_control and transients the job list depends on solved states (the screened reactions,
the starting surface), so their counts are estimates and their cache check goes by name.
"""
import glob
import os

import pandas as pd

from .. import conditions as C
from .. import runner as R
from .. import solver as S
from . import common


def _exact(ctx, conditions, cache='cache'):
    cache_dir = os.path.join(ctx.cfg.output_dir, cache)
    context = (ctx.model.fingerprint_context(), True)
    jobs = S.condition_jobs(conditions)
    cached = sum(os.path.exists(R.cache_path(cache_dir, j, context)) for j in jobs)
    return len(jobs), cached


def _by_name(ctx, cache, keys):
    """Jobs whose name already has a cache file (any fingerprint: an upper bound)."""
    folder = os.path.join(ctx.cfg.output_dir, cache)
    return sum(bool(glob.glob(os.path.join(folder, R.safe_name(k) + '_' + '?' * 12 + '.pkl'))) for k in keys)


def _where(case, pt):
    return f"{case} | {pt['T']:.0f} K | {pt['P'] / 1e5:g} bar"


def counts(ctx, names):
    """DataFrame: analysis, points, solves, cached, to_solve, note."""
    from . import curves, orders, pressure_sweep, sweep
    cfg = ctx.cfg
    rows = []

    def add(name, n_points, solves, cached, note=''):
        rows.append({'analysis': name, 'points': n_points, 'solves': solves, 'cached': cached,
                     'to_solve': max(solves - cached, 0), 'note': note})

    main_pts = C.series_points(cfg)
    if 'steady_state' in names or any(n in names for n in ('composition', 'diagnostics', 'paths', 'transients')):
        add('steady_state', len(main_pts), *_exact(ctx, C.main_conditions(ctx.model)), 'two starts each')
    if 'diagnostics' in names:
        s = cfg.diagnostics
        species = s['switch_off'] or ctx.model.molecular_adsorbate
        if species:
            pts = common.points(ctx, 'diagnostics')
            cases = s['switch_off_cases'] or C.case_labels(cfg)
            keys = [f"switch_off {species} | {_where(c, pt)}" for pt in pts for c in cases]
            add('diagnostics', len(pts), len(keys), _by_name(ctx, 'cache_diagnostics', keys), f'{species} switched off')
    if 'curves' in names:
        add('curves', len(common.points(ctx, 'curves')), *_exact(ctx, curves.curve_conditions(ctx)), 'two starts each')
    if 'sweep' in names:
        add('sweep', len(sweep.sweep_temperatures(ctx)), *_exact(ctx, sweep.sweep_conditions(ctx)),
            'temperature only, at reference_pressure')
    if 'arrhenius' in names:
        add('arrhenius', 0, 0, 0, 'no solving: rate constants from the mechanism, fits of the temperature sweep')
    if 'pressure_sweep' in names:
        n_pts = len(pressure_sweep.sweep_pressures(ctx)) * len(pressure_sweep.sweep_temperatures(ctx))
        add('pressure_sweep', n_pts, *_exact(ctx, pressure_sweep.sweep_conditions(ctx)),
            'pressure only, at pressure_temperatures')
    if 'equilibrium' in names:
        add('equilibrium', 0, 0, 0, 'no solving: equilibrium coverage from the mechanism thermochemistry')
    if 'convergence' in names:
        conds = cfg.convergence['conditions']
        n = len(conds) if conds else f"<= {cfg.convergence['max_conditions']}"
        add('convergence', 0, 0, 0, f'{n} traces to reactor_t_max (failed main conditions by default)')
    if 'orders' in names:
        add('orders', len(common.points(ctx, 'orders')), *_exact(ctx, orders.order_conditions(ctx)), 'two starts each')
    if 'rate_control' in names:
        s = cfg.rate_control
        pts = common.points(ctx, 'rate_control')
        cases = s['cases'] or [cfg.reference_label] + list(cfg.impurities)
        n_drc = s['drc_top_n'] + (2 if ctx.model.i_molec is not None else 1)
        n_trc = 4                                          # the role species plus about one covered species
        per_base = 1 + 2 * n_drc + 2 * n_trc
        keys = [f'base | {_where(c, pt)}' for pt in pts for c in cases]
        add('rate_control', len(pts), len(keys) * per_base, _by_name(ctx, 'cache_sensitivity', keys) * per_base,
            f'about {per_base} solves per case and point (one start each); estimate')
    if 'transients' in names:
        s = cfg.transients
        pts = common.points(ctx, 'transients')
        imps = common.impurity_labels(ctx, 'transients')
        keys = [f"transient r={r:.3g} | {imp} | {pt['T']:.0f} K | {pt['P'] / 1e5:g} bar"
                for pt in pts for imp in imps for r in s['ratios']]
        add('transients', len(pts), len(keys), _by_name(ctx, 'cache_transients', keys),
            'time integrations; cache checked by name')
    if 'paths' in names:
        add('paths', len(common.points(ctx, 'paths')), 0, 0, 'no solving: diagrams of main steady states')
    if 'flux' in names:
        add('flux', 0, 0, 0, 'one batch reactor per RMG surfaceReactor; videos skipped if present')
    return pd.DataFrame(rows)
