"""Degree of rate control (DRC) and thermodynamic rate control (TRC) (documentation:
Stage 11d).

  X_DRC,i = d ln r / d ln k_i          (k_f and k_r of reaction i scaled together)
  X_TRC,n = d ln r / d(-G_n / RT)      (species n's enthalpy shifted, transition states
                                        held fixed via rate multipliers)
Every perturbed steady state starts from a clean surface: a warm start from the base
looks converged at once and misses slow re-equilibration.
Runs at every point of the temperature and pressure series.
Outputs: sensitivity.csv; rate_control_vs_T.png and, with a pressure series,
rate_control_vs_P.png (TRC per case and point, plus the reaction orders when the 'orders'
analysis has run).
"""
import os

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from .. import solver as S
from ..conditions import case_labels
from ..labels import probe_label, species_label
from . import common


def _solve(model, job):
    r = job['row']
    theta, info = S.solve_steady_state(model, r['T_K'], r['P_total_bar'] * 1e5, r['X'], S.clean_surface(model),
                                       multipliers=job['multipliers'], h_shifts=job['h_shifts'])
    out = {'theta': theta, 'info': info}
    if theta is not None:      # model.surf holds the solved, perturbed state
        out['tof_fwd'] = model.probe_rates()[0] / model.gamma
        out['probe_consumed'] = -model.surf.net_production_rates[model.probe_kin_idx] / model.gamma
    return out


def drc_reactions(ctx, row, theta):
    """Probe and molecular-adsorption reactions plus the drc_top_n largest |net rates|."""
    model = ctx.model
    surf = common.load_state(model, row, theta)
    order = np.argsort(-np.abs(surf.net_rates_of_progress))
    first = [i for i in (model.i_probe, model.i_molec) if i is not None]
    return list(dict.fromkeys(first + [int(i) for i in order[:ctx.cfg.rate_control['drc_top_n']]]))


def trc_species(ctx, theta):
    model = ctx.model
    roles = [n for n in (model.site, model.probe_adsorbate, model.molecular_adsorbate) if n]
    covered = [n for n, c in zip(model.species_names, theta) if c >= ctx.cfg.rate_control['trc_min_coverage']]
    return list(dict.fromkeys(roles + covered))


def trc_multipliers(model, name, dG, T):
    """Hold transition states fixed: reactions with `name` as a reactant get their rate
    constants scaled by exp(nu' dG / RT)."""
    return {i: float(np.exp(r.reactants[name] * dG / (ct.gas_constant * T)))
            for i, r in enumerate(model.surf.reactions()) if name in r.reactants}


def run(ctx):
    model, cfg, s = ctx.model, ctx.cfg, ctx.cfg.rate_control
    ctx.log('\n== rate_control: DRC and TRC ==')
    ratio = common.setting(ctx, 'rate_control', 'ratio', cfg.impurity_ratio)
    cases = s['cases'] or [cfg.reference_label] + list(cfg.impurities)
    pts = common.points(ctx, 'rate_control')
    func = lambda job: _solve(model, job)                                   # noqa: E731
    note = lambda res: S.describe(res['info'])                              # noqa: E731

    def where(case, pt):
        return f"{case} | {pt['T']:.0f} K | {pt['P'] / 1e5:g} bar"

    rows = {(i, case): common.fixed_probe_row(ctx, 'sensitivity', case, pt, impurity_ratio=ratio)
            for i, pt in enumerate(pts) for case in cases}
    base_jobs = [common.job(f'base | {where(case, pts[i])}', row, multipliers={}, h_shifts={})
                 for (i, case), row in rows.items()]
    base = ctx.run_jobs(func, base_jobs, 'cache_sensitivity', describe=note)

    jobs, plan = [], []
    for (i, case), row in rows.items():
        pt = pts[i]
        b = base[f'base | {where(case, pt)}']
        if b['status'] != 'done' or not b['info']['converged']:
            ctx.log(f'  skipping {where(case, pt)}: base steady state not converged')
            continue
        for k in drc_reactions(ctx, row, b['theta']):
            for sign, m in (('+', 1 + s['drc_delta']), ('-', 1 - s['drc_delta'])):
                jobs.append(common.job(f'drc {k} {sign} | {where(case, pt)}', row, multipliers={k: m}, h_shifts={}))
            plan.append(('DRC', pt, case, k))
        for name in trc_species(ctx, b['theta']):
            for sign, dG in (('+', s['trc_delta']), ('-', -s['trc_delta'])):
                jobs.append(common.job(f'trc {name} {sign} | {where(case, pt)}', row,
                                       multipliers=trc_multipliers(model, name, dG, pt['T']), h_shifts={name: dG}))
            plan.append(('TRC', pt, case, name))
    ctx.log(f'  {len(jobs)} perturbed solves')
    pert = ctx.run_jobs(func, jobs, 'cache_sensitivity', describe=note)

    def log_rate(res, col):
        if res.get('status') != 'done' or not res['info']['converged']:
            return np.nan
        v = res[col]
        return np.log(v) if v > 0 else np.nan

    out = []
    for kind, pt, case, item in plan:
        tag = f'drc {item}' if kind == 'DRC' else f'trc {item}'
        up, dn = pert.get(f'{tag} + | {where(case, pt)}', {}), pert.get(f'{tag} - | {where(case, pt)}', {})
        b = base[f'base | {where(case, pt)}']
        if kind == 'DRC':
            span = np.log(1 + s['drc_delta']) - np.log(1 - s['drc_delta'])
            label = model.surf.reaction(item).equation
        else:
            span = -2 * s['trc_delta'] / (ct.gas_constant * pt['T'])
            label = f'{species_label(model, item)}  [{item}]'
        row = {'kind': kind, 'series': pt['series'], 'T_K': pt['T'], 'P_ref_bar': pt['P'] / 1e5, 'case': case,
               'item': item, 'label': label,
               'methods': f"{up.get('info', {}).get('method')}/{dn.get('info', {}).get('method')}"}
        for col, name in (('tof_fwd', 'X_tof_fwd'), ('probe_consumed', 'X_probe_consumed')):
            row[name] = (log_rate(up, col) - log_rate(dn, col)) / span
        if not (b.get('probe_consumed', 0) > 1e-6 * b.get('tof_fwd', 1)):    # near equilibrium: no meaning
            row['X_probe_consumed'] = np.nan
        out.append(row)
    sens = pd.DataFrame(out)
    ctx.save_table(sens, 'sensitivity.csv')
    if len(sens):
        for (kind, T, Pb, case), part in sens.groupby(['kind', 'T_K', 'P_ref_bar', 'case'], sort=False):
            part = part.reindex(part['X_tof_fwd'].abs().sort_values(ascending=False).index)
            ctx.log(f'  --- {kind}, {case}, {T:.0f} K, {Pb:g} bar '
                    f'(sum of X on forward TOF = {part["X_tof_fwd"].sum():.3f}) ---')
            ctx.log(part[['X_tof_fwd', 'X_probe_consumed', 'methods', 'label']].head(8)
                    .to_string(index=False, float_format=lambda v: f'{v:+.3f}'))
    orders = orders_table(ctx)
    for series in C.SERIES:
        fig = plot_rate_control(ctx, sens, orders, series)
        if fig is not None:
            P.save(fig, ctx.path(common.series_file('rate_control', series)), cfg)
            plt.close(fig)
    ctx.outputs['rate_control'] = {'base': base, 'perturbed': pert, 'sensitivity': sens}
    return ctx.outputs['rate_control']


def orders_table(ctx):
    """Reaction orders from this run, or from an earlier run's reaction_orders.csv."""
    if 'orders' in ctx.outputs:
        return ctx.outputs['orders']['orders']
    path = ctx.path('reaction_orders.csv')
    if os.path.exists(path):
        ctx.log('  reaction orders for the figure: from reaction_orders.csv of an earlier run')
        return pd.read_csv(path)
    return None


def plot_rate_control(ctx, sens, orders=None, series='temperature', min_abs=0.01):
    """For one series. Top: TRC of the forward probe rate, a row of panels (one per case) per
    point, species with |X| >= min_abs. Bottom (if orders are given): the reaction orders in
    the probe gas and the impurity, one panel per point."""
    if orders is not None and 'P_ref_bar' not in orders:
        orders = None                          # a table from before the pressure series: no point to match
    if sens.empty:
        return None
    trc = C.in_series(sens[sens['kind'] == 'TRC'], series)
    if trc.empty:
        return None
    order = case_labels(ctx.cfg)
    pts = list(trc[['T_K', 'P_ref_bar']].drop_duplicates().sort_values(['T_K', 'P_ref_bar']).itertuples(index=False))
    cases = [c for c in order if c in set(trc['case'])]
    ord_pts = [] if orders is None else [p for p in pts
                                         if ((orders['T_K'] == p[0]) & (orders['P_ref_bar'] == p[1])).any()]
    ncol = len(cases)
    mosaic = ([[f't{i}_{j}' for j in range(ncol)] for i in range(len(pts))]
              + [[f'o{i}'] * ncol for i in range(len(ord_pts))])
    fig, axd = plt.subplot_mosaic(mosaic, figsize=(max(3.4 * ncol, 8.0), 2.7 * len(pts) + 3.2 * len(ord_pts)),
                                  height_ratios=[1.0] * len(pts) + [1.15] * len(ord_pts))
    xmin = min(-2.2, trc['X_tof_fwd'].min() * 1.15)
    xmax = max(0.75, trc['X_tof_fwd'].max() * 1.15 + 0.6)
    for i, (T, Pb) in enumerate(pts):
        for j, case in enumerate(cases):
            ax = axd[f't{i}_{j}']
            part = trc[(trc['T_K'] == T) & (trc['P_ref_bar'] == Pb) & (trc['case'] == case)]
            show = part[part['X_tof_fwd'].abs() >= min_abs].sort_values('X_tof_fwd', key=abs, ascending=False)
            y = np.arange(len(show))[::-1]
            ax.barh(y, show['X_tof_fwd'], height=0.6, color=P.SERIES_COLORS[0], zorder=2)
            for yi, x in zip(y, show['X_tof_fwd']):              # values beside the zero line, clear of the names
                ax.annotate(P.signed(x), (max(x, 0), yi), xytext=(4, 0), textcoords='offset points',
                            ha='left', va='center', fontsize=8, color=P.SECONDARY_INK)
            ax.set_yticks(y, [l.split('  [')[0] for l in show['label']])
            ax.set_ylim(-0.6, max(len(show), 2) - 0.4)
            ax.axvline(0, color=P.MUTED, lw=0.8, zorder=1)
            ax.set_xlim(xmin, xmax)
            ax.set_xlabel(f'$X_{{TRC}}$ of the forward rate (sum {P.signed(part["X_tof_fwd"].sum())})',
                          fontsize=9, color=P.INK)
            P.style_axis(ax)
            P.keep_title(ax, f'{case}, {T:.0f} K, {Pb:g} bar')
            ax.grid(axis='y', visible=False)
            ax.grid(axis='x', color=P.GRID, linewidth=0.6)
    probe = probe_label(ctx.model)
    for i, (T, Pb) in enumerate(ord_pts):
        ax = axd[f'o{i}']
        o = orders[(orders['T_K'] == T) & (orders['P_ref_bar'] == Pb)].set_index('case')
        cs = [c for c in order if c in o.index]
        x = np.arange(len(cs))
        for col, dx, color, marker, lab in (('order_probe', -0.12, P.SERIES_COLORS[0], 'o', f'order in {probe}'),
                                           ('order_impurity', 0.12, P.SERIES_COLORS[1], 's', 'order in the impurity')):
            vals = o.loc[cs, col].to_numpy(dtype=float)
            ax.plot(x + dx, vals, ls='', marker=marker, ms=8, color=color, mec='white', mew=0.8, zorder=3,
                    label=lab if i == 0 else None)
            for xi, v in zip(x + dx, vals):
                if np.isfinite(v) and abs(v) >= 0.1:                   # label only the orders that matter
                    ax.annotate(P.signed(v), (xi, v), xytext=(0, 7 if v > 0 else -7), textcoords='offset points',
                                ha='center', va='bottom' if v > 0 else 'top', fontsize=8, color=P.SECONDARY_INK)
        ax.axhline(0, color=P.MUTED, lw=0.8, zorder=1)
        ax.set_xticks(x, cs)
        ax.set_xlim(-0.6, len(cs) - 0.4)
        vals_all = o[['order_probe', 'order_impurity']].to_numpy(dtype=float)
        lo, hi = (np.nanmin(vals_all), np.nanmax(vals_all)) if np.isfinite(vals_all).any() else (0.0, 0.0)
        ax.set_ylim(min(lo, -0.5) - 0.5, max(hi, 0.5) + 0.5)
        P.style_axis(ax, ylabel='reaction order  d ln TOF / d ln p')
        P.keep_title(ax, f'reaction orders, {T:.0f} K, {Pb:g} bar')
    if ord_pts:
        handles, labels = axd['o0'].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(labels),
                   frameon=False, fontsize=9)
    fig.suptitle(f'What controls the rate ({series} series): thermodynamic rate control and reaction orders',
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
