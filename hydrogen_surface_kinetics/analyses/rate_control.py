"""Degree of rate control (DRC) and thermodynamic rate control (TRC) (documentation:
Stage 11d).

  X_DRC,i = d ln r / d ln k_i          (k_f and k_r of reaction i scaled together)
  X_TRC,n = d ln r / d(-G_n / RT)      (species n's enthalpy shifted, transition states
                                        held fixed via rate multipliers)
Every perturbed steady state starts from a clean surface: a warm start from the base
looks converged at once and misses slow re-equilibration.
Outputs: sensitivity.csv; rate_control.png (TRC per case, plus the reaction orders when the
'orders' analysis has run).
"""
import os

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

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
    temps = common.temperatures(ctx, 'rate_control')
    func = lambda job: _solve(model, job)                                   # noqa: E731
    note = lambda res: S.describe(res['info'])                              # noqa: E731

    rows = {(T, case): common.fixed_probe_row(ctx, 'sensitivity', case, T, impurity_ratio=ratio) for T in temps for case in cases}
    base_jobs = [{'key': f'base | {case} | {T:.0f} K', 'row': row, 'multipliers': {}, 'h_shifts': {}}
                 for (T, case), row in rows.items()]
    base = ctx.run_jobs(func, base_jobs, 'cache_sensitivity', describe=note)

    jobs, plan = [], []
    for (T, case), row in rows.items():
        b = base[f'base | {case} | {T:.0f} K']
        if b['status'] != 'done' or not b['info']['converged']:
            ctx.log(f'  skipping {case} {T:.0f} K: base steady state not converged')
            continue
        for i in drc_reactions(ctx, row, b['theta']):
            for sign, m in (('+', 1 + s['drc_delta']), ('-', 1 - s['drc_delta'])):
                jobs.append({'key': f'drc {i} {sign} | {case} | {T:.0f} K', 'row': row,
                             'multipliers': {i: m}, 'h_shifts': {}})
            plan.append(('DRC', T, case, i))
        for name in trc_species(ctx, b['theta']):
            for sign, dG in (('+', s['trc_delta']), ('-', -s['trc_delta'])):
                jobs.append({'key': f'trc {name} {sign} | {case} | {T:.0f} K', 'row': row,
                             'multipliers': trc_multipliers(model, name, dG, T), 'h_shifts': {name: dG}})
            plan.append(('TRC', T, case, name))
    ctx.log(f'  {len(jobs)} perturbed solves')
    pert = ctx.run_jobs(func, jobs, 'cache_sensitivity', describe=note)

    def log_rate(res, col):
        if res.get('status') != 'done' or not res['info']['converged']:
            return np.nan
        v = res[col]
        return np.log(v) if v > 0 else np.nan

    out = []
    for kind, T, case, item in plan:
        tag = f'drc {item}' if kind == 'DRC' else f'trc {item}'
        up, dn = pert.get(f'{tag} + | {case} | {T:.0f} K', {}), pert.get(f'{tag} - | {case} | {T:.0f} K', {})
        b = base[f'base | {case} | {T:.0f} K']
        if kind == 'DRC':
            span = np.log(1 + s['drc_delta']) - np.log(1 - s['drc_delta'])
            label = model.surf.reaction(item).equation
        else:
            span = -2 * s['trc_delta'] / (ct.gas_constant * T)
            label = f'{species_label(model, item)}  [{item}]'
        row = {'kind': kind, 'T_K': T, 'case': case, 'item': item, 'label': label,
               'methods': f"{up.get('info', {}).get('method')}/{dn.get('info', {}).get('method')}"}
        for col, name in (('tof_fwd', 'X_tof_fwd'), ('probe_consumed', 'X_probe_consumed')):
            row[name] = (log_rate(up, col) - log_rate(dn, col)) / span
        if not (b.get('probe_consumed', 0) > 1e-6 * b.get('tof_fwd', 1)):    # near equilibrium: no meaning
            row['X_probe_consumed'] = np.nan
        out.append(row)
    sens = pd.DataFrame(out)
    ctx.save_table(sens, 'sensitivity.csv')
    for (kind, T, case), part in sens.groupby(['kind', 'T_K', 'case'], sort=False):
        part = part.reindex(part['X_tof_fwd'].abs().sort_values(ascending=False).index)
        ctx.log(f'  --- {kind}, {case}, {T:.0f} K (sum of X on forward TOF = {part["X_tof_fwd"].sum():.3f}) ---')
        ctx.log(part[['X_tof_fwd', 'X_probe_consumed', 'methods', 'label']].head(8)
                .to_string(index=False, float_format=lambda v: f'{v:+.3f}'))
    orders = orders_table(ctx)
    fig = plot_rate_control(ctx, sens, orders)
    if fig is not None:
        P.save(fig, ctx.path('rate_control.png'), cfg)
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


def plot_rate_control(ctx, sens, orders=None, min_abs=0.01):
    """Top: TRC of the forward probe rate, one panel per (T, case), species with |X| >= min_abs.
    Bottom (if orders are given): reaction orders in the probe gas and the impurity, per T."""
    trc = sens[sens['kind'] == 'TRC']
    if trc.empty:
        return None
    order = case_labels(ctx.cfg)
    temps = sorted(trc['T_K'].unique())
    cases = [c for c in order if c in set(trc['case'])]
    ord_T = [] if orders is None else [T for T in temps if T in set(orders['T_K'])]
    ncol = len(cases)
    mosaic = ([[f't{i}_{j}' for j in range(ncol)] for i in range(len(temps))]
              + [[f'o{i}'] * ncol for i in range(len(ord_T))])
    fig, axd = plt.subplot_mosaic(mosaic, figsize=(max(3.4 * ncol, 8.0), 2.7 * len(temps) + 3.2 * len(ord_T)),
                                  height_ratios=[1.0] * len(temps) + [1.15] * len(ord_T))
    xmin = min(-2.2, trc['X_tof_fwd'].min() * 1.15)
    xmax = max(0.75, trc['X_tof_fwd'].max() * 1.15 + 0.6)
    for i, T in enumerate(temps):
        for j, case in enumerate(cases):
            ax = axd[f't{i}_{j}']
            part = trc[(trc['T_K'] == T) & (trc['case'] == case)]
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
            P.keep_title(ax, f'{case}, {T:.0f} K')
            ax.grid(axis='y', visible=False)
            ax.grid(axis='x', color=P.GRID, linewidth=0.6)
    probe = probe_label(ctx.model)
    for i, T in enumerate(ord_T):
        ax = axd[f'o{i}']
        o = orders[orders['T_K'] == T].set_index('case')
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
        lo = np.nanmin(o[['order_probe', 'order_impurity']].to_numpy(dtype=float))
        hi = np.nanmax(o[['order_probe', 'order_impurity']].to_numpy(dtype=float))
        ax.set_ylim(min(lo, -0.5) - 0.5, max(hi, 0.5) + 0.5)
        P.style_axis(ax, ylabel='reaction order  d ln TOF / d ln p')
        P.keep_title(ax, f'reaction orders, {T:.0f} K')
    if ord_T:
        handles, labels = axd['o0'].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(labels),
                   frameon=False, fontsize=9)
    fig.suptitle('What controls the rate: thermodynamic rate control and reaction orders', x=0.01, ha='left',
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
