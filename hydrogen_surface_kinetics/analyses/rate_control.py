"""Degree of rate control (DRC) and thermodynamic rate control (TRC) (documentation:
Stage 11d).

  X_DRC,i = d ln r / d ln k_i          (k_f and k_r of reaction i scaled together)
  X_TRC,n = d ln r / d(-G_n / RT)      (species n's enthalpy shifted, transition states
                                        held fixed via rate multipliers)
Every perturbed steady state starts from a clean surface: a warm start from the base
looks converged at once and misses slow re-equilibration.
Runs at every point of the temperature and pressure series.
Outputs: sensitivity.csv; per point, rate_control_<T>K_<P>bar.png (TRC per case) and, where
the DRC is not trivial, rate_control_<T>K_<P>bar_DRC.png (the reaction orders have their own
figure, from the 'orders' analysis).
"""
import re

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
    if len(sens):
        for T, Pb in sens[['T_K', 'P_ref_bar']].drop_duplicates().sort_values(['T_K', 'P_ref_bar']).itertuples(index=False):
            trc_figs, drc_figs = plot_rate_control(ctx, sens, T, Pb)
            stem = f'rate_control_{T:.0f}K_{Pb:g}bar'
            P.save_parts(trc_figs, ctx.path(f'{stem}.png'), cfg)
            P.save_parts(drc_figs, ctx.path(f'{stem}_DRC.png'), cfg)
    ctx.outputs['rate_control'] = {'base': base, 'perturbed': pert, 'sensitivity': sens}
    return ctx.outputs['rate_control']


def drc_nontrivial(part, tol=0.02):
    """Is this (point, case)'s DRC worth a panel? Yes if the net probe-gas consumption has a
    DRC (it is far from equilibrium), or if the forward-rate DRC is not simply 1 on one step."""
    if part.empty:
        return False
    if part['X_probe_consumed'].notna().any():
        return True
    x = part['X_tof_fwd'].abs().dropna().sort_values(ascending=False).to_numpy()
    return len(x) > 0 and not (abs(x[0] - 1) < tol and (len(x) == 1 or x[1] < tol))


def short_equation(model, eq, width=42):
    """A reaction equation with readable species names, shortened to `width` characters."""
    surface = set(model.species_names)

    def name(term):
        n, _, sp = term.partition(' ') if term[:1].isdigit() and ' ' in term else ('', '', term)
        label = species_label(model, sp) if sp in surface else re.sub(r'\(\d+\)$', '', sp)
        return f'{n} {label}' if n else label
    sides = [' + '.join(name(t.strip()) for t in side.split(' + ')) for side in re.split(r'\s*<=>\s*|\s*=>\s*', eq)]
    text = ' <=> '.join(sides)
    return text if len(text) <= width else text[:width - 1] + '…'


def plot_rate_control(ctx, sens, T, Pb, min_abs=0.01):
    """One point (T, P), as (TRC figures, DRC figures). TRC of the forward probe rate: one panel
    per case, species with |X| >= min_abs, at most 3 x 3 per figure. DRC of the forward rate
    and of the net probe-gas consumption: one panel per case where it is not trivial (see
    drc_nontrivial), at most 3 x 2 per figure (the reaction labels are long)."""
    at = sens[(sens['T_K'] == T) & (sens['P_ref_bar'] == Pb)]
    trc, drc = at[at['kind'] == 'TRC'], at[at['kind'] == 'DRC']
    if trc.empty:
        return [], []
    cases = [c for c in case_labels(ctx.cfg) if c in set(trc['case'])]
    drc_cases = [c for c in cases if drc_nontrivial(drc[drc['case'] == c])]
    where = f'{T:.0f} K, {Pb:g} bar'
    xlim = (min(-2.2, trc['X_tof_fwd'].min() * 1.15), max(0.75, trc['X_tof_fwd'].max() * 1.15 + 0.6))
    return ([_trc_page(ctx, trc, page, where, xlim, min_abs) for page in P.pages(cases, 3 * P.MAX_ROWS)],
            [_drc_page(ctx, drc, page, where, min_abs) for page in P.pages(drc_cases, 2 * P.MAX_ROWS)])


def _trc_page(ctx, trc, cases, where, xlim, min_abs):
    nrow, ncol = P.grid_shape(len(cases), max_cols=3)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.2 * ncol, 2.8 * nrow + 0.6), squeeze=False)
    for ax, case in zip(axes.flat, cases):
        part = trc[trc['case'] == case]
        show = part[part['X_tof_fwd'].abs() >= min_abs].sort_values('X_tof_fwd', key=abs, ascending=False)
        y = np.arange(len(show))[::-1]
        ax.barh(y, show['X_tof_fwd'], height=0.6, color=P.SERIES_COLORS[0], zorder=2)
        for yi, x in zip(y, show['X_tof_fwd']):                  # values beside the zero line, clear of the names
            ax.annotate(P.signed(x), (max(x, 0), yi), xytext=(4, 0), textcoords='offset points',
                        ha='left', va='center', fontsize=8, color=P.SECONDARY_INK)
        ax.set_yticks(y, [l.split('  [')[0] for l in show['label']])
        ax.set_ylim(-0.6, max(len(show), 2) - 0.4)
        ax.axvline(0, color=P.MUTED, lw=0.8, zorder=1)
        ax.set_xlim(*xlim)
        ax.set_xlabel(f'$X_{{TRC}}$ (sum {P.signed(part["X_tof_fwd"].sum())})', fontsize=9, color=P.INK)
        P.style_axis(ax)
        P.keep_title(ax, f'{case}, {where}')
        ax.grid(axis='y', visible=False)
        ax.grid(axis='x', color=P.GRID, linewidth=0.6)
    for ax in axes.flat[len(cases):]:
        ax.set_visible(False)
    fig.suptitle(f'Thermodynamic rate control of the forward {probe_label(ctx.model)} rate, {where}',
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig


def _drc_page(ctx, drc, cases, where, min_abs):
    probe = probe_label(ctx.model)
    nrow, ncol = P.grid_shape(len(cases), max_cols=2)
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.6 * ncol, 2.9 * nrow + 0.8), squeeze=False)
    for ax, case in zip(axes.flat, cases):
        part = drc[drc['case'] == case]
        size = part[['X_tof_fwd', 'X_probe_consumed']].abs().max(axis=1)
        show = part[size >= min_abs].assign(_size=size).sort_values('_size', ascending=False).head(6)
        y = np.arange(len(show))[::-1]
        for col, dy, color in (('X_tof_fwd', 0.19, P.SERIES_COLORS[0]), ('X_probe_consumed', -0.19, P.SERIES_COLORS[1])):
            ax.barh(y + dy, np.nan_to_num(show[col].to_numpy(dtype=float)), height=0.36, color=color, zorder=2)
        ax.set_yticks(y, [f"#{int(r.item)}: {short_equation(ctx.model, r.label)}" for r in show.itertuples()],
                      fontsize=7)
        ax.set_ylim(-0.6, max(len(show), 2) - 0.4)
        ax.axvline(0, color=P.MUTED, lw=0.8, zorder=1)
        net_sum = part['X_probe_consumed'].sum(min_count=1)
        ax.set_xlabel('$X_{DRC}$' + (f' (net sum {P.signed(net_sum)})' if np.isfinite(net_sum) else ''),
                      fontsize=9, color=P.INK)
        P.style_axis(ax)
        P.keep_title(ax, f'{case}, {where}: DRC')
        ax.grid(axis='y', visible=False)
        ax.grid(axis='x', color=P.GRID, linewidth=0.6)
    for ax in axes.flat[len(cases):]:
        ax.set_visible(False)
    handles = [plt.Rectangle((0, 0), 1, 1, color=P.SERIES_COLORS[0]), plt.Rectangle((0, 0), 1, 1, color=P.SERIES_COLORS[1])]
    fig.legend(handles, ['DRC of the forward rate', f'DRC of net {probe} consumption'], loc='upper center',
               bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig
