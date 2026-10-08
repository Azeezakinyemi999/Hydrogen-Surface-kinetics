"""Degree of rate control (DRC) and thermodynamic rate control (TRC) (documentation:
Stage 11d).

  X_DRC,i = d ln r / d ln k_i          (k_f and k_r of reaction i scaled together)
  X_TRC,n = d ln r / d(-G_n / RT)      (species n's enthalpy shifted, transition states
                                        held fixed via rate multipliers)
Every perturbed steady state starts from a clean surface: a warm start from the base
looks converged at once and misses slow re-equilibration.
Output: sensitivity.csv.
"""
import cantera as ct
import numpy as np
import pandas as pd

from .. import solver as S
from ..labels import species_label
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
    ctx.outputs['rate_control'] = {'base': base, 'perturbed': pert, 'sensitivity': sens}
    return ctx.outputs['rate_control']
