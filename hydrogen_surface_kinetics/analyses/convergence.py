"""Convergence traces of the reactor method (documentation: Stage 4, Figure 4 on real data).

For chosen conditions, the fixed-gas surface is integrated from a clean surface exactly as
the reactor method does (same reactor, tolerances and checkpoints), but every checkpoint is
recorded and the integration continues to reactor_t_max. It shows how the coverages and the
forward TOF settle, how the normalized residual max_k |dtheta_k/dt| / allowed_k falls, and at
which checkpoint the reactor method accepts the state: the first one where the residual test
passes and TOF and coverages stopped changing since the previous checkpoint. For a
condition that never converged it shows why (for example a slow conversion still going on).

Default conditions: the main conditions whose clean start did not converge; if all
converged, the reference case at the first series point.
Outputs: convergence.csv (one row per checkpoint), convergence.png.
"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import plotting as P
from .. import solver as S
from ..labels import species_label
from . import common


def checkpoints(ctx):
    s, t_max = ctx.cfg.convergence, ctx.cfg.solver['reactor_t_max']
    n = int(round(s['per_decade'] * (np.log10(t_max) + 8))) + 1
    return np.logspace(-8, np.log10(t_max), n)


def trace(model, job):
    """Integrate from a clean surface, recording every checkpoint (the reactor method's test)."""
    s, r = model.cfg.solver, job['row']
    gas_X_start = S.set_state(model, r['T_K'], r['P_total_bar'] * 1e5, r['X'], S.clean_surface(model))
    net = S.fixed_gas_reactor(model)
    model.gas.set_multiplier(0.0)
    rows, thetas, previous, accepted, err = [], [], None, None, None
    try:
        for t in job['times']:
            net.advance(t)
            d = S.check_convergence(model, gas_X_start)
            theta = model.surf.coverages.copy()
            plateau = False
            if previous is not None:
                tof_change = abs(d['tof_fwd_1_per_s'] - previous[0]) / max(abs(previous[0]), 1e-300)
                plateau = bool(tof_change < s['plateau_tof_rtol']
                               and np.abs(theta - previous[1]).max() < s['plateau_theta_atol'])
            ok = bool(d['residual_ok'] and plateau)
            if ok and accepted is None:
                accepted = float(t)
            rows.append({'t_s': float(t), 'tof_fwd': d['tof_fwd_1_per_s'],
                         'norm_residual': d['residual_1_per_s'] / d['residual_allowed_1_per_s'],
                         'worst_species': d['residual_species'], 'residual_ok': bool(d['residual_ok']),
                         'plateau_ok': plateau, 'accepted_here': ok and accepted == float(t)})
            thetas.append(theta)
            previous = (d['tof_fwd_1_per_s'], theta)
    except ct.CanteraError as e:
        err = S._cantera_error(e)
    finally:
        model.gas.set_multiplier(1.0)
    return {'trajectory': pd.DataFrame(rows), 'coverages': np.array(thetas), 'accepted_s': accepted, 'error': err}


def default_conditions(ctx):
    summary = ctx.main()['summary']
    main = summary[summary['basis'] == 'fixed_probe']
    failed = main[~main['converged'].astype(bool)]
    if len(failed):
        return [(r.case, r.T_K, r.P_ref_bar * 1e5) for r in failed.itertuples()][:ctx.cfg.convergence['max_conditions']]
    pt = common.points(ctx)[0]
    return [(ctx.cfg.reference_label, pt['T'], pt['P'])]


def run(ctx):
    model, cfg = ctx.model, ctx.cfg
    ctx.log('\n== convergence: reactor-method traces ==')
    main = ctx.main()
    summary = main['summary']
    conds = cfg.convergence['conditions'] or default_conditions(ctx)
    times = checkpoints(ctx)
    jobs = []
    for cond in conds:
        case, T = cond[0], float(cond[1])
        Pa = float(cond[2]) if len(cond) > 2 else cfg.reference_pressure
        hit = summary[(summary['basis'] == 'fixed_probe') & (summary['case'] == case) & (summary['T_K'] == T)
                      & np.isclose(summary['P_ref_bar'].astype(float), Pa / 1e5)]
        if hit.empty:
            ctx.log(f'  {case}, {T:.0f} K, {Pa / 1e5:g} bar: not a main condition, skipped')
            continue
        res = main['results'][hit['key'].iloc[0]]
        jobs.append(common.job(f"convergence | {case} | {T:.0f} K | {Pa / 1e5:g} bar", res['row'], times=times))
    ctx.log(f'  {len(jobs)} traces')
    out = ctx.run_jobs(lambda job: trace(model, job), jobs, 'cache_convergence',
                       describe=lambda r: f"accepted at {r['accepted_s']:.3g} s" if r['accepted_s'] else
                       (r['error'] or 'never accepted'))
    tables = []
    for job in jobs:
        res = out[job['key']]
        if res['status'] != 'done':
            continue
        r = job['row']
        tables.append(res['trajectory'].assign(case=r['case'], T_K=r['T_K'], P_ref_bar=job['_labels']['P_ref_bar'],
                                               accepted_s=res['accepted_s'], error=res['error']))
    if tables:
        ctx.save_table(pd.concat(tables), 'convergence.csv')
        fig = plot_convergence(ctx, jobs, out)
        P.save(fig, ctx.path('convergence.png'), cfg)
        plt.close(fig)
    ctx.outputs['convergence'] = out
    return out


def plot_convergence(ctx, jobs, results, max_named=6, min_peak=0.01):
    """One row per condition: (a) coverages, (b) forward TOF, (c) normalized residual, against
    time; a dashed line where the reactor method accepts the state."""
    model = ctx.model
    jobs = [j for j in jobs if results[j['key']]['status'] == 'done' and len(results[j['key']]['trajectory'])]
    fig, axes = plt.subplots(len(jobs), 3, figsize=(15, 3.9 * len(jobs)), squeeze=False)
    for row, job in zip(axes, jobs):
        res, r = results[job['key']], job['row']
        traj, theta = res['trajectory'], res['coverages']
        t = traj['t_s'].to_numpy()
        peak = theta.max(axis=0)
        order = [k for k in np.argsort(-peak) if peak[k] >= min_peak][:max_named]
        for i, k in enumerate(order):
            name = model.species_names[k]
            color = P.EMPTY_COLOR if name == model.site else P.SERIES_COLORS[i % len(P.SERIES_COLORS)]
            row[0].plot(t, theta[:, k], color=P.MUTED if name == model.site else color, lw=1.6,
                        ls='--' if name == model.site else '-', label=species_label(model, name))
        row[0].set_ylim(0, 1.02)
        row[0].legend(loc='center left', frameon=False, fontsize=8)
        row[1].plot(t, traj['tof_fwd'], color=P.SERIES_COLORS[0], lw=1.6)
        row[1].set_yscale('log')
        row[2].plot(t, traj['norm_residual'], color=P.SERIES_COLORS[0], lw=1.6)
        row[2].axhline(1.0, color=P.MUTED, lw=1, ls='--')
        row[2].text(t[0], 1.0, ' residual limit', color=P.MUTED, fontsize=8, va='bottom')
        row[2].set_yscale('log')
        name = f"{r['case']}, {r['T_K']:.0f} K, {job['_labels']['P_ref_bar']:g} bar"
        for ax, desc, ylabel in zip(row, ('coverages', 'forward TOF', 'normalized residual'),
                                    ('fraction of sites', 'TOF (1/s)', 'max |dθ/dt| / allowed')):
            ax.set_xscale('log')
            ax.set_xlabel('surface time (s)')
            if res['accepted_s']:
                ax.axvline(res['accepted_s'], color=P.INK, lw=1.0, ls=':')
            P.style_axis(ax, ylabel=ylabel)
            P.keep_title(ax, f'{name}: {desc}')
        note = (f"accepted at {res['accepted_s']:.3g} s" if res['accepted_s'] else
                (f"stopped: {res['error'][:60]}" if res['error'] else 'never accepted (no plateau)'))
        row[2].text(0.98, 0.95, note, transform=row[2].transAxes, ha='right', va='top', fontsize=8,
                    color=P.SECONDARY_INK)
    fig.suptitle('Reactor-method convergence from a clean surface', x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig
