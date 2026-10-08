"""Poisoning over time (documentation: Stage 11f).

Start from the steady state in the probe gas alone (a working catalyst), switch the gas to
probe + impurity at t = 0, and follow the forward rate as the impurity builds up.
Outputs: poisoning_over_time.csv, poisoning_over_time_trajectories.csv, poisoning_over_time.png.
"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import plotting as P
from .. import solver as S
from ..labels import probe_label
from . import common


def times(ctx):
    s = ctx.cfg.transients
    n = int(round(s['per_decade'] * (np.log10(s['t_stop']) - np.log10(s['t_start'])))) + 1
    return np.logspace(np.log10(s['t_start']), np.log10(s['t_stop']), n)


def run_transient(model, job):
    """Fixed-gas surface integration (the Method 3 reactor), recording the state at job['times']."""
    r = job['row']
    gas_X_start = S.set_state(model, r['T_K'], r['P_total_bar'] * 1e5, r['X'], job['theta0'])
    net = S.fixed_gas_reactor(model)
    model.gas.set_multiplier(0.0)
    G = model.gamma
    tof0 = model.probe_rates()[0] / G                    # at t = 0, before any impurity adsorbs
    names = model.species_names
    i_site, i_probe = names.index(model.site), names.index(model.probe_adsorbate)
    i_molec = names.index(model.molecular_adsorbate) if model.molecular_adsorbate else None
    rows, err = [], None
    try:
        for t in job['times']:
            net.advance(t)
            th = model.surf.coverages
            molec = th[i_molec] if i_molec is not None else 0.0
            rows.append({'t_s': t, 'tof_fwd': model.probe_rates()[0] / G, 'theta_free': th[i_site],
                         'theta_probe': th[i_probe], 'theta_molec': molec,
                         'theta_impurity': max(0.0, 1 - th[i_site] - th[i_probe] - molec),
                         'top_species': common.top_species(model, th)})
    except ct.CanteraError as e:
        err = S._cantera_error(e)
    finally:
        model.gas.set_multiplier(1.0)
    return {'trajectory': pd.DataFrame(rows), 'tof0': tof0, 'error': err,
            'gas_drift': np.abs(model.gas.X - gas_X_start).max()}


def crossing_time(traj, tof0, frac):
    """First time the forward TOF falls to `frac` of its t = 0 value (log-interpolated)."""
    t, y = traj['t_s'].to_numpy(), traj['tof_fwd'].to_numpy() / tof0
    j = np.nonzero(y <= frac)[0]
    if len(j) == 0:
        return np.nan
    j = j[0]
    if j == 0:
        return t[0]
    w = (np.log10(frac) - np.log10(y[j - 1])) / (np.log10(y[j]) - np.log10(y[j - 1]))
    return 10 ** (np.log10(t[j - 1]) + w * (np.log10(t[j]) - np.log10(t[j - 1])))


def settle_time(traj, rtol=0.1):
    """First time after which the TOF stays within rtol of its final value."""
    y = traj['tof_fwd'].to_numpy()
    off = np.abs(y / y[-1] - 1) > rtol
    if not off.any():
        return traj['t_s'].iloc[0]
    last = np.nonzero(off)[0][-1]
    return traj['t_s'].iloc[last + 1] if last + 1 < len(traj) else np.nan


def reference_start(ctx, T):
    """Coverages of the probe-only steady state at T: from the main run if it has T."""
    cfg = ctx.cfg
    if float(T) in cfg.temperatures:
        summary = ctx.main()['summary']
        hit = summary[(summary['basis'] == 'fixed_probe') & (summary['case'] == cfg.reference_label)
                      & (summary['T_K'] == float(T))]
        res = ctx.main()['results'][hit['key'].iloc[0]]
    else:
        conds = pd.DataFrame([common.fixed_probe_row(ctx, 'fixed_probe', cfg.reference_label, T)])
        res = next(iter(ctx.run_conditions(conds, prefix='').values()))
    return S.coverages_to_dict(ctx.model, res['theta_a']) if common.converged(res) else None


def run(ctx):
    model, cfg, s = ctx.model, ctx.cfg, ctx.cfg.transients
    ctx.log('\n== transients: poisoning over time ==')
    tt = times(ctx)
    imps = common.impurity_labels(ctx, 'transients')
    jobs = []
    for T in common.temperatures(ctx, 'transients'):
        start = reference_start(ctx, T)
        if start is None:
            ctx.log(f'  skipping {T:.0f} K: no converged {cfg.reference_label} steady state')
            continue
        for imp in imps:
            for ratio in s['ratios']:
                row = common.fixed_probe_row(ctx, 'transient', imp, T, impurity_ratio=ratio, ratio=ratio)
                jobs.append({'key': f'transient r={ratio:.3g} | {imp} | {T:.0f} K', 'row': row,
                             'theta0': start, 'times': tt})
    ctx.log(f'  {len(jobs)} transient runs')
    res_all = ctx.run_jobs(lambda job: run_transient(model, job), jobs, 'cache_transients',
                           describe=lambda res: res['error'] or '')
    rows, trajs = [], []
    for job in jobs:
        res, r = res_all[job['key']], job['row']
        out = {'T_K': r['T_K'], 'case': r['case'], 'ratio': r['ratio'],
               'error': res.get('error') if res['status'] == 'done' else res['status']}
        traj = res.get('trajectory')
        if traj is not None and len(traj):
            tof0 = res['tof0']
            traj = traj.assign(T_K=r['T_K'], case=r['case'], ratio=r['ratio'])
            res['trajectory'] = traj
            j = int(traj['tof_fwd'].to_numpy().argmin())
            out.update({'tof_t0': tof0, 'tof_end': traj['tof_fwd'].iloc[-1],
                        'end_over_t0': traj['tof_fwd'].iloc[-1] / tof0,
                        'min_over_t0': traj['tof_fwd'].iloc[j] / tof0, 't_min_s': traj['t_s'].iloc[j],
                        't_end_s': traj['t_s'].iloc[-1],
                        't_half_s': crossing_time(traj, tof0, 0.5), 't_tenth_s': crossing_time(traj, tof0, 0.1),
                        't_settle_s': settle_time(traj), 'gas_drift': res['gas_drift'],
                        'top_species_end': traj['top_species'].iloc[-1]})
            trajs.append(traj)
        rows.append(out)
    summary = pd.DataFrame(rows)
    ctx.save_table(summary, 'poisoning_over_time.csv')
    if trajs:
        ctx.save_table(pd.concat(trajs), 'poisoning_over_time_trajectories.csv')
    ctx.log(summary.drop(columns=['top_species_end'], errors='ignore')
            .to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    if trajs:
        fig = plot_poisoning_over_time(ctx, jobs, res_all)
        P.save(fig, ctx.path('poisoning_over_time.png'), cfg)
        plt.close(fig)
    ctx.outputs['transients'] = {'results': res_all, 'summary': summary}
    return ctx.outputs['transients']


def plot_poisoning_over_time(ctx, jobs, results):
    s = ctx.cfg.transients
    probe = probe_label(ctx.model)
    shades = ['#9ec5f4', '#5598e7', '#256abf', '#0d366b', '#061b36']        # light -> dark = more impurity
    ratiocolors = dict(zip(sorted(s['ratios']), shades if len(s['ratios']) <= len(shades) else P.SERIES_COLORS))
    temps = sorted({j['row']['T_K'] for j in jobs})
    imps = [i for i in common.impurity_labels(ctx, 'transients') if any(j['row']['case'] == i for j in jobs)]
    fig, axes = plt.subplots(2 * len(temps), len(imps), figsize=(3.4 * len(imps), 5.6 * len(temps)),
                             squeeze=False, sharex=True)
    for ti, T in enumerate(temps):
        for ci, imp in enumerate(imps):
            ax_r, ax_c = axes[2 * ti][ci], axes[2 * ti + 1][ci]
            for job in jobs:
                r = job['row']
                if r['T_K'] != T or r['case'] != imp:
                    continue
                res = results[job['key']]
                traj = res.get('trajectory')
                if traj is None or not len(traj):
                    continue
                c = ratiocolors.get(r['ratio'], P.INK)
                ax_r.plot(traj['t_s'], traj['tof_fwd'] / res['tof0'], color=c, lw=1.8, label=f"{r['ratio']:g}")
                ax_c.plot(traj['t_s'], P.log_floor(traj['theta_impurity']), color=c, lw=1.8)
            ax_r.axhline(0.5, color=P.MUTED, lw=0.8, ls='--')
            for ax in (ax_r, ax_c):
                ax.set_xscale('log')
                ax.set_yscale('log')
            lo, hi = ax_r.get_ylim()
            ax_r.set_ylim(min(lo, 0.1), max(hi, 2.0))      # at least a decade: whole-decade ticks only
            P.style_axis(ax_r, f'{imp}, {T:.0f} K', 'TOF / TOF(t=0)' if ci == 0 else '')
            P.style_axis(ax_c, '', 'sites taken by impurity' if ci == 0 else '')
    fig.supxlabel('time after the impurity arrives (s)', fontsize=10, color=P.INK)
    handles, labels = axes[0][0].get_legend_handles_labels()
    fig.legend(handles, labels, title=f'p(impurity) / p({probe})', loc='upper center', bbox_to_anchor=(0.5, 1.03),
               ncol=max(1, len(labels)), frameon=False, fontsize=9, title_fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    return fig
