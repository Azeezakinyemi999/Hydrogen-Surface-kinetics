"""Poisoning over time (documentation: Stage 11.6).

At every point of the temperature and pressure series: start from the steady state in the
probe gas alone (a working catalyst), switch the gas to probe + impurity at t = 0, and
follow the forward rate as the impurity builds up.
Outputs: poisoning_over_time.csv, poisoning_over_time_trajectories.csv, and figures
poisoning_over_time_<T>K_<P>bar.png (one per point; or poisoning_over_time_<impurity>_vs_T/_vs_P.png,
one per impurity, when a series has more points than impurities).
"""
import re

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
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


def reference_start(ctx, pt):
    """Coverages of the probe-only steady state at the point: from the main run if it has the
    point, otherwise solved here."""
    cfg = ctx.cfg
    summary = ctx.main()['summary']
    hit = summary[(summary['basis'] == 'fixed_probe') & (summary['case'] == cfg.reference_label)
                  & (summary['T_K'] == pt['T']) & np.isclose(summary['P_ref_bar'].astype(float), pt['P'] / 1e5)]
    if len(hit):
        res = ctx.main()['results'][hit['key'].iloc[0]]
    else:
        conds = pd.DataFrame([common.fixed_probe_row(ctx, 'fixed_probe', cfg.reference_label, pt)])
        res = next(iter(ctx.run_conditions(conds, prefix='').values()))
    return S.coverages_to_dict(ctx.model, res['theta_a']) if common.converged(res) else None


def run(ctx):
    model, cfg, s = ctx.model, ctx.cfg, ctx.cfg.transients
    ctx.log('\n== transients: poisoning over time ==')
    tt = times(ctx)
    imps = common.impurity_labels(ctx, 'transients')
    jobs = []
    for pt in common.points(ctx, 'transients'):
        start = reference_start(ctx, pt)
        if start is None:
            ctx.log(f'  skipping {common.point_name(pt)}: no converged {cfg.reference_label} steady state')
            continue
        for imp in imps:
            for ratio in s['ratios']:
                row = common.fixed_probe_row(ctx, 'transient', imp, pt, impurity_ratio=ratio, ratio=ratio)
                key = f"transient r={ratio:.3g} | {imp} | {pt['T']:.0f} K | {pt['P'] / 1e5:g} bar"
                jobs.append(common.job(key, row, theta0=start, times=tt))
    ctx.log(f'  {len(jobs)} transient runs')
    res_all = ctx.run_jobs(lambda job: run_transient(model, job), jobs, 'cache_transients',
                           describe=lambda res: res['error'] or '')
    rows, trajs = [], []
    for job in jobs:
        res = res_all[job['key']]
        r = {**job['row'], **job['_labels']}
        where = dict(series=r['series'], T_K=r['T_K'], P_ref_bar=r['P_ref_bar'], case=r['case'], ratio=r['ratio'])
        out = {**where, 'error': res.get('error') if res['status'] == 'done' else res['status']}
        traj = res.get('trajectory')
        if traj is not None and len(traj):
            tof0 = res['tof0']
            traj = traj.assign(**where)
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
        drawn = set()
        for series in C.SERIES:
            for stem, figs in plot_poisoning_over_time(ctx, jobs, res_all, series, skip=drawn):
                P.save_parts(figs, ctx.path(f'{stem}.png'), cfg)
                drawn.add(stem)
    ctx.outputs['transients'] = {'results': res_all, 'summary': summary}
    return ctx.outputs['transients']


def plot_poisoning_over_time(ctx, jobs, results, series='temperature', skip=()):
    """The figures of one series as [(file stem, figures)]: one per point with the impurities
    across, or one per impurity with the points across when the series has more points than
    impurities. Each column (at most P.MAX_COLS per figure) is one (point, impurity): the rate
    relative to t = 0 (top) and the sites the impurity takes (bottom). Stems in `skip` (points
    of both series, already drawn) are left out."""
    mine = [j for j in jobs if j['_labels']['series'] in (series, 'both')]
    pts = sorted({(j['row']['T_K'], j['_labels']['P_ref_bar']) for j in mine})
    imps = [i for i in common.impurity_labels(ctx, 'transients') if any(j['row']['case'] == i for j in mine)]
    if not pts or not imps:
        return []
    if len(imps) >= len(pts):
        groups = [(f'poisoning_over_time_{T:.0f}K_{Pb:g}bar', [(T, Pb, imp) for imp in imps]) for T, Pb in pts]
    else:
        slug = lambda text: re.sub(r'[^A-Za-z0-9.+-]+', '_', text).strip('_')     # noqa: E731
        groups = [(f'poisoning_over_time_{slug(imp)}_{common.SERIES_SUFFIX[series]}', [(T, Pb, imp) for T, Pb in pts])
                  for imp in imps]
    return [(stem, [_poisoning_over_time_page(ctx, mine, results, page) for page in P.pages(panels, P.MAX_COLS)])
            for stem, panels in groups if stem not in skip]


def _poisoning_over_time_page(ctx, jobs, results, panels):
    s = ctx.cfg.transients
    probe = probe_label(ctx.model)
    shades = ['#9ec5f4', '#5598e7', '#256abf', '#0d366b', '#061b36']        # light -> dark = more impurity
    ratio_colors = dict(zip(sorted(s['ratios']), shades if len(s['ratios']) <= len(shades) else P.SERIES_COLORS))
    fig, axes = plt.subplots(2, len(panels), figsize=(3.4 * len(panels) + 1.0, 5.8), squeeze=False, sharex=True)
    for c_, (T, Pb, imp) in enumerate(panels):
        ax_r, ax_c = axes[0][c_], axes[1][c_]
        for job in jobs:
            r = job['row']
            if r['T_K'] != T or job['_labels']['P_ref_bar'] != Pb or r['case'] != imp:
                continue
            res = results[job['key']]
            traj = res.get('trajectory')
            if traj is None or not len(traj):
                continue
            c = ratio_colors.get(r['ratio'], P.INK)
            ax_r.plot(traj['t_s'], traj['tof_fwd'] / res['tof0'], color=c, lw=1.8)
            ax_c.plot(traj['t_s'], P.log_floor(traj['theta_impurity']), color=c, lw=1.8)
        ax_r.axhline(0.5, color=P.MUTED, lw=0.8, ls='--')
        for ax in (ax_r, ax_c):
            ax.set_xscale('log')
            ax.set_yscale('log')
        lo, hi = ax_r.get_ylim()
        ax_r.set_ylim(min(lo, 0.1), max(hi, 2.0))      # at least a decade: whole-decade ticks only
        P.style_axis(ax_r, f'{imp}, {T:.0f} K, {Pb:g} bar', 'TOF / TOF(t=0)' if c_ == 0 else '')
        P.style_axis(ax_c, '', 'sites taken by impurity' if c_ == 0 else '')
        ax_c.set_xlabel('time after the impurity arrives (s)')
    handles = [plt.Line2D([], [], color=c, lw=1.8) for c in ratio_colors.values()]
    fig.legend(handles, [f'{r:g}' for r in ratio_colors], title=f'p(impurity) / p({probe})', loc='upper center',
               bbox_to_anchor=(0.5, 1.03), ncol=len(handles), frameon=False, fontsize=9, title_fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig
