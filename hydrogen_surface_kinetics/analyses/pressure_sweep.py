"""Pressure sweep: the rate against the probe-gas pressure on a fine log grid, at fixed
temperature (documentation: Stage 11.2b). The pressure counterpart of the temperature sweep.

At each temperature (default: pressure_temperatures) the probe-gas partial pressure runs over
a log grid and every impurity keeps its ratio to it (impurity_ratio), so the gas has the same
purity at every pressure and only the total pressure changes.

Reports the rate against p (log-log); the local slope d ln TOF / d ln p (for the reference
case the order in the probe gas; for impurity cases every partial pressure moves together,
so it is the order in the total pressure at fixed purity); the poisoning factor against p,
and the half-rate pressure where it crosses 0.5.
Outputs: pressure_sweep.csv, half_rate_pressure.csv, pressure_sweep.png,
pressure_sweep_coverages_<T>K.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import probe_label
from . import common
from .composition import coverage_table, plot_coverage_bars


def sweep_pressures(ctx):
    """Pa: the section's pressures, or a log grid with per_decade points spanning the
    pressure series (plus reference_pressure); without a pressure series, two decades either
    side of reference_pressure."""
    cfg, s = ctx.cfg, ctx.cfg.pressure_sweep
    if s['pressures'] is not None:
        return sorted(float(p) for p in s['pressures'])
    span = set(cfg.pressures) | {cfg.reference_pressure}
    lo, hi = min(span), max(span)
    if lo == hi:
        lo, hi = lo / 100.0, hi * 100.0
    n = int(round(s['per_decade'] * np.log10(hi / lo))) + 1
    return list(np.logspace(np.log10(lo), np.log10(hi), max(n, 2)))


def sweep_temperatures(ctx):
    s = ctx.cfg.pressure_sweep
    return [float(t) for t in (ctx.cfg.pressure_temperatures if s['temperatures'] is None else s['temperatures'])]


def sweep_conditions(ctx):
    ratio = common.setting(ctx, 'pressure_sweep', 'ratio', ctx.cfg.impurity_ratio)
    cases = common.sweep_cases(ctx, 'pressure_sweep')
    rows = []
    for T in sweep_temperatures(ctx):
        for p in sweep_pressures(ctx):
            pt = {'T': T, 'P': p, 'series': 'pressure'}
            rows += [common.fixed_probe_row(ctx, f'Psweep r={C.ratio_tag(ratio)}', case, pt, impurity_ratio=ratio)
                     for case in cases]
    return pd.DataFrame(rows)


def local_slope(df):
    """d ln TOF / d ln p by central differences on each (case, T) curve of converged points."""
    df['slope_dlnTOF_dlnp'] = np.nan
    for _, part in df[df['converged'].astype(bool)].groupby(['case', 'T_K']):
        part = part.sort_values('P_ref_bar')
        if len(part) >= 3:
            slope = np.gradient(np.log(part['tof_fwd'].to_numpy()), np.log(part['P_ref_bar'].to_numpy()))
            df.loc[part.index, 'slope_dlnTOF_dlnp'] = slope
    return df


def half_rate_pressure(part):
    """First pressure (going up) at which the poisoning factor crosses 0.5, by log-log
    interpolation, with the direction of the crossing."""
    part = part[part['converged'].astype(bool)].sort_values('P_ref_bar')
    if part.empty:
        return pd.Series({'half_rate_P_bar': np.nan, 'direction': '', 'note': 'no converged points'})
    x, y = np.log10(part['P_ref_bar'].to_numpy()), np.log10(part['poisoning_factor'].to_numpy())
    above = y >= np.log10(0.5)
    if above.all():
        return pd.Series({'half_rate_P_bar': np.nan, 'direction': '', 'note': 'never below half in the range'})
    if not above.any():
        return pd.Series({'half_rate_P_bar': np.nan, 'direction': '', 'note': 'below half over the whole range'})
    j = int(np.nonzero(above[1:] != above[:-1])[0][0]) + 1
    frac = (np.log10(0.5) - y[j - 1]) / (y[j] - y[j - 1])
    return pd.Series({'half_rate_P_bar': 10 ** (x[j - 1] + frac * (x[j] - x[j - 1])),
                      'direction': 'recovers above half' if above[j] else 'falls below half', 'note': ''})


def run(ctx):
    cfg = ctx.cfg
    ctx.log('\n== pressure_sweep: pressure sweep ==')
    results = ctx.run_conditions(sweep_conditions(ctx), prefix='')
    df = common.results_frame(ctx, results)
    df = common.add_poisoning(df, cfg.reference_label, ctx.model.probe_sites, by=('T_K', 'P_ref_bar'))
    df = local_slope(df.sort_values(['case', 'T_K', 'P_ref_bar']).reset_index(drop=True))
    ctx.save_table(df, 'pressure_sweep.csv')
    half = (df[df['case'] != cfg.reference_label].groupby(['T_K', 'case'], sort=False)
            .apply(half_rate_pressure, include_groups=False).reset_index())
    ctx.save_table(half, 'half_rate_pressure.csv')
    ctx.log('  probe-gas pressure at which the rate crosses half of the reference:\n'
            + half.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    bad = common.problems(df)
    if len(bad):
        ctx.log('  needs attention:\n' + bad[['T_K', 'P_ref_bar', 'case', 'converged', 'starts_agree']]
                .to_string(index=False))
    P.save_parts(plot_pressure_sweep(ctx, df, half), ctx.path('pressure_sweep.png'), cfg)
    cov = coverage_table(ctx, results).droplevel(['basis', 'series'])
    cases = [c for c in common.sweep_cases(ctx, 'pressure_sweep') if c in set(cov.index.get_level_values('case'))]
    for T in sweep_temperatures(ctx):
        part = cov[cov.index.get_level_values('T_K') == T].droplevel('T_K')
        if part.empty:
            continue
        figs = plot_coverage_bars(ctx, part, cases, [p / 1e5 for p in sweep_pressures(ctx)],
                                  f'p({probe_label(ctx.model)}) (bar)',
                                  f'What covers the surface against pressure, {T:.0f} K', log_x=True)
        P.save_parts(figs, ctx.path(f'pressure_sweep_coverages_{T:.0f}K.png'), cfg)
    ctx.outputs['pressure_sweep'] = {'results': results, 'sweep': df, 'half_rate': half, 'coverages': cov}
    return ctx.outputs['pressure_sweep']


def plot_pressure_sweep(ctx, sweep, half=None):
    """One row per temperature: (a) forward TOF against p (log-log), (b) local slope,
    (c) poisoning factor with the half-rate pressures. Returns the figures (at most three
    rows each)."""
    return [_pressure_sweep_page(ctx, sweep, half, temps)
            for temps in P.pages(sorted(sweep['T_K'].unique()), P.rows_per_page(3))]


def _pressure_sweep_page(ctx, sweep, half, temps):
    cfg = ctx.cfg
    probe = probe_label(ctx.model)
    cases = common.sweep_cases(ctx, 'pressure_sweep')
    colors = P.colors_for(cases, first=cfg.reference_label)
    fig, axes = plt.subplots(len(temps), 3, figsize=(15, 4.4 * len(temps)), squeeze=False)
    for row, T in zip(axes, temps):
        for case in cases:
            d = sweep[(sweep['case'] == case) & (sweep['T_K'] == T)].sort_values('P_ref_bar')
            if not d['converged'].astype(bool).any():
                continue
            two = (d['starts_agree'] == False).to_numpy()              # noqa: E712
            style = dict(color=colors[case], ms=4.5, lw=1.6, label=case)
            row[0].plot(d['P_ref_bar'], d['tof_fwd'], '-o', **style)
            row[0].plot(d['P_ref_bar'].to_numpy()[two], d['tof_fwd'].to_numpy()[two], 'o', ms=9, mfc='none',
                        mec=P.INK, mew=1.2)
            row[1].plot(d['P_ref_bar'], d['slope_dlnTOF_dlnp'], '-o', **style)
            if case != cfg.reference_label:
                row[2].plot(d['P_ref_bar'], d['poisoning_factor'], '-o', **style)
                if half is not None:
                    h = half[(half['T_K'] == T) & (half['case'] == case) & half['half_rate_P_bar'].notna()]
                    if len(h):
                        row[2].plot(h['half_rate_P_bar'].iloc[0], 0.5, **P.HALF_MARK, mec=colors[case])
        for ax in row:
            ax.set_xscale('log')
            ax.set_xlabel(f'p({probe}) (bar)')
        row[0].set_yscale('log')
        P.style_axis(row[0], f'{T:.0f} K: forward {probe} dissociation rate', 'TOF (1/s)')
        row[1].axhline(0, color=P.MUTED, lw=0.8)
        P.style_axis(row[1], f'{T:.0f} K: local slope', 'd ln TOF / d ln p')
        row[2].axhline(0.5, color=P.MUTED, lw=1, ls='--')
        row[2].set_yscale('log')
        P.style_axis(row[2], f'{T:.0f} K: poisoning factor', f'TOF / TOF({cfg.reference_label})')
    handles, labels = axes[0][0].get_legend_handles_labels()
    if half is not None and half['half_rate_P_bar'].notna().any():
        handles.append(P.half_mark_handle('half-rate pressure'))
        labels.append('half-rate pressure')
    handles.append(P.two_states_handle())
    labels.append('two steady states')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(labels), frameon=False,
               fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94 if len(temps) == 1 else 0.97))
    return fig
