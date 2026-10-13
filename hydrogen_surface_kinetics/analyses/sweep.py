"""Temperature sweep: Arrhenius plot, apparent activation energy, release temperature
(documentation: Stage 11.2). Always at reference_pressure: it is the fine-grained
temperature series (pressure_sweep is its pressure counterpart). Outputs: temperature_sweep.csv, release_temperature.csv,
temperature_sweep.png, and sweep_coverages.png (what covers the surface against T).
"""
import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import probe_label, species_label
from . import common
from .composition import coverage_table, plot_coverage_bars


def sweep_temperatures(ctx):
    s = ctx.cfg.sweep
    if s['temperatures'] is not None:
        return [float(t) for t in s['temperatures']]
    lo, hi = min(ctx.cfg.temperatures), max(ctx.cfg.temperatures)
    return list(np.arange(lo, hi + 0.5 * s['step'], s['step']))


def sweep_conditions(ctx):
    ratio = common.setting(ctx, 'sweep', 'ratio', ctx.cfg.impurity_ratio)
    cases = common.sweep_cases(ctx, 'sweep')
    return pd.DataFrame([common.fixed_probe_row(ctx, f'Tsweep r={C.ratio_tag(ratio)}', case, common.reference_point(ctx, T),
                                                impurity_ratio=ratio)
                         for T in sweep_temperatures(ctx) for case in cases])


def apparent_activation_energy(df):
    """E_app = -R d(ln TOF)/d(1/T), central differences, per case (kJ/mol)."""
    df['E_app_kJ_per_mol'] = np.nan
    for case, part in df[df['converged'].astype(bool)].groupby('case'):
        if len(part) >= 3:
            dlnr = np.gradient(np.log(part['tof_fwd'].to_numpy()), 1.0 / part['T_K'].to_numpy())
            df.loc[part.index, 'E_app_kJ_per_mol'] = -ct.gas_constant * dlnr / 1e6
    return df


def release_temperature(part):
    """Lowest T at which a poisoned rate recovers to half of the reference rate."""
    part = part[part['converged'].astype(bool)].sort_values('T_K')
    T, y = part['T_K'].to_numpy(), np.log10(part['poisoning_factor'].to_numpy())
    half = np.log10(0.5)
    if len(T) == 0:
        return pd.Series({'release_T_K': np.nan, 'note': 'no converged points'})
    if y[0] >= half:
        return pd.Series({'release_T_K': np.nan, 'note': 'rate never below half (not poisoning)'})
    above = np.nonzero(y >= half)[0]
    if len(above) == 0:
        return pd.Series({'release_T_K': np.nan, 'note': f'still poisoned at {T[-1]:.0f} K'})
    j = above[0]
    frac = (half - y[j - 1]) / (y[j] - y[j - 1])
    return pd.Series({'release_T_K': T[j - 1] + frac * (T[j] - T[j - 1]), 'note': ''})


def run(ctx):
    cfg = ctx.cfg
    ctx.log('\n== sweep: temperature sweep ==')
    results = ctx.run_conditions(sweep_conditions(ctx), prefix='')
    df = common.results_frame(ctx, results)
    df = common.add_poisoning(df, cfg.reference_label, ctx.model.probe_sites, by=('T_K', 'P_ref_bar'))
    df = apparent_activation_energy(df.sort_values(['case', 'T_K']).reset_index(drop=True))
    ctx.save_table(df, 'temperature_sweep.csv')
    rel = (df[df['case'] != cfg.reference_label].groupby('case', sort=False)
           .apply(release_temperature, include_groups=False).reset_index())
    ctx.save_table(rel, 'release_temperature.csv')
    ctx.log('  temperature at which the rate recovers to half of the reference:\n'
            + rel.to_string(index=False, float_format=lambda v: f'{v:.0f}'))
    bad = common.problems(df)
    if len(bad):
        ctx.log('  needs attention:\n' + bad[['T_K', 'case', 'converged', 'starts_agree']].to_string(index=False))
    fig = plot_temperature_sweep(ctx, df, rel)
    P.save(fig, ctx.path('temperature_sweep.png'), cfg)
    plt.close(fig)
    cov = coverage_table(ctx, results).droplevel(['basis', 'series', 'P_ref_bar'])
    fig = plot_sweep_coverages(ctx, cov)
    P.save(fig, ctx.path('sweep_coverages.png'), cfg)
    plt.close(fig)
    ctx.outputs['sweep'] = {'results': results, 'sweep': df, 'release': rel, 'coverages': cov}
    return ctx.outputs['sweep']


def plot_temperature_sweep(ctx, sweep, release=None):
    cfg = ctx.cfg
    probe = probe_label(ctx.model)
    cases = common.sweep_cases(ctx, 'sweep')
    colors = P.colors_for(cases, first=cfg.reference_label)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.4))
    for case in cases:
        d = sweep[sweep['case'] == case].sort_values('T_K')        # unconverged rows are NaN: gaps
        if not d['converged'].astype(bool).any():
            continue
        two = (d['starts_agree'] == False).to_numpy()              # noqa: E712
        style = dict(color=colors[case], ms=4.5, lw=1.6, label=case)
        axes[0].plot(1000.0 / d['T_K'], d['tof_fwd'], '-o', **style)
        axes[0].plot(1000.0 / d['T_K'].to_numpy()[two], d['tof_fwd'].to_numpy()[two], 'o', ms=9, mfc='none',
                     mec=P.INK, mew=1.2)
        axes[1].plot(d['T_K'], d['E_app_kJ_per_mol'], '-o', **style)
        if case != cfg.reference_label:
            axes[2].plot(d['T_K'], d['poisoning_factor'], '-o', **style)
            if release is not None:                       # where the rate recovers to half
                r = release[(release['case'] == case) & release['release_T_K'].notna()]
                if len(r):
                    axes[2].plot(r['release_T_K'].iloc[0], 0.5, **P.HALF_MARK, mec=colors[case])
    axes[0].set_yscale('log')
    axes[0].set_xlabel('1000 / T  (1/K)')
    P.style_axis(axes[0], 'a) Arrhenius plot', f'forward {probe} dissociation TOF (1/s)')
    axes[1].axhline(0, color=P.MUTED, lw=0.8)
    axes[1].set_xlabel('T (K)')
    P.style_axis(axes[1], 'b) Local apparent activation energy', '$E_{app}$ (kJ/mol)')
    axes[2].axhline(0.5, color=P.MUTED, lw=1, ls='--')
    axes[2].text(sweep['T_K'].max(), 0.45, 'rate halved', color=P.MUTED, fontsize=8, va='top', ha='right')
    axes[2].set_yscale('log')
    axes[2].set_xlabel('T (K)')
    P.style_axis(axes[2], 'c) Poisoning factor vs temperature', f'TOF / TOF({cfg.reference_label})')
    handles, labels = axes[0].get_legend_handles_labels()
    if release is not None and release['release_T_K'].notna().any():
        handles.append(P.half_mark_handle('release temperature'))
        labels.append('release temperature')
    handles.append(P.two_states_handle())
    labels.append('two steady states')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.04), ncol=len(labels), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig


def plot_sweep_coverages(ctx, cov):
    """Stacked coverages against temperature, one bar per solved temperature and one panel per
    case (temperatures without a converged steady state are gaps)."""
    cases = [c for c in common.sweep_cases(ctx, 'sweep') if c in set(cov.index.get_level_values('case'))]
    return plot_coverage_bars(ctx, cov, cases, sweep_temperatures(ctx), 'T (K)',
                              'What covers the surface against temperature')
