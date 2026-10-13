"""Poisoning curves: rate vs impurity level (documentation: Stage 11a).

The probe-gas pressure stays fixed and each impurity is added at a ratio p(imp)/p(probe).
Outputs: poisoning_curves.csv, poisoning_tolerance.csv, poisoning_curves.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import plotting as P
from ..labels import probe_label
from . import common


def curve_conditions(ctx):
    s, ref = ctx.cfg.curves, ctx.cfg.reference_label
    rows = []
    for T in common.temperatures(ctx, 'curves'):
        rows.append(common.fixed_probe_row(ctx, 'curve r=0', ref, T, impurity_ratio=0.0, ratio=0.0))
        for imp in common.impurity_labels(ctx, 'curves'):
            for r in s['ratios']:
                rows.append(common.fixed_probe_row(ctx, f'curve r={r:.3g}', imp, T, impurity_ratio=r, ratio=r))
    return pd.DataFrame(rows)


def half_poisoning_ratio(part):
    """Impurity ratio at which the forward TOF halves (log-log interpolation)."""
    part = part[part['converged'].astype(bool)].sort_values('ratio')
    x, y = np.log10(part['ratio'].to_numpy()), np.log10(part['poisoning_factor'].to_numpy())
    below = np.nonzero(y <= np.log10(0.5))[0]
    if len(below) == 0:
        return pd.Series({'ratio_at_half_rate': np.nan, 'note': 'never halved in the range tested'})
    j = below[0]
    if j == 0:
        return pd.Series({'ratio_at_half_rate': 10 ** x[0], 'note': 'already halved at the lowest level: at most this'})
    frac = (np.log10(0.5) - y[j - 1]) / (y[j] - y[j - 1])
    return pd.Series({'ratio_at_half_rate': 10 ** (x[j - 1] + frac * (x[j] - x[j - 1])), 'note': ''})


def run(ctx):
    cfg = ctx.cfg
    ctx.log('\n== curves: poisoning curves ==')
    conds = curve_conditions(ctx)
    results = ctx.run_conditions(conds, prefix='')
    df = common.results_frame(ctx, results, extra_cols=('ratio',))
    df = common.add_poisoning(df, cfg.reference_label, ctx.model.probe_sites, by=('T_K',))
    df = df.sort_values(['T_K', 'case', 'ratio']).reset_index(drop=True)
    ctx.save_table(df, 'poisoning_curves.csv')
    tol = (df[df['case'] != cfg.reference_label].groupby(['T_K', 'case'], sort=False)
           .apply(half_poisoning_ratio, include_groups=False).reset_index())
    tol['ppm_of_probe'] = tol['ratio_at_half_rate'] * 1e6
    ctx.save_table(tol, 'poisoning_tolerance.csv')
    ctx.log('  impurity level that halves the forward rate:\n'
            + tol.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    bad = common.problems(df)
    if len(bad):
        ctx.log('  needs attention:\n' + bad[['T_K', 'case', 'ratio', 'converged', 'starts_agree']].to_string(index=False))
    fig = plot_poisoning_curves(ctx, df, tol)
    P.save(fig, ctx.path('poisoning_curves.png'), cfg)
    plt.close(fig)
    ctx.outputs['curves'] = {'results': results, 'curves': df, 'tolerance': tol}
    return ctx.outputs['curves']


def plot_poisoning_curves(ctx, curves, tolerance=None):
    cfg = ctx.cfg
    probe = probe_label(ctx.model)
    imps = common.impurity_labels(ctx, 'curves')
    colors = P.colors_for(imps)
    temps = sorted(curves['T_K'].unique())
    fig, axes = plt.subplots(len(temps), 2, figsize=(11, 3.9 * len(temps)), squeeze=False)
    for row, T in zip(axes, temps):
        part = curves[(curves['T_K'] == T) & (curves['case'] != cfg.reference_label)]
        for imp in imps:
            d = part[(part['case'] == imp) & part['converged'].astype(bool)].sort_values('ratio')
            if d.empty:
                continue
            two = (d['starts_agree'] == False).to_numpy()            # noqa: E712
            for ax, col in zip(row, ['poisoning_factor', 'theta_impurity']):
                y = P.log_floor(d[col]) if col == 'theta_impurity' else d[col].to_numpy()
                ax.plot(d['ratio'], y, '-o', color=colors[imp], ms=4.5, lw=1.6, label=imp)
                ax.plot(d['ratio'].to_numpy()[two], y[two], 'o', ms=9, mfc='none', mec=P.INK, mew=1.2)
            if tolerance is not None:                     # interpolated half-rate level (not the bounds)
                t = tolerance[(tolerance['T_K'] == T) & (tolerance['case'] == imp) & (tolerance['note'].fillna('') == '')]
                if len(t) and np.isfinite(t['ratio_at_half_rate'].iloc[0]):
                    row[0].plot(t['ratio_at_half_rate'].iloc[0], 0.5, **P.HALF_MARK, mec=colors[imp])
        row[0].axhline(0.5, color=P.MUTED, lw=1, ls='--')
        row[0].text(min(cfg.curves['ratios']), 0.45, 'rate halved', color=P.MUTED, fontsize=8, ha='left', va='top')
        P.style_axis(row[0], f'{T:.0f} K: forward {probe} dissociation rate', f'poisoning factor  TOF / TOF({cfg.reference_label})')
        P.style_axis(row[1], f'{T:.0f} K: sites taken by the impurity', f'fraction of sites (not empty, {probe} or its adsorbates)')
        for ax in row:
            ax.set_xscale('log')
            ax.set_yscale('log')
            ax.set_xlabel(f'impurity level  p(impurity) / p({probe})')
    handles, labels = axes[0][0].get_legend_handles_labels()
    if tolerance is not None and (tolerance['note'].fillna('') == '').any():
        handles.append(P.half_mark_handle('half-rate level'))
        labels.append('half-rate level')
    handles.append(P.two_states_handle())
    labels.append('two steady states (starts disagree)')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=len(labels), frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return fig
