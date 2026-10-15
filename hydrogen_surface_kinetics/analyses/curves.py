"""Poisoning curves: rate vs impurity level (documentation: Stage 11.1).

At every point of the temperature and pressure series, the probe-gas pressure stays at the
point's value and each impurity is added at a ratio p(imp)/p(probe).
Where the two starting surfaces reach different steady states, the figures draw the start-B
branch too (dashed) and shade the window: the two-steady-state picture of Stage 5.
site_blocking.png splits the poisoning factor into its coverage and site-count parts
(Stage 7) at one point.
Outputs: poisoning_curves.csv, poisoning_tolerance.csv, poisoning_curves_vs_T.png and,
with a pressure series, poisoning_curves_vs_P.png; site_blocking.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import probe_label
from . import common


def curve_conditions(ctx):
    s, ref = ctx.cfg.curves, ctx.cfg.reference_label
    rows = []
    for pt in common.points(ctx, 'curves'):
        rows.append(common.fixed_probe_row(ctx, 'curve r=0', ref, pt, impurity_ratio=0.0, ratio=0.0))
        for imp in common.impurity_labels(ctx, 'curves'):
            for r in s['ratios']:
                rows.append(common.fixed_probe_row(ctx, f'curve r={r:.3g}', imp, pt, impurity_ratio=r, ratio=r))
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
    df = common.add_poisoning(df, cfg.reference_label, ctx.model.probe_sites, by=('T_K', 'P_ref_bar'))
    df = df.sort_values(['T_K', 'P_ref_bar', 'case', 'ratio']).reset_index(drop=True)
    ctx.save_table(df, 'poisoning_curves.csv')
    tol = (df[df['case'] != cfg.reference_label].groupby(['series', 'T_K', 'P_ref_bar', 'case'], sort=False)
           .apply(half_poisoning_ratio, include_groups=False).reset_index())
    tol['ppm_of_probe'] = tol['ratio_at_half_rate'] * 1e6
    ctx.save_table(tol, 'poisoning_tolerance.csv')
    ctx.log('  impurity level that halves the forward rate:\n'
            + tol.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    bad = common.problems(df)
    if len(bad):
        ctx.log('  needs attention:\n' + bad[['T_K', 'P_ref_bar', 'case', 'ratio', 'converged', 'starts_agree']]
                .to_string(index=False))
    for series in C.SERIES:
        P.save_parts(plot_poisoning_curves(ctx, df, tol, series),
                     ctx.path(common.series_file('poisoning_curves', series)), cfg)
    P.save_parts(plot_site_blocking(ctx, df), ctx.path('site_blocking.png'), cfg)
    ctx.outputs['curves'] = {'results': results, 'curves': df, 'tolerance': tol}
    return ctx.outputs['curves']


def plot_poisoning_curves(ctx, curves, tolerance=None, series='temperature'):
    """One column per point of the series, at most P.MAX_COLS per figure: the rate relative
    to the reference (top) and the sites the impurity takes (bottom), against the impurity
    level. Returns the figures."""
    sub = C.in_series(curves, series)
    pts = list(sub[['T_K', 'P_ref_bar']].drop_duplicates().sort_values(['T_K', 'P_ref_bar']).itertuples(index=False))
    return [_poisoning_curves_page(ctx, sub, tolerance, page) for page in P.pages(pts, P.MAX_COLS)]


def _poisoning_curves_page(ctx, sub, tolerance, pts):
    cfg = ctx.cfg
    probe = probe_label(ctx.model)
    imps = common.impurity_labels(ctx, 'curves')
    colors = P.colors_for(imps)
    fig, axes = plt.subplots(2, len(pts), figsize=(3.6 * len(pts) + 1.2, 7.4), squeeze=False, sharex=True)
    for j, (T, Pb) in enumerate(pts):
        row = (axes[0][j], axes[1][j])                    # this point's column: rate, then sites
        part = sub[(sub['T_K'] == T) & (sub['P_ref_bar'] == Pb) & (sub['case'] != cfg.reference_label)]
        for imp in imps:
            d = part[(part['case'] == imp) & part['converged'].astype(bool)].sort_values('ratio')
            if d.empty:
                continue
            two = (d['starts_agree'] == False).to_numpy()            # noqa: E712
            for ax, col in zip(row, ['poisoning_factor', 'theta_impurity']):
                y = P.log_floor(d[col]) if col == 'theta_impurity' else d[col].to_numpy()
                ax.plot(d['ratio'], y, '-o', color=colors[imp], ms=4.5, lw=1.6, label=imp)
                ax.plot(d['ratio'].to_numpy()[two], y[two], 'o', ms=9, mfc='none', mec=P.INK, mew=1.2)
                col_b = col + '_B' if col == 'theta_impurity' else 'poisoning_factor_B'
                if two.any() and col_b in d:                  # the covered start's branch, and the window
                    yb = d[col_b].to_numpy(dtype=float)
                    yb = P.log_floor(yb) if col == 'theta_impurity' else yb
                    ax.plot(d['ratio'], yb, '--', color=colors[imp], lw=1.4)
                    ax.plot(d['ratio'].to_numpy()[two], yb[two], 's', color=colors[imp], ms=4.5, mfc='white')
                    r2 = d['ratio'].to_numpy()[two]
                    ax.axvspan(r2.min() / 1.3, r2.max() * 1.3, color=colors[imp], alpha=0.10, lw=0, zorder=0)
            if tolerance is not None:                     # interpolated half-rate level (not the bounds)
                t = tolerance[(tolerance['T_K'] == T) & (tolerance['P_ref_bar'] == Pb) & (tolerance['case'] == imp)
                              & (tolerance['note'].fillna('') == '')]
                if len(t) and np.isfinite(t['ratio_at_half_rate'].iloc[0]):
                    row[0].plot(t['ratio_at_half_rate'].iloc[0], 0.5, **P.HALF_MARK, mec=colors[imp])
        row[0].axhline(0.5, color=P.MUTED, lw=1, ls='--')
        P.style_axis(row[0], f'{T:.0f} K, {Pb:g} bar',
                     f'poisoning factor  TOF / TOF({cfg.reference_label})' if j == 0 else '')
        P.style_axis(row[1], '', 'impurity-held sites' if j == 0 else '')
        for ax in row:
            ax.set_xscale('log')
            ax.set_yscale('log')
        row[1].set_xlabel(f'p(impurity) / p({probe})')
    shown = [imp for imp in imps if (sub['case'] == imp).any()]
    handles = [plt.Line2D([], [], color=colors[imp], marker='o', ms=4.5, lw=1.6) for imp in shown]
    handles.append(plt.Line2D([], [], color=P.MUTED, lw=1, ls='--'))
    labels = list(shown) + ['rate halved']
    if tolerance is not None and (tolerance['note'].fillna('') == '').any():
        handles.append(P.half_mark_handle('half-rate level'))
        labels.append('half-rate level')
    handles.append(P.two_states_handle())
    labels.append('two steady states (starts disagree)')
    disagree = sub['starts_agree'] == False                         # noqa: E712
    if 'poisoning_factor_B' in sub and (disagree & sub['poisoning_factor_B'].notna()).any():
        handles.append(plt.Line2D([], [], ls='--', color=P.INK, marker='s', ms=4.5, mfc='white'))
        labels.append('covered start (start B)')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=min(len(labels), 5),
               frameon=False, fontsize=9)
    fig.suptitle(f'Poisoning curves: forward {probe} dissociation rate (top) and impurity-held sites (bottom)',
                 x=0.01, y=0.9, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.86))
    return fig


def blocking_rows(ctx, curves):
    """Rows of the site-blocking figure: the curves at curves['blocking_point'] (default: the
    first series point), for the impurities that change the rate or take sites there."""
    cfg = ctx.cfg
    bp = cfg.curves['blocking_point']
    pt = common.points(ctx, 'curves')[0] if bp is None else {'T': float(bp[0]), 'P': float(bp[1])}
    at = curves[(curves['T_K'] == pt['T']) & np.isclose(curves['P_ref_bar'].astype(float), pt['P'] / 1e5)
                & curves['converged'].astype(bool) & (curves['case'] != cfg.reference_label)]
    shown = [imp for imp in common.impurity_labels(ctx, 'curves')
             if len(g := at[at['case'] == imp]) and (g['poisoning_factor'].min() < 0.99 or g['theta_impurity'].max() > 0.01)]
    return pt, at, shown or [imp for imp in common.impurity_labels(ctx, 'curves') if (at['case'] == imp).any()]


def plot_site_blocking(ctx, curves):
    """Stage 7 on real data, one row per impurity at one point: (a) empty sites and sites the
    impurity takes, (b) the site-count factor f, (c) the poisoning factor from rates (markers)
    against the full site-blocking prediction (line) and the coverage-only one (dashed).
    Returns the figures (at most three rows each)."""
    pt, at, shown = blocking_rows(ctx, curves)
    return [_site_blocking_page(ctx, pt, at, page) for page in P.pages(shown, P.rows_per_page(3))]


def _site_blocking_page(ctx, pt, at, shown):
    probe = probe_label(ctx.model)
    colors = P.colors_for(common.impurity_labels(ctx, 'curves'))
    fig, axes = plt.subplots(len(shown), 3, figsize=(13, 3.2 * len(shown)), squeeze=False)
    for row, imp in zip(axes, shown):
        d = at[at['case'] == imp].sort_values('ratio')
        c = colors[imp]
        row[0].plot(d['ratio'], P.log_floor(d['theta_free']), '-o', color=P.MUTED, ms=4, lw=1.4, label='empty sites')
        row[0].plot(d['ratio'], P.log_floor(d['theta_impurity']), '-o', color=c, ms=4, lw=1.6, label='impurity-held')
        row[0].set_yscale('log')
        row[0].legend(loc='best', frameon=False, fontsize=8)
        row[1].plot(d['ratio'], d['site_count_factor'], '-o', color=c, ms=4, lw=1.6)
        row[1].axhline(1.0, color=P.MUTED, lw=0.8, ls='--')
        row[2].plot(d['ratio'], d['blocking_full'], '-', color=c, lw=1.6, label='site blocking, with f')
        row[2].plot(d['ratio'], d['blocking_from_theta'], '--', color=c, lw=1.4, label='coverages only')
        row[2].plot(d['ratio'], d['poisoning_factor'], 'o', color=c, ms=5, mec='white', mew=0.6, label='from rates')
        row[2].set_yscale('log')
        row[2].legend(loc='lower left', frameon=False, fontsize=8)
        for ax, desc, ylabel in zip(row, ('sites', 'site-count factor', 'poisoning factor'),
                                    ('fraction of sites', 'f = sum(sigma C) / Gamma', 'TOF / TOF(reference)')):
            ax.set_xscale('log')
            ax.set_xlabel(f'impurity level  p(impurity) / p({probe})')
            P.style_axis(ax, ylabel=ylabel)
            P.keep_title(ax, f"{imp}, {pt['T']:.0f} K, {pt['P'] / 1e5:g} bar: {desc}")
    fig.suptitle('Site blocking: does the poisoning come from lost sites?', x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig
