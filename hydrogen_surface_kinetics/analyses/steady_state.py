"""Main steady states: rates, coverages and poisoning factors (documentation: Stages 6-9).

Outputs: steady_state_results.csv; rates_summary_vs_T.png and, with a pressure series,
rates_summary_vs_P.png (basis_comparison.csv with include_input_basis = True).
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from .. import solver as S
from ..labels import probe_label
from . import common


def build_summary(ctx, results):
    cfg, model = ctx.cfg, ctx.model
    df = common.results_frame(ctx, results)
    df = common.add_poisoning(df, cfg.reference_label, model.probe_sites)
    df['case'] = pd.Categorical(df['case'], categories=C.case_labels(cfg), ordered=True)
    return df.sort_values(['basis', 'T_K', 'P_ref_bar', 'case']).reset_index(drop=True)


def run(ctx):
    cfg, model = ctx.cfg, ctx.model
    ctx.log('\n== steady_state: main conditions ==')
    conds = C.main_conditions(model)
    results = ctx.run_conditions(conds, prefix='')
    status = S.status_table(results)
    summary = build_summary(ctx, results)
    ctx.outputs['steady_state'] = {'conditions': conds, 'results': results, 'status': status, 'summary': summary}
    ctx.save_table(summary, 'steady_state_results.csv')

    bad = summary[~summary['converged'] | (summary['starts_agree'] != True)        # noqa: E712
                  | ((summary['blocking_check'] - 1).abs() > 1e-6)]
    ctx.log(f'  {len(summary) - len(bad)} / {len(summary)} conditions fully OK')
    if len(bad):
        ctx.log(bad[['basis', 'case', 'T_K', 'P_ref_bar', 'converged', 'starts_agree', 'blocking_check']]
                .to_string(index=False))
    show = ['basis', 'series', 'T_K', 'P_ref_bar', 'case', 'tof_fwd', 'tof_net', 'tof_probe_consumed', 'theta_free',
            'theta_probe', 'theta_molec', 'site_count_factor', 'poisoning_factor', 'slowdown_x']
    ctx.log(summary[[c for c in show if c in summary]].to_string(index=False, float_format=lambda v: f'{v:.4g}'))

    if {'fixed_probe', 'input_X'} <= set(summary['basis']):
        cmp = summary.pivot_table(index=['T_K', 'P_ref_bar', 'case'], columns='basis', values='tof_fwd', observed=True)
        cmp['ratio_input_over_fixed'] = cmp['input_X'] / cmp['fixed_probe']
        ctx.save_table(cmp.reset_index(), 'basis_comparison.csv')

    for series in C.SERIES:
        fig = plot_rates(ctx, summary, series)
        if fig is not None:
            P.save(fig, ctx.path(common.series_file('rates_summary', series)), cfg)
            plt.close(fig)
    ok = summary[(summary['basis'] == 'fixed_probe') & summary['converged']].copy()
    ok['point'] = [common.point_name(r) for r in ok.to_dict('records')]
    slow = ok.pivot_table(index='case', columns='point', values='slowdown_x', observed=True, sort=False)
    ctx.log('  times slower than the reference:\n' + slow.map(lambda v: f'{v:.3g}').to_string())
    return ctx.outputs['steady_state']


def plot_rates(ctx, summary, series='temperature', basis='fixed_probe'):
    """Forward TOF, slowdown and empty sites per case. One row of panels per group of the
    series (temperature series: one row, a marker per T; pressure series: a row per
    pressure-series temperature, a marker per P)."""
    cfg, model = ctx.cfg, ctx.model
    probe = probe_label(model)
    df = summary[(summary['basis'] == basis) & summary['converged']].copy()
    groups = common.series_groups(df, series)
    if not groups:
        return None
    ref = df[df['case'] == cfg.reference_label].set_index(['T_K', 'P_ref_bar'])['site_count_factor']
    n = model.probe_sites
    f_ref = [ref.get((T, Pb), np.nan) for T, Pb in zip(df['T_K'], df['P_ref_bar'])]
    df['slowdown_fixed_sites'] = df['slowdown_x'] * (df['site_count_factor'] / np.array(f_ref)) ** n

    cases = [c for c in C.case_labels(cfg) if c in set(df['case'])]
    x_of = {c: i for i, c in enumerate(cases)}
    fig, axes = plt.subplots(len(groups), 3, figsize=(15, 4.6 * len(groups)), squeeze=False)
    panels = [('tof_fwd', f'{probe} dissociation rate (forward)', 'TOF  (events per site per s)', True),
              ('slowdown_x', f'Slowdown vs. {cfg.reference_label}', f'times slower than {cfg.reference_label}', True),
              ('theta_free', 'Empty sites', 'fraction of sites empty', False)]
    for row_axes, (title, gdf, vary, label_of) in zip(axes, groups):
        gdf = df.loc[gdf.index]
        values = sorted(gdf[vary].unique())
        colors = dict(zip(values, P.SERIES_COLORS))
        markers = dict(zip(values, P.MARKERS))
        offset = {v: (i - (len(values) - 1) / 2) * min(0.22, 0.7 / max(len(values), 1)) for i, v in enumerate(values)}
        for ax, (col, desc, ylabel, log) in zip(row_axes, panels):
            for v in values:
                part = gdf[gdf[vary] == v]
                xs = [x_of[c] + offset[v] for c in part['case']]
                ax.scatter(xs, part[col], s=46, marker=markers[v], color=colors[v],
                           edgecolor='white', linewidth=0.8, zorder=3, label=label_of(v))
                if col == 'slowdown_x':
                    ax.scatter(xs, part['slowdown_fixed_sites'], s=46, marker=markers[v], facecolor='none',
                               edgecolor=colors[v], linewidth=1.2, zorder=2)
            if log:
                ax.set_yscale('log')
            if col == 'slowdown_x':
                ax.axhline(1.0, color=P.MUTED, linewidth=1.0, linestyle='--', zorder=1)
                ax.text(len(cases) - 0.5, 1.0, ' no effect', va='bottom', ha='right', fontsize=8, color=P.MUTED)
            if col == 'theta_free':
                ax.set_ylim(0, max(0.4, gdf[col].max() * 1.1))
            ax.set_xticks(range(len(cases)), cases, rotation=35, ha='right', fontsize=9)
            ax.set_xlim(-0.6, len(cases) - 0.4)
            P.style_axis(ax, ylabel=ylabel)
            P.keep_title(ax, f'{desc}, {title}' if len(groups) > 1 or series == 'pressure' else desc)
        handles, labels = row_axes[0].get_legend_handles_labels()
        if len(groups) > 1:                                    # each row has its own values: legend per row
            row_axes[2].legend(handles, labels, loc='upper right', frameon=False, fontsize=8)
    if len(groups) == 1:
        handles, labels = axes[0][0].get_legend_handles_labels()
        handles.append(plt.Line2D([], [], marker='o', linestyle='', markerfacecolor='none',
                                  markeredgecolor=P.MUTED, markersize=7))
        labels.append('panel b, open: if the site count were fixed')
        fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(labels),
                   frameon=False, fontsize=9)
    fig.suptitle(f'Main steady states: {series} series', x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.93 if len(groups) == 1 else 0.97))
    return fig
