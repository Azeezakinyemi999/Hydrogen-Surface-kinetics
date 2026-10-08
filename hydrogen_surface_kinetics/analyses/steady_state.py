"""Main steady states: rates, coverages and poisoning factors (documentation: Stages 6-8).

Outputs: steady_state_results.csv, rates_summary.png (and basis_comparison.csv with
include_input_basis = True).
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
    return df.sort_values(['basis', 'T_K', 'case']).reset_index(drop=True)


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
        ctx.log(bad[['basis', 'case', 'T_K', 'converged', 'starts_agree', 'blocking_check']].to_string(index=False))
    show = ['basis', 'T_K', 'case', 'tof_fwd', 'tof_net', 'tof_probe_consumed', 'theta_free', 'theta_probe',
            'theta_molec', 'site_count_factor', 'poisoning_factor', 'slowdown_x']
    ctx.log(summary[[c for c in show if c in summary]].to_string(index=False, float_format=lambda v: f'{v:.4g}'))

    if {'fixed_probe', 'input_X'} <= set(summary['basis']):
        cmp = summary.pivot_table(index=['T_K', 'case'], columns='basis', values='tof_fwd', observed=True)
        cmp['ratio_input_over_fixed'] = cmp['input_X'] / cmp['fixed_probe']
        ctx.save_table(cmp.reset_index(), 'basis_comparison.csv')

    fig, rate_data = plot_rates(ctx, summary)
    P.save(fig, ctx.path('rates_summary.png'), cfg)
    plt.close(fig)
    slow = rate_data.pivot_table(index='case', columns='T_K', values='slowdown_x', observed=True)
    ctx.log('  times slower than the reference:\n' + slow.map(lambda v: f'{v:.3g}').to_string())
    return ctx.outputs['steady_state']


def plot_rates(ctx, summary, basis='fixed_probe'):
    cfg, model = ctx.cfg, ctx.model
    probe = probe_label(model)
    df = summary[(summary['basis'] == basis) & summary['converged']].copy()
    ref = df[df['case'] == cfg.reference_label].set_index('T_K')['site_count_factor']
    n = model.probe_sites
    df['slowdown_fixed_sites'] = df['slowdown_x'] * (df['site_count_factor'] / df['T_K'].map(ref)) ** n

    cases = [c for c in C.case_labels(cfg) if c in set(df['case'])]
    x_of = {c: i for i, c in enumerate(cases)}
    temps = sorted(df['T_K'].unique())
    t_colors = dict(zip(temps, P.SERIES_COLORS))
    t_markers = dict(zip(temps, P.MARKERS))
    offset = {T: (i - (len(temps) - 1) / 2) * 0.22 for i, T in enumerate(temps)}

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    panels = [('tof_fwd', f'a) {probe} dissociation rate (forward)', 'TOF  (events per site per s)', True),
              ('slowdown_x', f'b) Slowdown vs. {cfg.reference_label}', f'times slower than {cfg.reference_label}', True),
              ('theta_free', 'c) Empty sites', 'fraction of sites empty', False)]
    for ax, (col, title, ylabel, log) in zip(axes, panels):
        for T in temps:
            part = df[df['T_K'] == T]
            xs = [x_of[c] + offset[T] for c in part['case']]
            ax.scatter(xs, part[col], s=46, marker=t_markers[T], color=t_colors[T],
                       edgecolor='white', linewidth=0.8, zorder=3, label=f'{T:.0f} K')
            if col == 'slowdown_x':
                ax.scatter(xs, part['slowdown_fixed_sites'], s=46, marker=t_markers[T], facecolor='none',
                           edgecolor=t_colors[T], linewidth=1.2, zorder=2)
        if log:
            ax.set_yscale('log')
        if col == 'slowdown_x':
            ax.axhline(1.0, color=P.MUTED, linewidth=1.0, linestyle='--', zorder=1)
            ax.text(len(cases) - 0.5, 1.0, ' no effect', va='bottom', ha='right', fontsize=8, color=P.MUTED)
        if col == 'theta_free':
            ax.set_ylim(0, max(0.4, df[col].max() * 1.1))
        ax.set_xticks(range(len(cases)), cases, rotation=35, ha='right', fontsize=9)
        ax.set_xlim(-0.6, len(cases) - 0.4)
        P.style_axis(ax, title, ylabel)
    handles, labels = axes[0].get_legend_handles_labels()
    handles.append(plt.Line2D([], [], marker='o', linestyle='', markerfacecolor='none',
                              markeredgecolor=P.MUTED, markersize=7))
    labels.append('panel b, open: if the site count were fixed')
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(labels),
               frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    return fig, df
