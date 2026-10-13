"""What covers the surface at steady state (documentation: Stage 9).

Outputs: steady_state_coverages.csv (every species with any coverage),
steady_state_surface_species.csv (species above min_listed), surface_composition_vs_T.png
and, with a pressure series, surface_composition_vs_P.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import species_label
from . import common


INDEX = ['basis', 'series', 'T_K', 'P_ref_bar', 'case']


def coverage_table(ctx, results):
    """Coverages of every converged condition: index (basis, series, T_K, P_ref_bar, case)."""
    names = ctx.model.species_names
    rows = [pd.Series(res['theta_a'], index=names, name=tuple(res['row'].get(k) for k in INDEX))
            for res in results.values() if common.converged(res)]
    cov = pd.DataFrame(rows)
    cov.index = pd.MultiIndex.from_tuples(cov.index, names=INDEX)
    return cov.loc[:, (cov > 0).any()]


def surface_species(ctx, cov, min_listed):
    rows = []
    for idx, theta in cov.iterrows():
        for name, c in theta[theta > min_listed].sort_values(ascending=False).items():
            rows.append({**dict(zip(INDEX, idx)), 'species': species_label(ctx.model, name),
                         'mechanism_name': name, 'coverage': c})
    return pd.DataFrame(rows)


def species_palette(cov, site, s, base=None):
    """Species that get their own color: those reaching s['min_plotted'] anywhere in `cov`,
    largest peak first, at most s['max_named']. Colors follow the species, not the rank:
    a species already in `base` ({species: color}, e.g. from the surface-composition figure)
    keeps that color, and the others take the unused palette slots in fixed order.
    Returns (named species, {species: color})."""
    base = dict(base or {})
    peak = cov.drop(columns=[site], errors='ignore').max()
    own = list(peak[peak >= s['min_plotted']].sort_values(ascending=False).index)
    free = [c for c in P.SERIES_COLORS if c not in base.values()]
    named, colors = [], {}
    for sp in own:
        if len(named) == s['max_named']:
            break
        if sp in base:
            colors[sp] = base[sp]
        elif free:
            colors[sp] = free.pop(0)
        else:
            continue                                   # out of slots: folds into "other"
        named.append(sp)
    return named, colors


def run(ctx):
    s = ctx.cfg.composition
    ctx.log('\n== composition ==')
    main = ctx.main()
    cov = coverage_table(ctx, main['results'])
    on_surface = surface_species(ctx, cov, s['min_listed'])
    ctx.save_table(cov, 'steady_state_coverages.csv', index=True)
    ctx.save_table(on_surface, 'steady_state_surface_species.csv')
    for (basis, T, Pb), part in on_surface.groupby(['basis', 'T_K', 'P_ref_bar'], sort=False):
        ctx.log(f'  --- {basis}, {T:.0f} K, {Pb:g} bar ---')
        for case, rows in part.groupby('case', sort=False):
            ctx.log(f'    {case:10s} ' + ', '.join(f'{r.species} {r.coverage:.3f}' for r in rows.itertuples()))
    colors = {}
    for basis in cov.index.get_level_values('basis').unique():
        data = cov.xs(basis, level='basis')
        named, colors[basis] = species_palette(data, ctx.model.site, s)      # one palette for both series
        for series in C.SERIES:
            fig = plot_surface_composition(ctx, data, series, named, colors[basis], basis)
            if fig is None:
                continue
            stem = 'surface_composition' + ('' if basis == 'fixed_probe' else '_' + basis)
            P.save(fig, ctx.path(common.series_file(stem, series)), ctx.cfg)
            plt.close(fig)
    ctx.outputs['composition'] = {'coverages': cov, 'on_surface': on_surface,
                                  'colors': colors.get('fixed_probe', {})}
    return ctx.outputs['composition']


def plot_surface_composition(ctx, data, series, named, colors, basis='fixed_probe'):
    """Stacked bars per case, one panel per point of the series (temperature series: a panel
    per T; pressure series: a row per pressure-series temperature, a panel per P)."""
    model = ctx.model
    site = model.site
    flat = data.reset_index()
    groups = common.series_groups(flat, series)
    if not groups:
        return None
    order = C.case_labels(ctx.cfg)
    ncol = max(g[1][g[2]].nunique() for g in groups)
    fig, axes = plt.subplots(len(groups), ncol, figsize=(4.2 * ncol, 4.6 * len(groups)), sharey=True, squeeze=False)
    for row_axes, (title, gdf, vary, label_of) in zip(axes, groups):
        values = sorted(gdf[vary].unique())
        for ax, v in zip(row_axes, values):
            rows = gdf[gdf[vary] == v].set_index('case')
            part = rows.reindex([c for c in order if c in rows.index])[data.columns]
            x = np.arange(len(part))
            bottom = np.zeros(len(part))
            empty = part[site] if site in part else pd.Series(0.0, index=part.index)
            segments = ([(empty, P.EMPTY_COLOR, 'empty site')]
                        + [(part[n], colors[n], species_label(model, n)) for n in named if n in part]
                        + [(part.drop(columns=[site] + named, errors='ignore').sum(axis=1), P.OTHER_COLOR, 'other')])
            for values_, color, label in segments:
                ax.bar(x, values_, bottom=bottom, width=0.72, color=color, label=label, edgecolor='white',
                       linewidth=1.0)
                bottom += values_.to_numpy()
            ax.set_xticks(x, list(part.index), rotation=35, ha='right', fontsize=9)
            ax.set_ylim(0, 1.0)
            T = rows['T_K'].iloc[0]
            Pb = rows['P_ref_bar'].iloc[0]
            P.style_axis(ax)
            P.keep_title(ax, f'{T:.0f} K, {Pb:g} bar')
        for ax in row_axes[len(values):]:
            ax.set_visible(False)
        row_axes[0].set_ylabel('Fraction of surface sites', fontsize=10, color=P.INK)
    handles = ([plt.Rectangle((0, 0), 1, 1, color=P.EMPTY_COLOR)] + [plt.Rectangle((0, 0), 1, 1, color=colors[n])
               for n in named] + [plt.Rectangle((0, 0), 1, 1, color=P.OTHER_COLOR)])
    labels = ['empty site'] + [species_label(model, n) for n in named] + ['other']
    fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=9, title='Surface species', title_fontsize=9)
    fig.suptitle(f'What covers the surface at steady state ({basis}, {series} series)', x=0.01, ha='left', fontsize=11)
    fig.tight_layout()
    return fig


def plot_coverage_bars(ctx, cov, cases, xs, xlabel, title, log_x=False):
    """Stacked coverages against one varied quantity: one bar per solved x (cov indexed by
    (case, x)), one panel per case. Species keep the colors of the surface-composition figure
    where it ran in this session. log_x: bars evenly spaced in log10(x) (pressures)."""
    model, s = ctx.model, ctx.cfg.composition
    site = model.site
    base = ctx.outputs.get('composition', {}).get('colors')
    named, colors = species_palette(cov, site, s, base)
    xs = sorted(float(x) for x in xs)
    pos = (lambda x: np.log10(np.asarray(x, dtype=float))) if log_x else (lambda x: np.asarray(x, dtype=float))
    steps = np.diff(pos(xs))
    width = 0.8 * (np.min(steps) if len(steps) else (0.5 if log_x else 50.0))
    ncol = min(3, len(cases))
    nrow = int(np.ceil(len(cases) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.4 * ncol, 3.3 * nrow), sharex=True, sharey=True, squeeze=False)
    layer_names = [site] + named + ['other']
    color_of = lambda n: colors.get(n, P.EMPTY_COLOR if n == site else P.OTHER_COLOR)   # noqa: E731
    for ax, case in zip(axes.flat, cases):
        part = cov.xs(case, level='case')
        part = part[part.index.isin(xs)].sort_index()                 # solved values only: gaps elsewhere
        stack = {site: part[site] if site in part else pd.Series(0.0, index=part.index),
                 'other': part.drop(columns=[site] + named, errors='ignore').sum(axis=1)}
        stack.update({n: part[n] for n in named if n in part})
        bottom = np.zeros(len(part))
        for n in layer_names:
            if n in stack:
                vals = stack[n].to_numpy(dtype=float)
                ax.bar(pos(part.index), vals, bottom=bottom, width=width, color=color_of(n), edgecolor='white',
                       linewidth=0.8)
                bottom = bottom + vals
        ax.set_ylim(0, 1)
        ax.set_xlim(pos(xs[0]) - width, pos(xs[-1]) + width)
        if log_x:                                                     # ticks at whole decades, labeled in x
            decades = np.arange(np.floor(pos(xs[0])), np.ceil(pos(xs[-1])) + 1)
            ax.set_xticks(decades, [f'{10 ** d:g}' for d in decades])
        P.style_axis(ax)
        P.keep_title(ax, case)
    for ax in axes.flat[len(cases):]:
        ax.set_visible(False)
    for ax in axes[:, 0]:
        ax.set_ylabel('fraction of sites', fontsize=10, color=P.INK)
    for j in range(ncol):                                             # x label under the lowest panel of each column
        col = [axes[i][j] for i in range(nrow) if axes[i][j].get_visible()]
        col[-1].set_xlabel(xlabel, fontsize=10, color=P.INK)
        col[-1].xaxis.set_tick_params(labelbottom=True)
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_of(n)) for n in layer_names]   # stack order, top first
    labels = ['empty site'] + [species_label(model, n) for n in named] + ['other']
    fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=9, title='Surface species', title_fontsize=9)
    fig.suptitle(title, x=0.01, ha='left', fontsize=11)
    fig.tight_layout()
    return fig
