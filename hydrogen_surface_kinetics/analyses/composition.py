"""What covers the surface at steady state (documentation: Stage 9).

Outputs: steady_state_coverages.csv (every species with any coverage),
steady_state_surface_species.csv (species above min_listed), surface_composition.png.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import species_label
from . import common


def coverage_table(ctx, results):
    names = ctx.model.species_names
    rows = [pd.Series(res['theta_a'], index=names, name=(res['row']['basis'], res['row']['T_K'], res['row']['case']))
            for res in results.values() if common.converged(res)]
    cov = pd.DataFrame(rows)
    cov.index = pd.MultiIndex.from_tuples(cov.index, names=['basis', 'T_K', 'case'])
    return cov.loc[:, (cov > 0).any()]


def surface_species(ctx, cov, min_listed):
    rows = []
    for (basis, T, case), theta in cov.iterrows():
        for name, c in theta[theta > min_listed].sort_values(ascending=False).items():
            rows.append({'basis': basis, 'T_K': T, 'case': case, 'species': species_label(ctx.model, name),
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
    for (basis, T), part in on_surface.groupby(['basis', 'T_K'], sort=False):
        ctx.log(f'  --- {basis}, {T:.0f} K ---')
        for case, rows in part.groupby('case', sort=False):
            ctx.log(f'    {case:10s} ' + ', '.join(f'{r.species} {r.coverage:.3f}' for r in rows.itertuples()))
    colors = {}
    for basis in cov.index.get_level_values('basis').unique():
        fig, colors[basis] = plot_surface_composition(ctx, cov, basis)
        P.save(fig, ctx.path(f'surface_composition{"" if basis == "fixed_probe" else "_" + basis}.png'), ctx.cfg)
        plt.close(fig)
    ctx.outputs['composition'] = {'coverages': cov, 'on_surface': on_surface,
                                  'colors': colors.get('fixed_probe', {})}
    return ctx.outputs['composition']


def plot_surface_composition(ctx, cov, basis='fixed_probe'):
    model, s = ctx.model, ctx.cfg.composition
    site = model.site
    data = cov.xs(basis, level='basis')
    named, colors = species_palette(data, site, s)
    order = C.case_labels(ctx.cfg)

    temps = sorted(data.index.get_level_values('T_K').unique())
    fig, axes = plt.subplots(1, len(temps), figsize=(4.2 * len(temps), 4.6), sharey=True, squeeze=False)
    axes = axes[0]
    for ax, T in zip(axes, temps):
        part = data.xs(T, level='T_K')
        part = part.reindex([c for c in order if c in part.index])
        x = np.arange(len(part))
        bottom = np.zeros(len(part))
        empty = part[site] if site in part else pd.Series(0.0, index=part.index)
        segments = ([(empty, P.EMPTY_COLOR, 'empty site')]
                    + [(part[n], colors[n], species_label(model, n)) for n in named]
                    + [(part.drop(columns=[site] + named, errors='ignore').sum(axis=1), P.OTHER_COLOR, 'other')])
        for values, color, label in segments:
            ax.bar(x, values, bottom=bottom, width=0.72, color=color, label=label, edgecolor='white', linewidth=1.0)
            bottom += values.to_numpy()
        ax.set_xticks(x, list(part.index), rotation=35, ha='right', fontsize=9)
        ax.set_ylim(0, 1.0)
        P.style_axis(ax, f'{T:.0f} K')
    axes[0].set_ylabel('Fraction of surface sites', fontsize=10, color=P.INK)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles[::-1], labels[::-1], loc='center left', bbox_to_anchor=(1.0, 0.5),
               frameon=False, fontsize=9, title='Surface species', title_fontsize=9)
    fig.suptitle(f'What covers the surface at steady state ({basis})', x=0.01, ha='left', fontsize=11)
    fig.tight_layout()
    return fig, colors
