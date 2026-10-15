"""Reaction orders n_j = d ln(TOF) / d ln(p_j) by central differences (documentation:
Stage 11.3), at every point of the temperature and pressure series. p_j is scaled by
(1 +/- delta) with every other partial pressure fixed.
Outputs: reaction_orders.csv, reaction_orders_vs_T.png and (with a pressure series)
reaction_orders_vs_P.png: one panel per point.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import probe_label
from . import common


def order_conditions(ctx):
    model, cfg, s = ctx.model, ctx.cfg, ctx.cfg.orders
    ratio = common.setting(ctx, 'orders', 'ratio', cfg.impurity_ratio)
    rows = []
    for pt in common.points(ctx, 'orders'):
        for case in [cfg.reference_label] + common.impurity_labels(ctx, 'orders'):
            base = C.partial_pressures(model, case, 'fixed_probe', pt['P'], ratio)
            points = [('base', base)]
            variables = [('probe', model.probe_gas)]
            if case != cfg.reference_label:
                variables.append(('impurity', model.impurities[case]))
            for var, sp in variables:
                for sign, factor in (('+', 1 + s['delta']), ('-', 1 - s['delta'])):
                    pp = dict(base)
                    pp[sp] = base[sp] * factor
                    points.append((f'{var}{sign}', pp))
            for point, pp in points:
                rows.append(common.fixed_probe_row(ctx, f'order {point}', case, pt, pressures=pp))
    return pd.DataFrame(rows)


def run(ctx):
    cfg, s = ctx.cfg, ctx.cfg.orders
    ctx.log('\n== orders: reaction orders ==')
    results = ctx.run_conditions(order_conditions(ctx), prefix='')
    tof, agree = {}, {}
    for res in results.values():
        r = res['row']
        k = (r['T_K'], r['P_ref_bar'], r['case'], r['basis'].split()[-1])
        tof[k] = res['info_a']['tof_fwd_1_per_s'] if common.converged(res) else np.nan
        agree[k] = res.get('starts_agree')
    step = np.log(1 + s['delta']) - np.log(1 - s['delta'])
    rows = []
    for pt in common.points(ctx, 'orders'):
        T, Pb = pt['T'], pt['P'] / 1e5
        for case in [cfg.reference_label] + common.impurity_labels(ctx, 'orders'):
            out = {'series': pt['series'], 'T_K': T, 'P_ref_bar': Pb, 'case': case,
                   'tof_fwd_base': tof.get((T, Pb, case, 'base'), np.nan)}
            for var in ('probe', 'impurity'):
                up, dn = tof.get((T, Pb, case, f'{var}+'), np.nan), tof.get((T, Pb, case, f'{var}-'), np.nan)
                out[f'order_{var}'] = (np.log(up) - np.log(dn)) / step
            keys = [k for k in agree if k[:3] == (T, Pb, case)]
            out['all_starts_agree'] = all(agree[k] is not None and bool(agree[k]) for k in keys)
            rows.append(out)
    df = pd.DataFrame(rows)
    ctx.save_table(df, 'reaction_orders.csv')
    ctx.log('  reaction orders of the forward rate (n = d ln TOF / d ln p):\n'
            + df.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    for series in C.SERIES:
        P.save_parts(plot_orders(ctx, df, series), ctx.path(common.series_file('reaction_orders', series)), cfg)
    ctx.outputs['orders'] = {'results': results, 'orders': df}
    return ctx.outputs['orders']


def plot_orders(ctx, orders, series='temperature'):
    """Reaction orders of the forward probe rate per case, one panel per point of the series
    (at most 3 x 3 per figure): circle = order in the probe gas, square = order in the
    impurity, labeled where |n| >= 0.1 and, for a case, where it differs from the reference
    case's order by 0.05 or more. Returns the figures."""
    part = C.in_series(orders, series)
    pts = list(part[['T_K', 'P_ref_bar']].drop_duplicates().sort_values(['T_K', 'P_ref_bar']).itertuples(index=False))
    vals = part[['order_probe', 'order_impurity']].to_numpy(dtype=float)
    lo, hi = (np.nanmin(vals), np.nanmax(vals)) if np.isfinite(vals).any() else (0.0, 0.0)
    ylim = (min(lo, -0.5) - 0.5, max(hi, 0.5) + 0.5)                 # the same scale in every panel
    return [_orders_page(ctx, part, page, ylim) for page in P.pages(pts, 3 * P.MAX_ROWS)]


def _orders_page(ctx, orders, pts, ylim):
    probe = probe_label(ctx.model)
    nrow, ncol = P.grid_shape(len(pts), max_cols=3)
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.4 * nrow + 0.5), squeeze=False)
    markers = (('order_probe', -0.12, P.SERIES_COLORS[0], 'o', f'order in {probe}', 'right'),
               ('order_impurity', 0.12, P.SERIES_COLORS[1], 's', 'order in the impurity', 'left'))
    for k, (ax, (T, Pb)) in enumerate(zip(axes.flat, pts)):
        o = orders[(orders['T_K'] == T) & np.isclose(orders['P_ref_bar'], Pb)].set_index('case')
        cs = [c for c in C.case_labels(ctx.cfg) if c in o.index]
        x = np.arange(len(cs))
        ref = ctx.cfg.reference_label
        for col, dx, color, marker, lab, ha in markers:
            v_all = o.loc[cs, col].to_numpy(dtype=float)
            v_ref = o.loc[ref, col] if ref in o.index else np.nan
            ax.plot(x + dx, v_all, ls='', marker=marker, ms=7, color=color, mec='white', mew=0.8, zorder=3, label=lab)
            for xi, v, case in zip(x + dx, v_all, cs):
                same = case != ref and np.isfinite(v_ref) and abs(v - v_ref) < 0.05     # as the reference: unlabeled
                if np.isfinite(v) and abs(v) >= 0.1 and not same:
                    ax.annotate(P.signed(v), (xi, v), xytext=(-3 if ha == 'right' else 3, 6 if v > 0 else -6),
                                textcoords='offset points', ha=ha, va='bottom' if v > 0 else 'top', fontsize=8,
                                color=P.SECONDARY_INK)
        ax.axhline(0, color=P.MUTED, lw=0.8, zorder=1)
        ax.set_xticks(x, cs, rotation=35, ha='right', fontsize=9)
        ax.set_xlim(-0.6, len(cs) - 0.4)
        ax.set_ylim(*ylim)
        P.style_axis(ax, ylabel='order  d ln TOF / d ln p' if k % ncol == 0 else '')
        P.keep_title(ax, f'{T:.0f} K, {Pb:g} bar')
    for ax in axes.flat[len(pts):]:
        ax.set_visible(False)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False, fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    return fig
