"""Reaction orders n_j = d ln(TOF) / d ln(p_j) by central differences (documentation:
Stage 11c). p_j is scaled by (1 +/- delta) with every other partial pressure fixed.
Output: reaction_orders.csv.
"""
import numpy as np
import pandas as pd

from .. import conditions as C
from . import common


def order_conditions(ctx):
    model, cfg, s = ctx.model, ctx.cfg, ctx.cfg.orders
    ratio = common.setting(ctx, 'orders', 'ratio', cfg.impurity_ratio)
    rows = []
    for T in common.temperatures(ctx, 'orders'):
        for case in [cfg.reference_label] + common.impurity_labels(ctx, 'orders'):
            base = C.partial_pressures(model, case, 'fixed_probe', common.probe_pressure(ctx, T), ratio)
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
                rows.append(common.fixed_probe_row(ctx, f'order {point}', case, T, pressures=pp))
    return pd.DataFrame(rows)


def run(ctx):
    cfg, s = ctx.cfg, ctx.cfg.orders
    ctx.log('\n== orders: reaction orders ==')
    results = ctx.run_conditions(order_conditions(ctx), prefix='')
    tof, agree = {}, {}
    for res in results.values():
        r = res['row']
        k = (r['T_K'], r['case'], r['basis'].split()[-1])
        tof[k] = res['info_a']['tof_fwd_1_per_s'] if common.converged(res) else np.nan
        agree[k] = res.get('starts_agree')
    step = np.log(1 + s['delta']) - np.log(1 - s['delta'])
    rows = []
    for T in common.temperatures(ctx, 'orders'):
        for case in [cfg.reference_label] + common.impurity_labels(ctx, 'orders'):
            out = {'T_K': T, 'case': case, 'tof_fwd_base': tof.get((T, case, 'base'), np.nan)}
            for var in ('probe', 'impurity'):
                up, dn = tof.get((T, case, f'{var}+'), np.nan), tof.get((T, case, f'{var}-'), np.nan)
                out[f'order_{var}'] = (np.log(up) - np.log(dn)) / step
            pts = [k for k in agree if k[0] == T and k[1] == case]
            out['all_starts_agree'] = all(agree[k] is not None and bool(agree[k]) for k in pts)
            rows.append(out)
    df = pd.DataFrame(rows)
    ctx.save_table(df, 'reaction_orders.csv')
    ctx.log('  reaction orders of the forward rate (n = d ln TOF / d ln p):\n'
            + df.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
    ctx.outputs['orders'] = {'results': results, 'orders': df}
    return ctx.outputs['orders']
