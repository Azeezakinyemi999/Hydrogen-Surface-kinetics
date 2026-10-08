"""Compare this run's output CSVs with those of the original notebook
(Surface-Impurity-Effect/analysis/steady_state_Ni.ipynb), which stays the reference
until it is retired. Run after `bash run.sh runs/ni_first_try`:

    python runs/ni_first_try/compare_with_notebook.py

Prints the worst difference per column (relative, or absolute for orders and DRC/TRC).
Expected (October 2026): <= 2e-7 everywhere except near-zero quantities (net rates at
equilibrium, tiny coverages), the DRC reactions screened for NH3 (all with X ~ 0), and
transient half-times (a few %, see docs Stage 11.6).
"""
import os
import numpy as np
import pandas as pd

OLD = '/Users/akinyemi.az/Desktop/PhD_Folder/research/MHI/Surface-Impurity-Effect/rmg_runs/Ni/first_try_saved/steady_state'
NEW = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'output')
CASE = {'pure H2': 'H2 only'}
COLS = {'tof_diss_fwd': 'tof_fwd', 'tof_diss_net': 'tof_net', 'tof_diss_rev': 'tof_rev', 'tof_H2_consumed': 'tof_probe_consumed', 'theta_H': 'theta_probe', 'theta_H2X': 'theta_molec', 'order_H2': 'order_probe',
        'X_h2_consumed': 'X_probe_consumed', 'tof_without_H2X': 'tof_switched_off', 'tof_with_H2X': 'tof_normal'}


def load(name_old, name_new=None):
    o = pd.read_csv(os.path.join(OLD, name_old)).rename(columns=COLS)
    n = pd.read_csv(os.path.join(NEW, name_new or name_old))
    if 'case' in o:
        o['case'] = o['case'].replace(CASE)
    if 'basis' in o:
        o['basis'] = o['basis'].replace({'fixed_PH2': 'fixed_probe'})
    return o, n


def compare(title, name_old, keys, values, name_new=None, rtol_log=False, query=None):
    o, n = load(name_old, name_new)
    if query:
        n = n.query(query)
    for k in keys:
        if k in ('T_K', 'ratio'):
            o[k], n[k] = o[k].astype(float).round(12), n[k].astype(float).round(12)
        if k == 'item':
            o[k], n[k] = o[k].astype(str), n[k].astype(str)
    m = o.merge(n, on=keys, suffixes=('_old', '_new'), how='outer', indicator=True)
    miss = m[m['_merge'] != 'both']
    print(f'\n=== {title}: {len(m)} rows, {len(miss)} unmatched')
    if len(miss):
        print(miss[keys + ['_merge']].head(10).to_string(index=False))
    m = m[m['_merge'] == 'both']
    for v in values:
        a, b = m[f'{v}_old'].astype(float), m[f'{v}_new'].astype(float)
        both_nan = a.isna() & b.isna()
        one_nan = a.isna() ^ b.isna()
        d = (b - a).abs() / np.maximum(a.abs(), 1e-300) if not rtol_log else (b - a).abs()
        d = d[~both_nan & ~one_nan]
        worst = d.idxmax() if len(d) else None
        where = '' if worst is None else ' at ' + ', '.join(str(m.loc[worst, k]) for k in keys)
        kind = 'abs' if rtol_log else 'rel'
        print(f'  {v:22s} worst {kind} diff {d.max() if len(d) else np.nan:9.2e}{where}   NaN mismatches {int(one_nan.sum())}')


compare('main steady states', 'steady_state_results.csv', ['basis', 'case', 'T_K'],
        ['tof_fwd', 'tof_net', 'theta_free', 'theta_probe', 'theta_molec', 'poisoning_factor', 'blocking_check'])
compare('poisoning curves', 'poisoning_curves.csv', ['T_K', 'case', 'ratio'],
        ['tof_fwd', 'theta_free', 'theta_impurity', 'poisoning_factor'])
compare('temperature sweep', 'temperature_sweep.csv', ['T_K', 'case'],
        ['tof_fwd', 'theta_free', 'poisoning_factor', 'E_app_kJ_per_mol'])
compare('reaction orders', 'reaction_orders.csv', ['T_K', 'case'], ['tof_fwd_base', 'order_probe', 'order_impurity'],
        rtol_log=True)
compare('sensitivity (DRC/TRC)', 'sensitivity.csv', ['kind', 'T_K', 'case', 'item'], ['X_tof_fwd', 'X_probe_consumed'],
        rtol_log=True)
compare('poisoning over time', 'poisoning_over_time.csv', ['T_K', 'case', 'ratio'],
        ['tof_t0', 'end_over_t0', 'min_over_t0', 't_half_s', 't_tenth_s'])
compare('binding check', 'mechanism_binding_check.csv', ['adsorption'], ['dH_298_mech', 'dG_700', 'dG_1000'])
compare('switch-off (no H2X)', 'mechanism_no_H2X.csv', ['case', 'T_K'], ['tof_normal', 'tof_switched_off'],
        name_new='mechanism_switch_off.csv')
