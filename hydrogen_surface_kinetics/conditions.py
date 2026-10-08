"""Operating conditions (documentation: Stage 2).

One row per (basis, case, temperature). Cases are the probe-gas-only reference, each
impurity alone, and all impurities together. Two bases:
  fixed_probe : the probe-gas partial pressure is held at the reference pressure and
                impurities are ADDED on top (main analysis: probe-gas supply unchanged)
  input_X     : fixed mole fractions at the reference total pressure (probe gas diluted)
"""
import numpy as np
import pandas as pd


def reference_pressure(cfg, T):
    """Pa: probe-gas partial pressure (fixed_probe) or total pressure (input_X) at T."""
    return cfg.reference_pressure if cfg.pressure_mode == 'fixed' else cfg.input_pressures[float(T)]


def case_labels(cfg):
    labels = [cfg.reference_label] + list(cfg.impurities)
    if cfg.mixture_label and len(cfg.impurities) > 1:
        labels.append(cfg.mixture_label)
    return labels


def impurities_in(cfg, case):
    if case == cfg.reference_label:
        return []
    if case == cfg.mixture_label:
        return list(cfg.impurities)
    return [case]


def partial_pressures(model, case, basis, P_ref, ratio=None):
    """{gas species: Pa} for one case. `ratio` (fixed_probe only) overrides cfg.impurity_ratio."""
    cfg = model.cfg
    imps = impurities_in(cfg, case)
    if basis == 'fixed_probe':
        r = cfg.impurity_ratio if ratio is None else ratio
        if case == cfg.mixture_label and cfg.mixture_mode == 'split':
            r = r / len(imps)
        pp = {model.probe_gas: P_ref}
        for label in imps:
            pp[model.impurities[label]] = r * P_ref
        return pp
    if basis == 'input_X':
        ib = cfg.input_basis
        if not imps:
            X = {model.probe_gas: 1.0}
        elif case == cfg.mixture_label:
            X = {model.probe_gas: ib['mixture']['probe'],
                 **{model.impurities[label]: ib['mixture']['each_impurity'] for label in imps}}
        else:
            X = {model.probe_gas: ib['single']['probe'], model.impurities[case]: ib['single']['impurity']}
        total = sum(X.values())
        return {sp: x / total * P_ref for sp, x in X.items()}
    raise ValueError(f'unknown basis {basis!r}')


def condition_row(model, basis, case, T, pp):
    P_total = sum(pp.values())
    P_probe = pp[model.probe_gas]
    return {'basis': basis, 'case': case, 'T_K': float(T),
            'P_total_bar': P_total / 1e5,
            'P_probe_bar': P_probe / 1e5,
            'P_impurity_bar': (P_total - P_probe) / 1e5,
            'X': {sp: p / P_total for sp, p in pp.items()}}


def main_conditions(model):
    """The main conditions table: every case at every temperature, on each basis."""
    cfg = model.cfg
    bases = ['fixed_probe'] + (['input_X'] if cfg.include_input_basis else [])
    rows = [condition_row(model, basis, case, T, partial_pressures(model, case, basis, reference_pressure(cfg, T)))
            for basis in bases for T in cfg.temperatures for case in case_labels(cfg)]
    df = pd.DataFrame(rows)
    assert all(np.isclose(sum(x.values()), 1.0) for x in df['X'])
    return df


def condition_key(row):
    """Readable name for one condition, e.g. 'fixed_probe | CO | 700 K'."""
    return f"{row['basis']} | {row['case']} | {row['T_K']:.0f} K"
