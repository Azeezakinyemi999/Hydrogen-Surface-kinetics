"""Operating conditions (documentation: Stage 3).

Temperature and pressure are varied separately, never together:
  temperature series : every T in `temperatures`, at `reference_pressure`
  pressure series    : every P in `pressures` (plus reference_pressure, so the series
                       meet), at each T in `pressure_temperatures`
A point (T, P) that is in both is solved once and labeled 'both'.

One row per (basis, case, point). Cases are the probe-gas-only reference, each impurity
alone, and all impurities together. Two bases:
  fixed_probe : the probe-gas partial pressure is P and impurities are ADDED on top at
                impurity_ratio x P (main analysis: probe-gas supply unchanged)
  input_X     : fixed mole fractions at total pressure P (probe gas diluted)

The point's labels (P_ref_bar, series) are not physics: they stay out of the cache
fingerprint, so relabeling a point never forces it to be solved again.
"""
import numpy as np
import pandas as pd

SERIES = ('temperature', 'pressure')
LABEL_COLUMNS = ('P_ref_bar', 'series')


def series_points(cfg, series=None, temperatures=None, pressures=None, pressure_temperatures=None):
    """[{'T': K, 'P': Pa, 'series': 'temperature' | 'pressure' | 'both'}], temperature series
    first (by T), then each pressure series (by T, then P). None = the main settings."""
    series = list(SERIES) if series is None else ([series] if isinstance(series, str) else list(series))
    temps = [float(t) for t in (cfg.temperatures if temperatures is None else temperatures)]
    pres = [float(p) for p in (cfg.pressures if pressures is None else pressures)]
    p_temps = [float(t) for t in (cfg.pressure_temperatures if pressure_temperatures is None else pressure_temperatures)]
    P_ref = float(cfg.reference_pressure)
    tags = {}
    if 'temperature' in series:
        for T in temps:
            tags.setdefault((T, P_ref), set()).add('temperature')
    if 'pressure' in series and pres:
        for T in p_temps:
            for P in sorted(set(pres) | {P_ref}):
                tags.setdefault((T, P), set()).add('pressure')
    return [{'T': T, 'P': P, 'series': 'both' if len(t) == 2 else next(iter(t))} for (T, P), t in tags.items()]


def in_series(df, name):
    """Rows of a results table that belong to one series ('both' counts for either)."""
    return df[df['series'].isin([name, 'both'])]


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


def point_row(model, basis, case, point, pp, **extra):
    """A condition row with the point's labels (see LABEL_COLUMNS)."""
    row = condition_row(model, basis, case, point['T'], pp)
    row.update(extra)
    row.update({'P_ref_bar': point['P'] / 1e5, 'series': point['series']})
    return row


def main_conditions(model):
    """The main conditions table: every case at every series point, on each basis."""
    cfg = model.cfg
    bases = ['fixed_probe'] + (['input_X'] if cfg.include_input_basis else [])
    rows = [point_row(model, basis, case, pt, partial_pressures(model, case, basis, pt['P']))
            for basis in bases for pt in series_points(cfg) for case in case_labels(cfg)]
    df = pd.DataFrame(rows)
    assert all(np.isclose(sum(x.values()), 1.0) for x in df['X'])
    return df


def split_labels(row):
    """(physics, labels): the part of a row that is hashed into the cache key, and the rest."""
    physics = {k: v for k, v in row.items() if k not in LABEL_COLUMNS}
    return physics, {k: row[k] for k in LABEL_COLUMNS if k in row}


def condition_key(row):
    """Readable name for one condition, e.g. 'fixed_probe | CO | 700 K | 1 bar'."""
    return f"{row['basis']} | {row['case']} | {row['T_K']:.0f} K | {row['P_ref_bar']:g} bar"
