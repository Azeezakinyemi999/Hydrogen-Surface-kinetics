"""Shared pieces for the analyses: the run context, rates at a solved state, condition rows."""
import os

import numpy as np
import pandas as pd

from .. import conditions as C
from .. import solver as S


class Context:
    """What one run shares between analyses: the model, its settings, and each analysis'
    outputs (ctx.outputs['steady_state'] = {'results': ..., 'summary': ...}, ...)."""

    def __init__(self, model, run_overrides=None, log=print):
        self.model, self.cfg = model, model.cfg
        self.run_overrides = {k: v for k, v in (run_overrides or {}).items() if v is not None}
        self.outputs = {}
        self.log = log
        os.makedirs(self.cfg.output_dir, exist_ok=True)

    def path(self, *parts):
        p = os.path.join(self.cfg.output_dir, *parts)
        os.makedirs(os.path.dirname(p), exist_ok=True)
        return p

    def save_table(self, df, name, **kwargs):
        p = self.path(name)
        df.to_csv(p, **{'index': False, **kwargs})
        self.log(f'  saved {os.path.relpath(p, self.cfg.run_dir)}')
        return p

    def run_conditions(self, conditions, prefix, both_starts=True, cache_subdir='cache'):
        return S.run_conditions(self.model, conditions, cache_subdir=cache_subdir, prefix=prefix,
                                both_starts=both_starts, log=self.log, **self.run_overrides)

    def run_jobs(self, func, jobs, cache_subdir, describe=None, extra_context=()):
        from ..runner import run_jobs
        run = self.cfg.section('run', **self.run_overrides)
        return run_jobs(func, jobs, cache_dir=self.path(cache_subdir, ''),
                        context=(self.model.fingerprint_context(), extra_context),
                        workers=run['workers'], timeout=run['job_timeout'], recompute=run['recompute'],
                        describe=describe, log=self.log)

    def main(self):
        """Main steady-state outputs, running that analysis first if needed."""
        if 'steady_state' not in self.outputs:
            from . import steady_state
            steady_state.run(self)
        return self.outputs['steady_state']


# --- settings helpers ------------------------------------------------------------------
def setting(ctx, section, key, fallback):
    v = ctx.cfg[section].get(key)
    return fallback if v is None else v


def points(ctx, section=None):
    """The series points of an analysis: its section's series / temperatures / pressures /
    pressure_temperatures, each None = the main settings (see conditions.series_points)."""
    sec = ctx.cfg[section] if section else {}
    return C.series_points(ctx.cfg, series=sec.get('series'), temperatures=sec.get('temperatures'),
                           pressures=sec.get('pressures'), pressure_temperatures=sec.get('pressure_temperatures'))


def reference_point(ctx, T, series='temperature'):
    return {'T': float(T), 'P': float(ctx.cfg.reference_pressure), 'series': series}


def point_name(pt):
    """'300 K, 1 bar' (pt: a series point, or a row with T_K and P_ref_bar)."""
    T = pt['T'] if 'T' in pt else pt['T_K']
    P = pt['P'] / 1e5 if 'P' in pt else pt['P_ref_bar']
    return f'{T:.0f} K, {P:g} bar'


SERIES_SUFFIX = {'temperature': 'vs_T', 'pressure': 'vs_P'}


def series_file(name, series, ext='png'):
    """'rates_summary_vs_T.png' / 'rates_summary_vs_P.png'."""
    return f'{name}_{SERIES_SUFFIX[series]}.{ext}'


def series_groups(df, series):
    """Split a results table (with T_K, P_ref_bar, series) into figure groups for one series:
    temperature series -> one group at the reference pressure, varying T;
    pressure series    -> one group per pressure-series temperature, varying P.
    Returns [(group title, rows, varying column, value -> label)]."""
    part = C.in_series(df, series)
    if part.empty:
        return []
    if series == 'temperature':
        return [(f'{P:g} bar', g, 'T_K', lambda v: f'{v:.0f} K') for P, g in part.groupby('P_ref_bar', sort=True)]
    return [(f'{T:.0f} K', g, 'P_ref_bar', lambda v: f'{v:g} bar') for T, g in part.groupby('T_K', sort=True)]


def impurity_labels(ctx, section):
    labels = setting(ctx, section, 'impurities', list(ctx.cfg.impurities))
    unknown = [l for l in labels if l not in ctx.cfg.impurities]
    if unknown:
        raise ValueError(f"{section}: unknown impurities {unknown}")
    return list(labels)


def fixed_probe_row(ctx, basis, case, point, impurity_ratio=None, pressures=None, **extra):
    """A fixed_probe condition row at a series point; `pressures` ({species: Pa}) overrides
    the case (e.g. a perturbed point). The row carries the point's labels."""
    pp = pressures or C.partial_pressures(ctx.model, case, 'fixed_probe', point['P'], impurity_ratio)
    return C.point_row(ctx.model, basis, case, point, pp, **extra)


def job(key, row, **fields):
    """A runner job: the row's physics is hashed, its labels go into '_labels'."""
    physics, labels = C.split_labels(row)
    return {'key': key, 'row': physics, '_labels': labels, **fields}


def job_row(res):
    """The full row (physics + labels) of a finished runner job."""
    return {**res['job']['row'], **res['job'].get('_labels', {})}


# --- rates at a solved state -----------------------------------------------------------
def load_state(model, row, theta, multipliers=None, h_shifts=None):
    """Fresh surface with the gas and the solved coverages (no solving)."""
    surf = model.new_surface(multipliers, h_shifts)
    model.gas.TPX = row['T_K'], row['P_total_bar'] * 1e5, row['X']
    surf.TP = row['T_K'], row['P_total_bar'] * 1e5
    surf.coverages = theta
    return surf


def evaluate(model, row, theta):
    """Rates (1/s per site) and coverages at a solved steady state."""
    surf = load_state(model, row, theta)
    G = model.gamma

    def rates(i, sign):
        fwd, rev = surf.forward_rates_of_progress[i], surf.reverse_rates_of_progress[i]
        if sign == -1:
            fwd, rev = rev, fwd
        return fwd / G, rev / G, (fwd - rev) / G

    fwd, rev, net = rates(model.i_probe, model.sign_probe)
    cov = dict(zip(surf.species_names, theta))
    out = {'tof_fwd': fwd, 'tof_rev': rev, 'tof_net': net,
           'tof_probe_consumed': -surf.net_production_rates[model.probe_kin_idx] / G}
    if model.i_molec is not None:
        m_fwd, _, m_net = rates(model.i_molec, model.sign_molec)
        out.update({'tof_molec_fwd': m_fwd, 'tof_molec_net': m_net})
    out['theta_free'] = cov[model.site]
    out['theta_probe'] = cov[model.probe_adsorbate]
    out['theta_molec'] = cov[model.molecular_adsorbate] if model.molecular_adsorbate else 0.0
    out['theta_impurity'] = max(0.0, 1.0 - out['theta_free'] - out['theta_probe'] - out['theta_molec'])
    # Sum C_k = Gamma for species, so multi-site species raise the site count (up to 2x)
    out['site_count_factor'] = (surf.concentrations * model.site_size).sum() / G
    out['top_species'] = top_species(model, theta)
    return out


def role_coverages(model, theta):
    """(empty sites, probe adsorbate, molecular adsorbate, everything else) of a coverage array."""
    names = model.species_names
    free = theta[names.index(model.site)]
    probe = theta[names.index(model.probe_adsorbate)]
    molec = theta[names.index(model.molecular_adsorbate)] if model.molecular_adsorbate else 0.0
    return free, probe, molec, max(0.0, 1.0 - free - probe - molec)


def sweep_cases(ctx, section):
    """Reference + the section's impurities + the all-impurities mixture (when the section
    uses every impurity and there is a mixture case)."""
    cfg = ctx.cfg
    cases = [cfg.reference_label] + impurity_labels(ctx, section)
    if cfg[section].get('impurities') is None and cfg.mixture_label and len(cfg.impurities) > 1:
        cases.append(cfg.mixture_label)
    return cases


def top_species(model, theta, n=3):
    order = np.argsort(-np.asarray(theta))[:n]
    return ', '.join(f'{model.species_names[i]} {theta[i]:.3f}' for i in order)


def converged(res):
    return res.get('status') == 'done' and bool(res['info_a']['converged'])


def results_frame(ctx, results, extra_cols=()):
    """One row per condition: row identifiers + evaluate() + convergence flags."""
    rows = []
    for key, res in results.items():
        r = res['row']
        out = {'key': key, 'basis': r['basis'], 'case': r['case'], 'T_K': r['T_K'],
               'P_ref_bar': r.get('P_ref_bar'), 'series': r.get('series'),
               'P_total_bar': r['P_total_bar'], 'P_probe_bar': r['P_probe_bar'],
               'P_impurity_bar': r['P_impurity_bar'], **{c: r.get(c) for c in extra_cols}}
        ok = converged(res)
        if ok:
            out.update(evaluate(ctx.model, r, res['theta_a']))
            info_b = res.get('info_b') or {}
            if info_b.get('converged') and res.get('theta_b') is not None:     # start B's steady state
                free, _, _, imp = role_coverages(ctx.model, res['theta_b'])
                out.update({'tof_fwd_B': info_b['tof_fwd_1_per_s'], 'theta_free_B': free, 'theta_impurity_B': imp})
        out.update({'converged': ok, 'starts_agree': res.get('starts_agree'),
                    'start_A_solved_by': res['info_a']['method'] if res.get('status') == 'done' else res.get('status'),
                    'start_B_solved_by': (res.get('info_b') or {}).get('method', res.get('status')),
                    'residual_1_per_s': res['info_a'].get('residual_1_per_s') if ok else np.nan})
        rows.append(out)
    return pd.DataFrame(rows)


def add_poisoning(df, ref_case, n_sites, by=('basis', 'T_K', 'P_ref_bar')):
    """poisoning_factor = TOF / TOF(reference) at the same `by` (the same point); with the site-blocking check
    TOF/TOF_ref = (theta*/theta*_ref)^n (p/p_ref) (f/f_ref)^n."""
    by = list(by)
    cols = ['tof_fwd', 'theta_free', 'P_probe_bar', 'site_count_factor']
    if 'tof_fwd' not in df:
        df = df.assign(**{c: np.nan for c in cols})
    ref = df[df['case'] == ref_case].set_index(by)[cols].add_suffix('_ref')
    df = df.join(ref, on=by)
    df['poisoning_factor'] = df['tof_fwd'] / df['tof_fwd_ref']
    df['slowdown_x'] = 1.0 / df['poisoning_factor']
    if 'tof_fwd_B' in df:
        df['poisoning_factor_B'] = df['tof_fwd_B'] / df['tof_fwd_ref']
    df['blocking_from_theta'] = (df['theta_free'] / df['theta_free_ref']) ** n_sites * df['P_probe_bar'] / df['P_probe_bar_ref']
    df['blocking_full'] = df['blocking_from_theta'] * (df['site_count_factor'] / df['site_count_factor_ref']) ** n_sites
    df['blocking_check'] = df['poisoning_factor'] / df['blocking_full']
    return df.drop(columns=[c + '_ref' for c in cols])


def problems(df):
    return df[~df['converged'].astype(bool) | (df['starts_agree'] == False)]   # noqa: E712
