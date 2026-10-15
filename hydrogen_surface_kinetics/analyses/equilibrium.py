"""Equilibrium coverage of the probe step (documentation: Stage 11.2c). No solving.

At equilibrium, A2 + n* <=> n A*  gives  theta_A / theta_v = (K p / P0)^(1/n), with the
dimensionless K = theta_A^n / (theta_v^n p/P0) and P0 = 1 bar. K comes from the mechanism's
thermochemistry, exp(-dG0/RT) converted from its reference pressure to P0; the rates give
the same K in a thermodynamically consistent mechanism (check_K_rates = 1). On a surface
holding only the probe adsorbate (theta_A + theta_v = 1), x = (K p/P0)^(1/n):

  theta_H,eq = x / (1 + x),  theta_v,eq = 1 / (1 + x),  C_H,eq = theta_H,eq Gamma (mol H / m2),
  J0 = kd p theta_v,eq^n: the exchange flux at equilibrium (mol H2 / m2 / s),
  T_1/2: the temperature at which the clean surface is half covered at p (K p/P0 = 1).

The equilibrium check q = (theta_H/theta_v) / (K p/P0)^(1/n) of every converged steady state
tells whether the probe step is at equilibrium there (q = 1), whatever covers the rest of
the surface: impurities lower theta_H and theta_v together but cannot change their ratio.
Outputs: equilibrium_coverage.csv, half_coverage_temperature.csv, equilibrium_check.csv,
equilibrium_coverage.png and equilibrium_coverage_per_area.png.
"""
import os

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import per_area as PA
from .. import plotting as P
from ..conditions import case_labels
from ..labels import probe_label
from .diagnostics import reaction_thermo

P0 = PA.P0
T_RANGE = (200.0, 3000.0)           # K: where T_1/2 is looked for
SOURCES = (('main', 'steady_state', 'summary', 'steady_state_results.csv'),
           ('temperature sweep', 'sweep', 'sweep', 'temperature_sweep.csv'),
           ('pressure sweep', 'pressure_sweep', 'sweep', 'pressure_sweep.csv'))


def K_thermo(model, T):
    """Dimensionless K of the probe step at T from the thermochemistry, at P0 = 1 bar."""
    reac = {model.probe_gas: 1, model.site: model.probe_sites}
    prod = {model.probe_adsorbate: model.probe_sites}
    _, dG = reaction_thermo(model, reac, prod, T)
    P_thermo = model.gas.species(model.probe_gas).thermo.reference_pressure
    return float(np.exp(-dG * 1e6 / (ct.gas_constant * T)) * P0 / P_thermo)


def coverage(K, p, n, G):
    """(theta_H, theta_v, C_H) of a surface holding only the probe adsorbate, at equilibrium."""
    x = (K * p / P0) ** (1.0 / n)
    return x / (1 + x), 1 / (1 + x), x / (1 + x) * G


def half_coverage_temperature(model, p, lo=T_RANGE[0], hi=T_RANGE[1]):
    """T (K) at which K(T) p/P0 = 1, by a scan and bisection of ln K + ln(p/P0); NaN with a
    note when it is outside [lo, hi]."""
    f = lambda T: np.log(K_thermo(model, T)) + np.log(p / P0)                  # noqa: E731
    Ts = np.arange(lo, hi + 1, 10.0)
    vals = np.array([f(T) for T in Ts])
    cross = np.nonzero(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if not len(cross):
        side = 'covered' if vals[-1] > 0 else 'less than half covered'
        return np.nan, f'{side} at every T from {lo:.0f} to {hi:.0f} K'
    a, b = Ts[cross[0]], Ts[cross[0] + 1]
    fa = f(a)
    for _ in range(60):
        m = 0.5 * (a + b)
        fm = f(m)
        if np.sign(fm) == np.sign(fa):
            a, fa = m, fm
        else:
            b = m
    return 0.5 * (a + b), ''


def grid(ctx):
    """(temperatures K, pressures Pa) of the coverage table."""
    cfg, s = ctx.cfg, ctx.cfg.equilibrium
    temps = s['temperatures'] or sorted(set(cfg.temperatures) | set(cfg.pressure_temperatures or []))
    if s['pressures'] is not None:
        return [float(t) for t in temps], sorted(float(p) for p in s['pressures'])
    span = set(cfg.pressures) | {cfg.reference_pressure}
    lo, hi = min(span), max(span)
    if lo == hi:
        lo, hi = lo / 100.0, hi * 100.0
    n = int(round(s['per_decade'] * np.log10(hi / lo))) + 1
    return [float(t) for t in temps], list(np.logspace(np.log10(lo), np.log10(hi), max(n, 2)))


def results_tables(ctx):
    """[(source, results table)] of the steady states of this run: this run's outputs, else
    the CSVs of an earlier run."""
    out = []
    for source, analysis, key, csv in SOURCES:
        if analysis in ctx.outputs:
            out.append((source, ctx.outputs[analysis][key]))
        elif os.path.exists(os.path.join(ctx.cfg.output_dir, csv)):
            out.append((source, pd.read_csv(os.path.join(ctx.cfg.output_dir, csv))))
    return out


def check_table(ctx):
    """q = (theta_H/theta_v) / (K p/P0)^(1/n) for every converged steady state."""
    model, n = ctx.model, ctx.model.probe_sites
    Ks, rows = {}, []
    for source, df in results_tables(ctx):
        df = df[df['converged'].astype(bool)] if 'converged' in df else df
        for r in df.to_dict('records'):
            T, p = float(r['T_K']), float(r['P_probe_bar']) * 1e5
            if T not in Ks:
                Ks[T] = K_thermo(model, T)
            ratio = r['theta_probe'] / r['theta_free'] if r['theta_free'] > 0 else np.nan
            expected = (Ks[T] * p / P0) ** (1.0 / n)
            rows.append({'source': source, 'basis': r.get('basis'), 'series': r.get('series'), 'case': r['case'], 'T_K': T,
                         'P_ref_bar': r.get('P_ref_bar'), 'theta_H': r['theta_probe'], 'theta_v': r['theta_free'],
                         'theta_H_over_theta_v': ratio, 'equilibrium_ratio': expected, 'q': ratio / expected})
    return pd.DataFrame(rows)


def run(ctx):
    cfg, model = ctx.cfg, ctx.model
    ctx.log('\n== equilibrium: equilibrium coverage of the probe step ==')
    G, n = PA.site_density_mol_m2(model), model.probe_sites
    temps, pressures = grid(ctx)
    rows = []
    for T in temps:
        K = K_thermo(model, T)
        clean = PA.clean_surface(model, T, cfg.reference_pressure)                # kd, and K from the rates
        for p in pressures:
            th, tv, C = coverage(K, p, n, G)
            rows.append({'T_K': T, 'p_bar': p / 1e5, 'K_eq': K, 'check_K_rates': clean['K_eq'] / K,
                         'theta_H_eq': th, 'theta_v_eq': tv, 'C_H_eq_mol_m2': C,
                         'J0_clean_mol_m2_s': clean['kd_clean_mol_m2_s_Pa'] * p * tv ** n})
    cov = pd.DataFrame(rows)
    ctx.save_table(cov, 'equilibrium_coverage.csv')
    half_p = sorted(set(cfg.pressures) | {cfg.reference_pressure})
    half = pd.DataFrame([{'p_bar': p / 1e5, 'T_half_K': t, 'note': note}
                         for p in half_p for t, note in [half_coverage_temperature(model, p)]])
    ctx.save_table(half, 'half_coverage_temperature.csv')
    ctx.log('  temperature at which the clean surface is half covered:\n'
            + half.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    check = check_table(ctx)
    if len(check):
        ctx.save_table(check, 'equilibrium_check.csv')
        ctx.log(f"  equilibrium check q = (theta_H/theta_v) / (K p/P0)^(1/n): {check['q'].min():.4g} to "
                f"{check['q'].max():.4g} over {len(check)} steady states (1 = the probe step is at equilibrium)")
    for stem, area in (('equilibrium_coverage', False), ('equilibrium_coverage_per_area', True)):
        fig = plot_equilibrium(ctx, temps, pressures, half, check, per_area=area)
        P.save(fig, ctx.path(f'{stem}.png'), cfg)
        plt.close(fig)
    ctx.outputs['equilibrium'] = {'coverage': cov, 'half': half, 'check': check}
    return ctx.outputs['equilibrium']


def plot_equilibrium(ctx, temps, pressures, half, check, per_area=False):
    """(a) isotherms: theta_H,eq against p, a line per temperature, with the steady-state
    theta_H of the probe gas alone; (b) theta_H,eq against T, a line per pressure, with T_1/2
    (diamonds on 0.5) and the steady-state theta_H; (c) the equilibrium check q per case and
    point. per_area: (a) and (b) as C_H (mol H / m2)."""
    cfg, model = ctx.cfg, ctx.model
    probe = probe_label(model)
    G, n = PA.site_density_mol_m2(model), model.probe_sites
    scale = G if per_area else 1.0
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
    ref = check[(check['case'] == cfg.reference_label) & (check['basis'] != 'input_X')] if len(check) else check
    p_all = list(pressures) + (list(ref['P_ref_bar'].astype(float) * 1e5) if len(ref) else [])
    p_fine = np.logspace(np.log10(min(p_all)), np.log10(max(p_all)), 80)          # covers every point drawn
    for i, T in enumerate(temps):
        K = K_thermo(model, T)
        color = P.SERIES_COLORS[i % len(P.SERIES_COLORS)]
        axes[0].plot(p_fine / 1e5, coverage(K, p_fine, n, G)[0] * scale, color=color, lw=1.6, label=f'{T:.0f} K')
        if len(ref):
            d = ref[np.isclose(ref['T_K'], T)]
            axes[0].plot(d['P_ref_bar'], d['theta_H'] * scale, 'o', color=color, ms=6, mfc='white', mew=1.4)
    T_all = list(temps) + [t for t in half['T_half_K'] if np.isfinite(t)] + (list(ref['T_K']) if len(ref) else [])
    T_lo, T_hi = min(T_all) * 0.8, max(T_all) * 1.2
    T_fine = np.linspace(T_lo, T_hi, 80)
    K_fine = np.array([K_thermo(model, T) for T in T_fine])
    for i, rec in enumerate(half.to_dict('records')):
        color = P.SERIES_COLORS[(i + 3) % len(P.SERIES_COLORS)]
        p = rec['p_bar'] * 1e5
        axes[1].plot(T_fine, coverage(K_fine, p, n, G)[0] * scale, color=color, lw=1.6, label=f"{rec['p_bar']:g} bar")
        if np.isfinite(rec['T_half_K']):
            axes[1].plot(rec['T_half_K'], 0.5 * scale, **P.HALF_MARK, mec=color)
        if len(ref):
            d = ref[np.isclose(ref['P_ref_bar'].astype(float), rec['p_bar'])]
            axes[1].plot(d['T_K'], d['theta_H'] * scale, 'o', color=color, ms=6, mfc='white', mew=1.4)
    if per_area:
        label = PA.axis_label('$C_\\mathrm{H,eq}$', 'mol H m$^{-2}$', model)
    else:
        label = 'θ(H*) at equilibrium' if probe == 'H2' else 'probe adsorbate coverage at equilibrium'
    axes[0].set_xscale('log')
    axes[0].set_xlabel(f'p({probe}) (bar)')
    axes[1].set_xlabel('T (K)')
    for ax, title in ((axes[0], 'a) Isotherms (clean surface)'), (axes[1], 'b) Against temperature (clean surface)')):
        ax.set_ylim(0, 1.02 * scale)
        ax.legend(frameon=False, fontsize=8, loc='best')
        P.style_axis(ax, title, label)
    main = check[(check['source'] == 'main') & (check['basis'] == 'fixed_probe')] if len(check) else check
    if len(main):
        cases = [c for c in case_labels(cfg) if c in set(main['case'])]
        pts = list(dict.fromkeys(zip(main['T_K'], main['P_ref_bar'])))
        colors = P.colors_for(cases, first=cfg.reference_label)
        for j, case in enumerate(cases):
            d = main[main['case'] == case]
            xs = [pts.index((T, Pb)) + (j - (len(cases) - 1) / 2) * 0.12 for T, Pb in zip(d['T_K'], d['P_ref_bar'])]
            axes[2].plot(xs, d['q'], 'o', color=colors[case], ms=6, mec='white', mew=0.8, label=case, ls='')
        axes[2].axhline(1.0, color=P.MUTED, lw=1, ls='--')
        axes[2].set_xticks(range(len(pts)), [f'{T:.0f} K\n{Pb:g} bar' for T, Pb in pts], fontsize=8)
        q = main['q'].to_numpy(dtype=float)
        spread = max(0.02, np.nanmax(np.abs(q - 1)) * 1.3) if np.isfinite(q).any() else 0.02
        axes[2].set_ylim(1 - spread, 1 + spread)
        axes[2].legend(frameon=False, fontsize=8, loc='best')
    P.style_axis(axes[2], 'c) Equilibrium check of the main steady states',
                 '$q = (\\theta_\\mathrm{H}/\\theta_v)\\,/\\,(K\\,p/P^\\circ)^{1/n}$')
    fig.tight_layout()
    return fig
