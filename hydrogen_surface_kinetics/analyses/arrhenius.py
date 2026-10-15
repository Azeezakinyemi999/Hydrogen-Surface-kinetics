"""Arrhenius form of the hydrogen rates (documentation: Stage 11.2c).

Two levels, each fitted to y = A exp(-Ea / RT) over the temperature-sweep range:

intrinsic  the metal's own H2 step, independent of impurities (no solving): per site,
           k_diss = forward rate on an EMPTY surface in the probe gas at reference_pressure,
           k_rec  = reverse rate on a surface FULL of the probe adsorbate,
           K      = k_diss / k_rec at p = 1 bar: the equilibrium constant
                    theta_A^2 / (theta*^2 p/P0) of A2 + 2* <=> 2A*, fitted in van 't Hoff form
                    ln K = -dH/RT + dS/R.
observed   the steady-state rates of the temperature sweep, per case (reference, each
           impurity, mixture): forward TOF (dissociation), reverse TOF (recombination) and
           their ratio (1 = the H2 step is at equilibrium). Impurities act through the
           coverages, so their effect shows here: dEa and A / A_ref against the reference case.

A rate that is not a single exponential (coverages changing with T) gets a low R^2; the
local Ea(T) = -R dlnk/d(1/T) is then the number to read.
Outputs: arrhenius_intrinsic.csv (values per T), arrhenius_fits.csv (one row per fit, with its
conditions), arrhenius_local.csv (local Ea per case and T), arrhenius.png.
"""
import os

import cantera as ct
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from ..labels import probe_label
from . import common
from .diagnostics import reaction_thermo
from .sweep import sweep_temperatures

P0 = 1e5                      # Pa: standard pressure of the equilibrium constant
QUANTITIES = (('tof_fwd', 'dissociation (forward TOF)'), ('tof_rev', 'recombination (reverse TOF)'),
              ('ratio', 'forward / reverse'))


def fit(T, y):
    """ln y = ln A - Ea/(RT) by least squares. Returns dict(A, Ea_kJ_per_mol, R2, n, T_min, T_max);
    NaNs when fewer than two positive finite points."""
    T, y = np.asarray(T, dtype=float), np.asarray(y, dtype=float)
    ok = np.isfinite(T) & np.isfinite(y) & (y > 0)
    T, y = T[ok], y[ok]
    out = {'n_points': int(len(T)), 'T_min_K': T.min() if len(T) else np.nan, 'T_max_K': T.max() if len(T) else np.nan}
    if len(T) < 2:
        return {**out, 'A': np.nan, 'Ea_kJ_per_mol': np.nan, 'R2': np.nan}
    x, ly = 1.0 / T, np.log(y)
    slope, intercept = np.polyfit(x, ly, 1)
    resid = ly - (slope * x + intercept)
    ss_tot = ((ly - ly.mean()) ** 2).sum()
    # a flat quantity (e.g. forward/reverse = 1 at equilibrium) has no meaningful R^2
    r2 = 1.0 - (resid ** 2).sum() / ss_tot if np.std(ly) > 1e-3 else np.nan
    return {**out, 'A': float(np.exp(intercept)), 'Ea_kJ_per_mol': float(-slope * ct.gas_constant / 1e6), 'R2': float(r2)}


def local_ea(T, y):
    """-R dln y / d(1/T) by central differences (kJ/mol); NaN where y is not positive."""
    T, y = np.asarray(T, dtype=float), np.asarray(y, dtype=float)
    out = np.full(len(T), np.nan)
    ok = np.isfinite(y) & (y > 0)
    if ok.sum() >= 3:
        out[ok] = -ct.gas_constant * np.gradient(np.log(y[ok]), 1.0 / T[ok]) / 1e6
    return out


def intrinsic_table(ctx, temps):
    """Per-site rates of the probe step on two reference surfaces, and K, at each T."""
    model = ctx.model
    G, i, sign = model.gamma, model.i_probe, model.sign_probe
    rows = []
    for T in temps:
        rates = {}
        for name, cov in (('clean', {model.site: 1.0}), ('full', {model.probe_adsorbate: 1.0})):
            surf = model.new_surface()
            model.gas.TPX = T, ctx.cfg.reference_pressure, {model.probe_gas: 1.0}
            surf.TP = T, ctx.cfg.reference_pressure
            surf.coverages = cov
            fwd, rev = surf.forward_rates_of_progress[i], surf.reverse_rates_of_progress[i]
            rates[name] = (fwd, rev) if sign == 1 else (rev, fwd)        # in the adsorption direction
        k_diss = rates['clean'][0] / G
        k_rec = rates['full'][1] / G
        K = k_diss / k_rec * P0 / ctx.cfg.reference_pressure if k_rec > 0 else np.nan
        reac = {model.probe_gas: 1, model.site: model.probe_sites}
        prod = {model.probe_adsorbate: model.probe_sites}
        dH, dG = reaction_thermo(model, reac, prod, T)
        rows.append({'T_K': T, 'k_diss_clean_1_per_s': k_diss, 'k_rec_full_1_per_s': k_rec, 'K_eq': K,
                     'dH_thermo_kJ_per_mol': dH, 'dS_thermo_J_per_mol_K': (dH - dG) / T * 1e3})
    return pd.DataFrame(rows)


def mechanism_parameters(model):
    """The probe reaction as the mechanism gives it, e.g. 'sticking-coefficient: A=0.05, b=0, Ea=0 J/mol'."""
    r = model.surf.reaction(model.i_probe)
    data = r.input_data
    for key in ('sticking-coefficient', 'rate-constant'):
        if key in data:
            p = data[key]
            return f"{r.equation}; {key}: " + ', '.join(f'{k}={v}' for k, v in p.items())
    return r.equation


def case_conditions(ctx, case, row):
    """'p(H2) = 1 bar; CO 0.0053 bar (r = 0.00526)': the gas of one case at its point."""
    cfg, model = ctx.cfg, ctx.model
    probe = probe_label(model)
    ratio = common.setting(ctx, 'sweep', 'ratio', cfg.impurity_ratio)
    parts = [f"p({probe}) = {row['P_probe_bar']:.3g} bar"]
    for label in C.impurities_in(cfg, case):
        r = C.ratio_for(ratio, label)
        parts.append(f"{label} {r * row['P_probe_bar']:.3g} bar (r = {r:.3g})")
    return '; '.join(parts)


def sweep_table(ctx):
    """The temperature sweep from this run, or from temperature_sweep.csv of an earlier one."""
    if 'sweep' in ctx.outputs:
        return ctx.outputs['sweep']['sweep']
    path = os.path.join(ctx.cfg.output_dir, 'temperature_sweep.csv')
    if os.path.exists(path):
        ctx.log('  observed rates: from temperature_sweep.csv of an earlier run')
        return pd.read_csv(path)
    return None


def run(ctx):
    cfg, model = ctx.cfg, ctx.model
    ctx.log('\n== arrhenius: hydrogen rates in Arrhenius form ==')
    temps = [float(t) for t in (cfg.arrhenius['temperatures'] or sweep_temperatures(ctx))]
    probe = probe_label(model)
    intr = intrinsic_table(ctx, temps)
    ctx.save_table(intr, 'arrhenius_intrinsic.csv')
    clean = f"empty surface, p({probe}) = {cfg.reference_pressure / 1e5:g} bar"
    fits = []
    for col, what, cond in (('k_diss_clean_1_per_s', 'dissociation k_diss (intrinsic)', clean),
                            ('k_rec_full_1_per_s', 'recombination k_rec (intrinsic)',
                             f'surface full of {probe} adsorbate'),
                            ('K_eq', 'equilibrium K = k_diss / k_rec (intrinsic)', 'per bar; van t Hoff form')):
        f = fit(intr['T_K'], intr[col])
        row = {'level': 'intrinsic', 'quantity': what, 'case': '(clean surface)', 'conditions': cond, **f}
        if col == 'K_eq':                  # ln K = -dH/RT + dS/R
            row.update({'dH_kJ_per_mol': f['Ea_kJ_per_mol'], 'dS_J_per_mol_K': ct.gas_constant / 1e3 * np.log(f['A'])
                        if np.isfinite(f['A']) else np.nan})
        fits.append(row)
    ctx.log(f'  probe reaction in the mechanism: {mechanism_parameters(model)}')

    sweep = sweep_table(ctx)
    local_rows = []
    if sweep is None:
        ctx.log('  no temperature sweep: the observed rates are skipped (add "sweep" to analyses)')
    else:
        ok = sweep[sweep['converged'].astype(bool)].copy()
        ok['ratio'] = ok['tof_fwd'] / ok['tof_rev']
        for case in [c for c in common.sweep_cases(ctx, 'sweep') if c in set(ok['case'])]:
            part = ok[ok['case'] == case].sort_values('T_K')
            cond = case_conditions(ctx, case, part.iloc[0])
            for col, what in QUANTITIES:
                fits.append({'level': 'observed', 'quantity': what, 'case': case, 'conditions': cond,
                             **fit(part['T_K'], part[col])})
            local_rows += [{'case': case, 'T_K': T, 'Ea_local_fwd_kJ_per_mol': a, 'Ea_local_rev_kJ_per_mol': b,
                            'fwd_over_rev': rr}
                           for T, a, b, rr in zip(part['T_K'], local_ea(part['T_K'], part['tof_fwd']),
                                                  local_ea(part['T_K'], part['tof_rev']), part['ratio'])]
    fits = pd.DataFrame(fits)
    obs = fits['level'] == 'observed'
    if obs.any():                                                       # impurity effect vs the reference
        ref = fits[obs & (fits['case'] == cfg.reference_label)].set_index('quantity')
        fits.loc[obs, 'dEa_vs_reference_kJ_per_mol'] = [
            r.Ea_kJ_per_mol - ref['Ea_kJ_per_mol'].get(r.quantity, np.nan) for r in fits[obs].itertuples()]
        fits.loc[obs, 'A_over_reference'] = [r.A / ref['A'].get(r.quantity, np.nan) for r in fits[obs].itertuples()]
    fits['curved'] = fits['R2'].notna() & (fits['R2'] < cfg.arrhenius['min_r2'])
    ctx.save_table(fits, 'arrhenius_fits.csv')
    local = pd.DataFrame(local_rows)
    if len(local):
        ctx.save_table(local, 'arrhenius_local.csv')
    show = fits[['level', 'quantity', 'case', 'A', 'Ea_kJ_per_mol', 'R2']]
    ctx.log(show.to_string(index=False, float_format=lambda v: f'{v:.4g}'))
    fig = plot_arrhenius(ctx, intr, fits, sweep)
    P.save(fig, ctx.path('arrhenius.png'), cfg)
    plt.close(fig)
    ctx.outputs['arrhenius'] = {'intrinsic': intr, 'fits': fits, 'local': local}
    return ctx.outputs['arrhenius']


def plot_arrhenius(ctx, intr, fits, sweep):
    """(a) intrinsic k_diss and k_rec, (b) K (van 't Hoff), (c) observed forward TOF per case,
    (d) apparent Ea of the forward and reverse TOF per case; fit lines dashed."""
    cfg = ctx.cfg
    probe = probe_label(ctx.model)
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.4))
    x = 1000.0 / intr['T_K']

    def line(ax, row, color, T):
        T = np.linspace(row['T_min_K'], row['T_max_K'], 50) if np.isfinite(row['T_min_K']) else T
        if np.isfinite(row['A']):
            ax.plot(1000.0 / T, row['A'] * np.exp(-row['Ea_kJ_per_mol'] * 1e6 / (ct.gas_constant * T)),
                    ls='--', color=color, lw=1.0)

    intr_fit = fits[fits['level'] == 'intrinsic'].set_index('quantity')
    for col, q, color, lab in (('k_diss_clean_1_per_s', 'dissociation k_diss (intrinsic)', P.SERIES_COLORS[0],
                                f'dissociation, empty surface, {cfg.reference_pressure / 1e5:g} bar'),
                               ('k_rec_full_1_per_s', 'recombination k_rec (intrinsic)', P.SERIES_COLORS[1],
                                'recombination, surface full of H*' if probe == 'H2' else 'recombination, full surface')):
        r = intr_fit.loc[q]
        axes[0][0].plot(x, intr[col], 'o', color=color, ms=5, label=f"{lab}: Ea = {r['Ea_kJ_per_mol']:.1f} kJ/mol")
        line(axes[0][0], r, color, intr['T_K'])
    axes[0][0].set_yscale('log')
    axes[0][0].legend(frameon=False, fontsize=8, loc='best')
    P.style_axis(axes[0][0], 'a) Intrinsic rates of the H2 step (clean metal)', 'rate per site (1/s)')
    k = intr_fit.loc['equilibrium K = k_diss / k_rec (intrinsic)']
    axes[0][1].plot(x, intr['K_eq'], 'o', color=P.INK, ms=5)
    line(axes[0][1], k, P.INK, intr['T_K'])
    axes[0][1].set_yscale('log')
    slope_up = len(intr) > 1 and intr['K_eq'].iloc[-1] < intr['K_eq'].iloc[0]       # K falls with T: rises to the right
    axes[0][1].text(0.03, 0.94 if slope_up else 0.06,
                    f"ΔH = {k['dH_kJ_per_mol']:.1f} kJ/mol, ΔS = {k['dS_J_per_mol_K']:.0f} J/mol/K",
                    transform=axes[0][1].transAxes, fontsize=8, color=P.SECONDARY_INK, va='top' if slope_up else 'bottom')
    P.style_axis(axes[0][1], 'b) Equilibrium constant (van t Hoff)', 'K = θ²(H*) / (θ*² p/P°)')
    obs = fits[fits['level'] == 'observed']
    if sweep is not None and len(obs):
        cases = [c for c in common.sweep_cases(ctx, 'sweep') if c in set(obs['case'])]
        colors = P.colors_for(cases, first=cfg.reference_label)
        ok = sweep[sweep['converged'].astype(bool)]
        for case in cases:
            d = ok[ok['case'] == case].sort_values('T_K')
            axes[1][0].plot(1000.0 / d['T_K'], d['tof_fwd'], 'o', color=colors[case], ms=4, label=case)
            r = obs[(obs['case'] == case) & (obs['quantity'] == QUANTITIES[0][1])].iloc[0]
            line(axes[1][0], r, colors[case], d['T_K'])
        axes[1][0].set_yscale('log')
        axes[1][0].legend(frameon=False, fontsize=7, ncol=2, loc='best')
        xs = np.arange(len(cases))
        for q, dx, color, lab in ((QUANTITIES[0][1], -0.18, P.SERIES_COLORS[0], 'forward (dissociation)'),
                                  (QUANTITIES[1][1], 0.18, P.SERIES_COLORS[1], 'reverse (recombination)')):
            vals = [obs[(obs['case'] == c) & (obs['quantity'] == q)]['Ea_kJ_per_mol'].iloc[0] for c in cases]
            axes[1][1].bar(xs + dx, vals, width=0.36, color=color, label=lab)
        axes[1][1].axhline(intr_fit.loc['dissociation k_diss (intrinsic)']['Ea_kJ_per_mol'], color=P.MUTED, ls='--', lw=1,
                           label='intrinsic dissociation (clean surface)')
        e_all = [v for v in obs['Ea_kJ_per_mol'] if np.isfinite(v)] + [intr_fit.loc['dissociation k_diss (intrinsic)']['Ea_kJ_per_mol']]
        lo, hi = min(0.0, min(e_all)), max(0.0, max(e_all))
        axes[1][1].set_ylim(lo - 0.08 * (hi - lo + 1), hi + 0.45 * (hi - lo + 1))     # headroom for the legend
        axes[1][1].set_xticks(xs, cases, rotation=35, ha='right')
        axes[1][1].legend(frameon=False, fontsize=8, loc='best')
    P.style_axis(axes[1][0], 'c) Observed forward TOF (on the real surface)', 'TOF (1/s)')
    P.style_axis(axes[1][1], 'd) Apparent activation energy per case', 'Ea (kJ/mol)')
    for ax in (axes[0][0], axes[0][1], axes[1][0]):
        ax.set_xlabel('1000 / T  (1/K)')
    fig.tight_layout()
    return fig
