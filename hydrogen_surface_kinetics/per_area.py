"""Format 2: the package's per-site rates (TOF, 1/s) per area of surface, and what that
makes possible (documentation: docs/rates_units_equilibrium.md). No solving.

  J_diss, J_rec, J_net = TOF x Gamma               mol H2 / m2 / s
  J_exchange           = min(J_diss, J_rec)        the part of the gross flux that cycles
                                                    between gas and surface (J0 at equilibrium)
  kd_eff               = J_diss / p                mol H2 / m2 / s / Pa: the literature's adsorption
                                                    coefficient, on the surface as it is
  gamma_eff            = J_diss / Z                sticking probability per collision, with the
                         Z = p / sqrt(2 pi M R T)  impingement flux (mol H2 / m2 / s)
  C_H                  = theta_H x Gamma           mol H / m2: hydrogen on the surface

On a clean surface (theta_v = 1) J_diss = kd p and gamma_eff is the mechanism's sticking
coefficient. Recombination is the surface step J_rec = kr (theta_H Gamma)^n, with n the
probe step's number of sites (2 for H2: kr in m2 / mol / s); nothing subsurface.
Fluxes count probe-gas molecules (H2), surface concentrations atoms (H). Gamma is the
mechanism's site density (an ideal facet), so 'per m2' is per atomic surface area.
"""
import cantera as ct
import numpy as np
import pandas as pd

P0 = 1e5                            # Pa: standard pressure of K
R = ct.gas_constant / 1e3           # J / (mol K)

# column -> unit, for tables and the report
COLUMNS = {
    'J_diss_mol_m2_s': 'mol H2 m-2 s-1', 'J_rec_mol_m2_s': 'mol H2 m-2 s-1', 'J_net_mol_m2_s': 'mol H2 m-2 s-1',
    'J_exchange_mol_m2_s': 'mol H2 m-2 s-1', 'kd_eff_mol_m2_s_Pa': 'mol H2 m-2 s-1 Pa-1', 'gamma_eff': '-',
    'C_H_mol_m2': 'mol H m-2',
}


def site_density_mol_m2(model):
    """Gamma in mol of sites per m2 (Cantera keeps kmol/m2)."""
    return model.gamma * 1e3


def gamma_label(model):
    """'$\\Gamma$ = 2.50$\\times$10$^{-5}$ mol m$^{-2}$' for axis labels: the site density a
    per-area figure was converted with."""
    G = site_density_mol_m2(model)
    e = int(np.floor(np.log10(G)))
    return rf'$\Gamma$ = {G / 10 ** e:.2f}$\times$10$^{{{e}}}$ mol m$^{{-2}}$'


def axis_label(quantity, unit, model):
    """Two-line axis label of a per-area panel: the quantity and unit, then the site density."""
    return f'{quantity} ({unit})\n{gamma_label(model)}'


def molar_mass_kg_mol(model):
    return model.gas.molecular_weights[model.gas.species_index(model.probe_gas)] / 1e3


def impingement_flux(T, p, M):
    """Z = p / sqrt(2 pi M R T): mol of gas hitting each m2 per second (p in Pa, M in kg/mol)."""
    return np.asarray(p, dtype=float) / np.sqrt(2 * np.pi * M * R * np.asarray(T, dtype=float))


def add_columns(df, model):
    """The Format-2 columns of a results table with tof_fwd, tof_rev, tof_net, theta_probe,
    P_probe_bar and T_K (NaN where the rates are missing)."""
    if df is None or 'tof_fwd' not in df:
        return df
    G, M = site_density_mol_m2(model), molar_mass_kg_mol(model)
    num = lambda col: pd.to_numeric(df[col], errors='coerce') if col in df else np.nan    # noqa: E731
    p = num('P_probe_bar') * 1e5
    J_diss, J_rec = num('tof_fwd') * G, num('tof_rev') * G
    return df.assign(J_diss_mol_m2_s=J_diss, J_rec_mol_m2_s=J_rec, J_net_mol_m2_s=num('tof_net') * G,
                     J_exchange_mol_m2_s=np.fmin(J_diss, J_rec).where(J_diss.notna() & J_rec.notna()),
                     kd_eff_mol_m2_s_Pa=J_diss / p, gamma_eff=J_diss / impingement_flux(num('T_K'), p, M),
                     C_H_mol_m2=num('theta_probe') * G)


def probe_rates_on(model, T, p, coverages):
    """Forward and reverse rate of progress of the probe step, in the adsorption direction
    (kmol/m2/s), on a fresh surface with these coverages under the pure probe gas at p."""
    surf = model.new_surface()
    model.gas.TPX = T, p, {model.probe_gas: 1.0}
    surf.TP = T, p
    surf.coverages = coverages
    i = model.i_probe
    fwd, rev = surf.forward_rates_of_progress[i], surf.reverse_rates_of_progress[i]
    return (fwd, rev) if model.sign_probe == 1 else (rev, fwd)


def clean_surface(model, T, p):
    """The probe step alone at T (K) and probe-gas pressure p (Pa), from the mechanism:
    kd (clean-surface adsorption coefficient), kr (surface recombination coefficient),
    gamma_clean (sticking probability on an empty surface), K (dimensionless equilibrium
    constant theta_H^n / (theta_v^n p/P0)), and the equilibrium of a surface holding only the
    probe adsorbate: theta_H_eq, theta_v_eq, C_H_eq and the exchange flux J0."""
    G, M, n = site_density_mol_m2(model), molar_mass_kg_mol(model), model.probe_sites
    fwd, _ = probe_rates_on(model, T, p, {model.site: 1.0})
    _, rev = probe_rates_on(model, T, p, {model.probe_adsorbate: 1.0})
    kd = fwd * 1e3 / p                              # J_diss = kd p on an empty surface
    kr = rev * 1e3 / G ** n                         # J_rec = kr (theta_H G)^n on a full surface
    K = kd / kr * P0 / G ** n if kr > 0 else np.nan
    x = (K * p / P0) ** (1.0 / n)                   # theta_H / theta_v at equilibrium
    theta_H = x / (1 + x)
    theta_v = 1 / (1 + x)
    return {'T_K': float(T), 'p_Pa': float(p), 'kd_clean_mol_m2_s_Pa': kd, 'kr_surface_m2_mol_s': kr,
            'gamma_clean': kd * np.sqrt(2 * np.pi * M * R * T), 'K_eq': K,
            'theta_H_eq': theta_H, 'theta_v_eq': theta_v, 'C_H_eq_mol_m2': theta_H * G,
            'J0_clean_mol_m2_s': kd * p * theta_v ** n}
