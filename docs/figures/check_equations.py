"""Check the equations in docs/simulation_workflow.md against Cantera itself.

Uses the toy mechanism (toy_mechanism.yaml, fictitious) so it runs in seconds and
needs nothing from a real RMG run. Run it after upgrading Cantera: every line
should print PASS. If one fails, the corresponding section of the documentation
no longer describes how your Cantera version behaves.

    python check_equations.py
"""
import os

import cantera as ct
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
YAML = os.path.join(HERE, 'toy_mechanism.yaml')
R = ct.gas_constant  # J/kmol/K

gas = ct.Solution(YAML, 'gas')
surf = ct.Interface(YAML, 'surface', [gas])
GAMMA = surf.site_density                                   # kmol/m^2
SIZE = np.array([s.size for s in surf.species()])
SURF_IDX = np.array([surf.kinetics_species_index(n) for n in surf.species_names])

T, P = 500.0, 2.0e5
gas.TPX = T, P, 'A2:0.6, I:0.3, AI:0.1'
surf.TP = T, P
surf.coverages = {'S': 0.4, 'A(s)': 0.3, 'I(s)': 0.2, 'I(ss)': 0.1}

results = []


def check(label, ok, detail=''):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {label}" + (f'   ({detail})' if detail else ''))


# 1. Sticking coefficient -> forward rate constant (Section 2.2)
for i, gas_sp, m in [(0, 'A2', 2), (1, 'I', 1)]:
    rate = surf.reaction(i).rate
    gamma = rate.pre_exponential_factor * T ** rate.temperature_exponent * np.exp(-rate.activation_energy / (R * T))
    W = gas.molecular_weights[gas.species_index(gas_sp)]
    kf_formula = gamma / GAMMA ** m * np.sqrt(R * T / (2 * np.pi * W))
    kf_cantera = surf.forward_rate_constants[i]
    check(f'sticking k_f, reaction {i} (m = {m})', np.isclose(kf_formula, kf_cantera, rtol=1e-10),
          f'formula {kf_formula:.6e}, Cantera {kf_cantera:.6e}')

# 2. Equilibrium constant in concentration units and reverse rate constant (Section 2.3)
#    Uses REFERENCE-state Gibbs energies (gas at its reference pressure P_ref) with
#    standard concentrations P_ref/RT (gas) and Gamma/sigma_k (surface).
def g_ref_RT(sp):
    return (sp.thermo.h(T) - T * sp.thermo.s(T)) / (R * T)


mu0_RT = {**{n: g_ref_RT(gas.species(n)) for n in gas.species_names},
          **{n: g_ref_RT(surf.species(n)) for n in surf.species_names}}
P_ref = gas.species('A2').thermo.reference_pressure
C0 = {**{n: P_ref / (R * T) for n in gas.species_names},
      **{n: GAMMA / s.size for n, s in zip(surf.species_names, surf.species())}}
for i, r in enumerate(surf.reactions()):
    nu = {k: -v for k, v in r.reactants.items()}
    for k, v in r.products.items():
        nu[k] = nu.get(k, 0) + v
    dG_RT = sum(v * mu0_RT[k] for k, v in nu.items())
    Kc_formula = np.exp(-dG_RT) * np.prod([C0[k] ** v for k, v in nu.items()])
    Kc_cantera = surf.equilibrium_constants[i]
    check(f'K_c, reaction {i} ({r.equation})', np.isclose(Kc_formula, Kc_cantera, rtol=1e-8),
          f'formula {Kc_formula:.6e}, Cantera {Kc_cantera:.6e}')
    kr_formula = surf.forward_rate_constants[i] / Kc_cantera
    check(f'k_r = k_f / K_c, reaction {i}', np.isclose(kr_formula, surf.reverse_rate_constants[i], rtol=1e-10))

# Cantera's gas "standard" Gibbs energies are at the CURRENT pressure, not P_ref
check('gas.standard_gibbs_RT = reference value + ln(P/P_ref)',
      np.isclose(gas.standard_gibbs_RT[0], mu0_RT['A2'] + np.log(P / P_ref), rtol=1e-12))

# 3. Mass-action rates of progress (Section 2.2)
conc = {**dict(zip(gas.species_names, gas.concentrations)), **dict(zip(surf.species_names, surf.concentrations))}
for i, r in enumerate(surf.reactions()):
    qf = surf.forward_rate_constants[i] * np.prod([conc[k] ** v for k, v in r.reactants.items()])
    qr = surf.reverse_rate_constants[i] * np.prod([conc[k] ** v for k, v in r.products.items()])
    check(f'mass action q_f, q_r, reaction {i}',
          np.isclose(qf, surf.forward_rates_of_progress[i], rtol=1e-10)
          and np.isclose(qr, surf.reverse_rates_of_progress[i], rtol=1e-10))

# 4. Coverages, concentrations and two-site species (Section 2.1)
theta = surf.coverages
C = surf.concentrations
S = np.sum(theta / SIZE)
check('sum of species concentrations equals the site density', np.isclose(C.sum(), GAMMA, rtol=1e-12),
      f'sum C_k / Gamma = {C.sum() / GAMMA:.6f}')
check('C_k = Gamma * theta_k / (sigma_k * S), with S = sum(theta_j / sigma_j)',
      np.allclose(C, GAMMA * theta / (SIZE * S), rtol=1e-12))
check('theta_k = sigma_k C_k / sum_j(sigma_j C_j)', np.allclose(theta, SIZE * C / np.sum(SIZE * C), rtol=1e-12))
f = np.sum(SIZE * C) / GAMMA
check('site-count factor f = sum(sigma_k C_k)/Gamma = 1/S', np.isclose(f, 1 / S, rtol=1e-12), f'f = {f:.4f}')

# 5. Rate multipliers scale forward AND reverse rates, so K is unchanged (Stages 10 and 11)
qf0, qr0 = surf.forward_rates_of_progress.copy(), surf.reverse_rates_of_progress.copy()
surf.set_multiplier(3.0, 2)
qf1, qr1 = surf.forward_rates_of_progress.copy(), surf.reverse_rates_of_progress.copy()
surf.set_multiplier(1.0, 2)
check('set_multiplier(3, i) triples forward and reverse rate of reaction i only',
      np.isclose(qf1[2], 3 * qf0[2]) and np.isclose(qr1[2], 3 * qr0[2])
      and np.allclose(np.delete(qf1, 2), np.delete(qf0, 2)))

# 6. Coverage equation used by Cantera's reactor: d(theta_k)/dt = sigma_k * sdot_k / Gamma (Section 2.4)
#    Compared with the integrator's own derivative during the fast transient. That
#    derivative is interpolated by CVODES, so agreement is only to ~0.1-1 %.
reactor = ct.IdealGasReactor(gas, energy='off')
reactor.volume = 1.0
ct.ReactorSurface(surf, reactor, A=1e-10)        # tiny surface: the gas cannot change
net = ct.ReactorNet([reactor])
net.rtol, net.atol = 1e-12, 1e-24
net.advance(1e-6)
surf_vars = [i for i in range(net.n_vars) if net.component_name(i).split(': ')[-1] in surf.species_names]
dtheta_cantera = net.get_derivative(1)[surf_vars]
dtheta_formula = surf.net_production_rates[SURF_IDX] * SIZE / GAMMA
dev = np.max(np.abs(dtheta_cantera - dtheta_formula)) / np.max(np.abs(dtheta_formula))
check('reactor integrates d(theta)/dt = sigma * sdot / Gamma', dev < 1e-2, f'max relative deviation {dev:.1e}')

print(f'\n{sum(results)} / {len(results)} checks passed with Cantera {ct.__version__}')
