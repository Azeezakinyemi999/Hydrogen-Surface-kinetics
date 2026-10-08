"""The mechanism and the roles the analysis is about (documentation: Stage 1).

A Model loads the gas and surface phases once, checks the site density, finds the
role species, the probe reaction (probe gas + n empty sites <=> probe adsorbate) and the
optional molecular adsorption, and resolves the impurities. It creates a fresh surface
object for every solve (see Model.new_surface).
"""
import hashlib
import warnings

import cantera as ct
import numpy as np

from .config import site_density_kmol_m2


class ModelError(ValueError):
    pass


def find_role_reaction(surf, gas_species, site, adsorbate):
    """The single reaction  gas_species + n site <=> m adsorbate  (either direction).
    Returns (index, sign, n_sites): sign = +1 if written in that direction, else -1."""
    left, right = {gas_species, site}, {adsorbate}
    hits = []
    for i, r in enumerate(surf.reactions()):
        if set(r.reactants) == left and set(r.products) == right:
            hits.append((i, 1, r.reactants[site]))
        elif set(r.products) == left and set(r.reactants) == right:
            hits.append((i, -1, r.products[site]))
    if len(hits) != 1:
        raise ModelError(f'expected exactly one reaction {gas_species} + site <=> {adsorbate}, '
                         f'found {len(hits)}: {[surf.reaction(h[0]).equation for h in hits]}')
    return hits[0]


class Model:
    def __init__(self, cfg):
        self.cfg = cfg
        self.warnings = []
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.gas = ct.Solution(cfg.mechanism, cfg.gas_phase, transport_model=None)
            self.surf = ct.Interface(cfg.mechanism, cfg.surface_phase, [self.gas])
        self.sticking_warnings = [str(w.message).strip() for w in caught if 'Sticking coefficient' in str(w.message)]

        # Site density
        self.gamma = self.surf.site_density                                  # kmol/m^2
        if cfg.expected_site_density is not None:
            expected = site_density_kmol_m2(cfg.expected_site_density)
            if not np.isclose(self.gamma, expected, rtol=1e-3):
                raise ModelError(f'site density {self.gamma:.4e} kmol/m2 does not match the expected '
                                 f'{expected:.4e} kmol/m2 ({cfg.expected_site_density}): wrong or mis-converted file?')

        # Role species
        for role in ('empty_site', 'probe_adsorbate', 'molecular_adsorbate'):
            name = cfg.get(role)
            if name and name not in self.surf.species_names:
                raise ModelError(f"{role} '{name}' is not a surface species in {cfg.mechanism}")
        if cfg.probe_gas not in self.gas.species_names:
            raise ModelError(f"probe_gas '{cfg.probe_gas}' is not a gas species in {cfg.mechanism}")
        self.probe_gas, self.site = cfg.probe_gas, cfg.empty_site
        self.probe_adsorbate, self.molecular_adsorbate = cfg.probe_adsorbate, cfg.molecular_adsorbate

        # Role reactions
        self.i_probe, self.sign_probe, self.probe_sites = find_role_reaction(
            self.surf, self.probe_gas, self.site, self.probe_adsorbate)
        self.i_molec = self.sign_molec = None
        if self.molecular_adsorbate:
            self.i_molec, self.sign_molec, _ = find_role_reaction(
                self.surf, self.probe_gas, self.site, self.molecular_adsorbate)
        others = [i for i, r in enumerate(self.surf.reactions())
                  if (self.probe_gas in r.reactants or self.probe_gas in r.products)
                  and i not in (self.i_probe, self.i_molec)]
        if others:
            self.warnings.append(f'{len(others)} other surface reactions use the probe gas directly; '
                                 f'net probe-gas consumption includes them (first: {self.surf.reaction(others[0]).equation})')

        # Impurities
        self.impurities = {label: self._resolve_gas(label, spec) for label, spec in cfg.impurities.items()}

        # Arrays used by every solve (same mechanism, so the same for every fresh object)
        self.site_size = np.array([s.size for s in self.surf.species()])
        self.surf_kin_idx = np.array([self.surf.kinetics_species_index(n) for n in self.surf.species_names])
        self.probe_kin_idx = self.surf.kinetics_species_index(self.probe_gas)
        with open(cfg.mechanism, 'rb') as f:
            self.mechanism_sha256 = hashlib.sha256(f.read()).hexdigest()

    # ----------------------------------------------------------------------------------
    def _resolve_gas(self, label, spec):
        if isinstance(spec, str):
            if spec not in self.gas.species_names:
                raise ModelError(f"impurity '{label}': '{spec}' is not a gas species")
            return spec
        target = {el: float(n) for el, n in spec.items()}
        hits = [s.name for s in self.gas.species() if s.composition == target]
        if len(hits) != 1:
            raise ModelError(f"impurity '{label}': expected one gas species with composition {spec}, found {hits}; "
                             'give the species name instead')
        return hits[0]

    @property
    def species_names(self):
        return self.surf.species_names

    def new_surface(self, multipliers=None, h_shifts=None):
        """A fresh Interface object. Cantera's steady-state solver keeps state attached to
        the object, so every solve gets its own. Optional what-if changes:
          multipliers: {reaction index: factor}  (scales forward and reverse rates)
          h_shifts:    {surface species: J/kmol} (enthalpy shift; NASA7 thermo only)"""
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            surf = ct.Interface(self.cfg.mechanism, self.cfg.surface_phase, [self.gas])
        for name, dh in (h_shifts or {}).items():
            k = surf.species_index(name)
            sp = surf.species(k)
            c = np.array(sp.thermo.coeffs)
            if isinstance(sp.thermo, ct.NasaPoly2):
                # [T_mid, a1..a7 high-T, a1..a7 low-T]: a6 (the enthalpy constant) of both ranges
                c[[6, 13]] += dh / ct.gas_constant
                sp.thermo = ct.NasaPoly2(sp.thermo.min_temp, sp.thermo.max_temp, sp.thermo.reference_pressure, c)
            elif isinstance(sp.thermo, ct.ConstantCp):
                c[1] += dh                      # [T0, h0, s0, cp0]
                sp.thermo = ct.ConstantCp(sp.thermo.min_temp, sp.thermo.max_temp, sp.thermo.reference_pressure, c)
            else:
                raise ModelError(f'{name}: only NASA7 or constant-cp thermo can be shifted')
            surf.modify_species(k, sp)
        for i, m in (multipliers or {}).items():
            surf.set_multiplier(m, i)
        self.surf = surf
        return surf

    def probe_rates(self, surf=None):
        """Forward, reverse, net rate of probe gas + n site -> n' probe adsorbate (kmol/m2/s)."""
        surf = surf or self.surf
        fwd, rev = surf.forward_rates_of_progress[self.i_probe], surf.reverse_rates_of_progress[self.i_probe]
        if self.sign_probe == -1:
            fwd, rev = rev, fwd
        return fwd, rev, fwd - rev

    def fingerprint_context(self):
        """Everything besides a job's own spec that affects a cached result."""
        return (self.mechanism_sha256, self.cfg.gas_phase, self.cfg.surface_phase,
                sorted(self.cfg.solver.items()), ct.__version__)

    def summary(self):
        lines = [f'mechanism: {self.cfg.mechanism}',
                 f'  gas: {self.gas.n_species} species, {self.gas.n_reactions} reactions; '
                 f'surface: {self.surf.n_species} species, {self.surf.n_reactions} reactions',
                 f'  site density: {self.gamma:.4e} kmol/m2',
                 f'  probe reaction #{self.i_probe}: {self.surf.reaction(self.i_probe).equation}  '
                 f'(direction {self.sign_probe:+d}, {self.probe_sites} empty sites)']
        if self.i_molec is not None:
            lines.append(f'  molecular adsorption #{self.i_molec}: {self.surf.reaction(self.i_molec).equation}')
        lines.append(f'  impurities: {self.impurities}')
        if self.sticking_warnings:
            lines.append(f'  {len(self.sticking_warnings)} sticking coefficients exceed 1 at some temperature')
        lines += [f'  warning: {w}' for w in self.warnings]
        return '\n'.join(lines)
