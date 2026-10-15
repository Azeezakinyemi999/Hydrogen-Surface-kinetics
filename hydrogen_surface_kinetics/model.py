"""The mechanism and the roles the analysis is about (documentation: Stage 1).

A Model loads the gas and surface phases once, checks the site density, finds the
role species, the probe reaction (probe gas + n empty sites <=> probe adsorbate) and the
optional molecular adsorption, and resolves the impurities. It creates a fresh surface
object for every solve (see Model.new_surface).
"""
import hashlib
import os
import re
import warnings

import cantera as ct
import numpy as np

from .config import site_density_kmol_m2

UNMATCHED_DUPLICATE = 'No duplicate found for declared duplicate reaction'


class ModelError(ValueError):
    pass


def equation_species(equation):
    """The species of a Cantera equation, grouped the way Kinetics::checkDuplicates groups
    reactions: each side's species once, both sides together (a reaction matches its reverse)."""
    eq = re.sub(r'\s\(\+\s*[^\s)]+\)', ' ', equation)                 # falloff third body (+M)
    sides = re.split(r'\s(?:<=>|=>|=)\s', eq, maxsplit=1)
    names = [n for side in sides for n in {t.split()[-1] for t in side.split(' + ') if t.strip()}]
    return tuple(sorted(n for n in names if n != 'M'))


def remove_unmatched_duplicates(source, target):
    """Copy the Cantera YAML `source` to `target`, without the `duplicate: true` of reactions
    that have no other duplicate-marked reaction with the same species (RMG output sometimes
    has them, and Cantera then refuses the file). The flag does not change any rate. Returns
    the number of flags removed; the first line of `target` records the source's sha256."""
    import yaml
    with open(source, 'rb') as f:
        data = f.read()
    text = data.decode()
    root = yaml.compose(text, Loader=getattr(yaml, 'CSafeLoader', yaml.SafeLoader))
    flag_lines = {}                                                    # species -> lines of the flags
    for _, section in root.value:
        if not isinstance(section, yaml.SequenceNode):
            continue
        for node in section.value:
            if not isinstance(node, yaml.MappingNode):
                continue
            fields = {k.value: (k, v) for k, v in node.value if isinstance(k, yaml.ScalarNode)}
            dup = fields.get('duplicate')
            if 'equation' in fields and dup and str(dup[1].value).lower() in ('true', 'yes', 'on'):
                key = equation_species(fields['equation'][1].value)
                flag_lines.setdefault(key, []).append(dup[0].start_mark.line)
    lines = text.splitlines(keepends=True)
    drop = [found[0] for found in flag_lines.values() if len(found) == 1]
    for i in drop:
        if not re.fullmatch(r'\s*duplicate:\s*\w+\s*', lines[i]):
            raise ModelError(f'{source} line {i + 1}: expected "duplicate: true" on a line of its own')
        lines[i] = ''
    with open(target, 'w') as f:
        f.write(f'# {source} (sha256 {hashlib.sha256(data).hexdigest()}) without the duplicate flag '
                f'of {len(drop)} reactions that have no duplicate\n')
        f.writelines(lines)
    return len(drop)


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
        with open(cfg.mechanism, 'rb') as f:
            self.mechanism_sha256 = hashlib.sha256(f.read()).hexdigest()
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            self.mechanism_file = self._load_mechanism()
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

    # ----------------------------------------------------------------------------------
    def _load_mechanism(self):
        """Load the gas and surface phases and return the file they came from: the mechanism,
        or, if Cantera refuses it only for duplicate flags without a duplicate, a copy without
        those flags in the output folder (reused while the mechanism is unchanged)."""
        cfg = self.cfg
        stem = os.path.splitext(os.path.basename(cfg.mechanism))[0]
        fixed = os.path.join(cfg.output_dir, f'{stem}.duplicates_fixed.yaml')
        header = ''
        if os.path.exists(fixed):
            with open(fixed) as f:
                header = f.readline()
        if self.mechanism_sha256 not in header:
            try:
                self._load_phases(cfg.mechanism)
                return cfg.mechanism
            except ct.CanteraError as e:
                if UNMATCHED_DUPLICATE not in str(e):
                    raise
            os.makedirs(cfg.output_dir, exist_ok=True)
            remove_unmatched_duplicates(cfg.mechanism, fixed)
            with open(fixed) as f:
                header = f.readline()
        self._load_phases(fixed)
        n = re.search(r'of (\d+) reactions', header).group(1)
        self.warnings.append(f'{n} reactions in {os.path.basename(cfg.mechanism)} are marked duplicate but have '
                             f'no duplicate; loaded {fixed}, the same mechanism without those flags (rates unchanged)')
        return fixed

    def _load_phases(self, path):
        self.gas = ct.Solution(path, self.cfg.gas_phase, transport_model=None)
        self.surf = ct.Interface(path, self.cfg.surface_phase, [self.gas])

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
            surf = ct.Interface(self.mechanism_file, self.cfg.surface_phase, [self.gas])
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
