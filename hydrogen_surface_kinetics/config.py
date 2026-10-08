"""Load and check a run's input.py.

An input file is plain Python: top-level variables become settings. Anything left out
takes the default below. Dictionary-valued sections (``solver``, ``run``, and one per
analysis) are merged key by key, so an input file only lists what it changes.
"""
import copy
import os
import runpy
import types

# --- Defaults -------------------------------------------------------------------------
SOLVER_DEFAULTS = dict(
    abs_tol=1e-8,                   # 1/s   floor on the allowed |d(theta)/dt|
    rel_gross_tol=1e-12,            # -     plus this fraction of each species' gross rate
    sum_tol=1e-8,                   # -     |sum(theta) - 1|
    neg_tol=-1e-12,                 # -     coverages above this count as numerical zero
    gas_tol=1e-9,                   # -     largest allowed change in gas mole fractions
    direct_min_T=400.0,             # K     below this, skip the direct (Newton) solve
    reactor_area_per_volume=1e-10,  # 1/m   tiny surface in a huge gas volume: fixed gas
    reactor_rtol=1e-10,
    reactor_atol=1e-22,
    reactor_t_max=1e8,              # s     longest surface time integrated
    plateau_tof_rtol=1e-4,          # -     TOF change allowed between checkpoints
    plateau_theta_atol=1e-4,        # -     coverage change allowed between checkpoints
    probe_covered_fraction=0.99,    # -     start B (probe gas only): this fraction of probe adsorbate
    pre_exposure_time=1.0,          # s     start B (impurities): time in the impurity gas alone
    multi_ss_theta_atol=1e-3,       # -     starts agree if coverages differ less than this ...
    multi_ss_tof_rtol=1e-2,         # -     ... and the forward TOF by less than this fraction
)

RUN_DEFAULTS = dict(
    workers=None,                   # parallel processes; None = half the logical CPUs
    job_timeout=1800.0,             # s     time limit per job
    recompute=False,                # True = ignore cached results
)

# One section per analysis. None = take it from the main settings (see each module).
ANALYSIS_DEFAULTS = dict(
    composition=dict(
        min_listed=1e-3,            # list species above this coverage
        min_plotted=0.02,           # a species gets its own bar segment if it reaches this anywhere
        max_named=7,                # at most this many named segments; the rest is "other"
    ),
    diagnostics=dict(
        adsorption_checks=[],       # (label, {reactants}, {products}, (lit. dH lo, hi) kJ/mol, note)
        switch_off=None,            # species whose reactions are switched off in a what-if; None = molecular_adsorbate
        switch_off_cases=None,      # cases for the what-if; None = reference + every impurity
        top_reactions=[],           # (case, T, species): largest sources/sinks of that species
        top_n=5,
    ),
    curves=dict(
        ratios=[1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.25],
        temperatures=None,          # None = the main temperatures
        impurities=None,            # None = every impurity
    ),
    sweep=dict(
        temperatures=None,          # None = from the lowest to the highest main temperature in 50 K steps
        step=50.0,
        ratio=None,                 # None = impurity_ratio
        impurities=None,
    ),
    orders=dict(temperatures=None, ratio=None, delta=0.10, impurities=None),
    rate_control=dict(
        temperatures=None,
        cases=None,                 # None = reference + every impurity
        ratio=None,
        drc_delta=0.05,             # relative change of each screened reaction's rate constants
        drc_top_n=20,               # reactions screened, by |net rate| at the base steady state
        trc_delta=1.0e6,            # J/kmol (1 kJ/mol) change of a species' free energy
        trc_min_coverage=1e-3,      # species screened for TRC (plus the role species)
    ),
    paths=dict(
        conditions=None,            # [(case, T), ...]; None = every main case at the lowest temperature
        elements=None,              # None = the elements of the probe gas
        threshold=1e-8,             # hide fluxes below this fraction of the largest
        flow_type='OneWayFlow',     # 'OneWayFlow' (needed near equilibrium) or 'NetFlow'
        formats=['png'],            # Graphviz output formats, e.g. ['png', 'pdf']
    ),
    transients=dict(
        temperatures=None,
        impurities=None,
        ratios=[1e-6, 1e-4, 1e-2, 0.25],
        t_start=1e-9, t_stop=1e8, per_decade=4,   # s: recording times
    ),
    # Optional, needs RMG-Py (+ RDKit, Graphviz, ffmpeg): flux-diagram videos of the RMG
    # input file's own batch reactors. Paths are relative to this input file.
    flux=dict(
        rmg_path=None,              # RMG-Py source folder, if rmgpy is not importable already
        rmg_input=None,             # the RMG input.py with the surfaceReactor blocks (required)
        seed_dir=None,              # seed mechanism folder with reactions.py + dictionary.txt (required)
        species_dir=None,           # structure images (drawn if missing); None = <output>/flux_species
        reactor_indices=None,       # None = every reactor in rmg_input; or e.g. [0, 6, 24]
        n_points=100,               # time points per reactor (~45 for a quick check)
        rtol=1e-10, atol=1e-22,     # Cantera's defaults fail on stiff 300 K reactors
        skip_existing=True,         # skip reactors whose folder already has flux_diagram.avi
    ),
)
ANALYSIS_ORDER = ['steady_state', *ANALYSIS_DEFAULTS]

DEFAULTS = dict(
    # files
    mechanism=None,                 # Cantera YAML (required)
    gas_phase='gas',
    surface_phase='surface',
    output_dir='output',            # relative to the input file
    expected_site_density=None,     # e.g. (3.148e-9, 'mol/cm2'); None = no check
    # roles (species names in the mechanism)
    probe_gas=None,                 # required, e.g. 'H2(6)'
    empty_site=None,                # required, e.g. 'site(10)'
    probe_adsorbate=None,           # required: the dissociated adsorbate, e.g. 'HX(11)'
    molecular_adsorbate=None,       # optional: molecularly adsorbed probe gas
    reference_label='probe only',   # case name of the probe-gas-only reference
    probe_label=None,               # name of the probe gas in plots; None = its formula, e.g. 'H2'
    species_labels={},              # readable plot names for surface species: {mechanism name: label}
    # impurities: label -> species name, or elemental composition dict
    impurities={},
    mixture_label='all',            # case name of all impurities together; None = no mixture case
    # operating conditions
    temperatures=[300.0],
    pressure_mode='fixed',          # 'fixed' or 'input'
    reference_pressure=1e5,         # Pa, probe-gas partial pressure in 'fixed' mode
    input_pressures={},             # {T: Pa} in 'input' mode
    impurity_ratio=0.25,            # p(impurity) / p(probe gas) on the fixed-probe basis
    mixture_mode='same_each',       # 'same_each' or 'split'
    include_input_basis=False,      # also run the fixed-mole-fraction basis
    input_basis=dict(single=dict(probe=0.8, impurity=0.2), mixture=dict(probe=0.75, each_impurity=0.05)),
    # what to run
    analyses=['steady_state'],
    solver=SOLVER_DEFAULTS,
    run=RUN_DEFAULTS,
    plot_style='screen',            # 'screen' or 'publication'
    publication_width_cm=17.8,      # 17.8 = double column, 8.5 = single column
    **ANALYSIS_DEFAULTS,
)

DICT_SECTIONS = ('solver', 'run', 'input_basis', *ANALYSIS_DEFAULTS)
PATH_KEYS = ('mechanism', 'output_dir')
REQUIRED = ('mechanism', 'probe_gas', 'empty_site', 'probe_adsorbate')
SITE_DENSITY_UNITS = {'mol/cm2': 10.0, 'mol/m2': 1e-3, 'kmol/m2': 1.0}   # factor to kmol/m2


class ConfigError(ValueError):
    pass


class Config(dict):
    """Settings with attribute access: cfg.temperatures or cfg['temperatures']."""

    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(f"config has no setting '{name}'") from None

    def __setattr__(self, name, value):
        self[name] = value

    def section(self, name, **overrides):
        """A dict-valued setting with keyword overrides, e.g. cfg.section('run', workers=4)."""
        merged = dict(self.get(name) or {})
        merged.update({k: v for k, v in overrides.items() if v is not None})
        return merged


def _merge(defaults, values):
    out = copy.deepcopy(defaults)
    for key, value in values.items():
        if key in DICT_SECTIONS and isinstance(value, dict):
            out[key] = {**out.get(key, {}), **value}
        else:
            out[key] = value
    return out


def load_config(path, **overrides):
    """Read an input.py, apply defaults and keyword overrides, resolve paths, and check it."""
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.exists(path):
        raise ConfigError(f'input file not found: {path}')
    namespace = runpy.run_path(path)
    values = {k: v for k, v in namespace.items()
              if not k.startswith('_') and not isinstance(v, (types.ModuleType, types.FunctionType, type))}
    cfg = Config(_merge(_merge(DEFAULTS, values), overrides))      # overrides merge into sections too
    cfg['input_file'] = path
    cfg['run_dir'] = os.path.dirname(path)
    for key in PATH_KEYS:
        if cfg.get(key):
            p = os.path.expanduser(str(cfg[key]))
            cfg[key] = p if os.path.isabs(p) else os.path.normpath(os.path.join(cfg.run_dir, p))
    check_config(cfg)
    return cfg


def site_density_kmol_m2(spec):
    """(value, unit) -> kmol/m2."""
    value, unit = spec
    return float(value) * SITE_DENSITY_UNITS[unit]


def check_config(cfg):
    problems = []
    for key in REQUIRED:
        if not cfg.get(key):
            problems.append(f"'{key}' is required")
    if cfg.get('mechanism') and not os.path.exists(cfg.mechanism):
        problems.append(f"mechanism file not found: {cfg.mechanism}")
    if not isinstance(cfg.impurities, dict):
        problems.append("'impurities' must be a dict: label -> species name or composition dict")
    else:
        for label, spec in cfg.impurities.items():
            if not isinstance(spec, (str, dict)):
                problems.append(f"impurity '{label}': give a species name or a composition dict")
    if cfg.reference_label in cfg.impurities or (cfg.mixture_label and cfg.mixture_label in cfg.impurities):
        problems.append("'reference_label' and 'mixture_label' must differ from the impurity labels")
    try:
        temps = [float(t) for t in cfg.temperatures]
        if not temps:
            raise ValueError
        cfg['temperatures'] = temps
    except (TypeError, ValueError):
        problems.append("'temperatures' must be a non-empty list of numbers (K)")
        temps = []
    if cfg.pressure_mode not in ('fixed', 'input'):
        problems.append("'pressure_mode' must be 'fixed' or 'input'")
    elif cfg.pressure_mode == 'input':
        missing = [t for t in temps if t not in {float(k) for k in cfg.input_pressures}]
        if missing:
            problems.append(f"'input_pressures' has no entry for T = {missing}")
        cfg['input_pressures'] = {float(k): float(v) for k, v in cfg.input_pressures.items()}
    if cfg.mixture_mode not in ('same_each', 'split'):
        problems.append("'mixture_mode' must be 'same_each' or 'split'")
    if cfg.expected_site_density is not None:
        spec = cfg.expected_site_density
        if not (isinstance(spec, (tuple, list)) and len(spec) == 2 and spec[1] in SITE_DENSITY_UNITS):
            problems.append(f"'expected_site_density' must be (value, unit) with unit in {list(SITE_DENSITY_UNITS)}")
    if cfg.plot_style not in ('screen', 'publication'):
        problems.append("'plot_style' must be 'screen' or 'publication'")
    analyses = [cfg.analyses] if isinstance(cfg.analyses, str) else list(cfg.analyses)
    bad = [a for a in analyses if a not in ANALYSIS_ORDER]
    if bad:
        problems.append(f"unknown analyses {bad}; choose from {ANALYSIS_ORDER}")
    cfg['analyses'] = analyses
    for name, defaults in ANALYSIS_DEFAULTS.items():
        unknown = set(cfg[name]) - set(defaults)
        if unknown:
            problems.append(f"unknown '{name}' settings: {sorted(unknown)}")
    unknown = set(cfg.solver) - set(SOLVER_DEFAULTS)
    if unknown:
        problems.append(f"unknown solver settings: {sorted(unknown)}")
    unknown = set(cfg.run) - set(RUN_DEFAULTS)
    if unknown:
        problems.append(f"unknown run settings: {sorted(unknown)}")
    if problems:
        raise ConfigError(f"problems in {cfg.get('input_file', 'input')}:\n  - " + '\n  - '.join(problems))
