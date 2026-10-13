# H2 on Ni(111): first RMG mechanism (Danielle's run, first_try_saved).
# Anything not set here takes the package default (hydrogen_surface_kinetics/config.py).

# --- Files -----------------------------------------------------------------------------
# chem_annotated.yaml, NOT chem_annotated_Ni1.yaml (its site density is 10x too low)
mechanism = '/Users/akinyemi.az/Desktop/PhD_Folder/research/MHI/Surface-Impurity-Effect/rmg_runs/Ni/first_try_saved/chem_annotated.yaml'
output_dir = 'output'                            # relative to this file
expected_site_density = (3.148e-9, 'mol/cm2')    # surfaceSiteDensity in the RMG input.py

# --- Species roles ---------------------------------------------------------------------
probe_gas = 'H2(6)'
empty_site = 'site(10)'
probe_adsorbate = 'HX(11)'                       # H2 + 2 site <=> 2 HX
molecular_adsorbate = '[H][H].[Pt](42)'          # H2 + site <=> H2X
reference_label = 'H2 only'
probe_label = 'H2'

impurities = {                                   # found by composition, so renumbering is harmless
    'H2O': {'H': 2, 'O': 1},
    'NH3': {'N': 1, 'H': 3},
    'CH4': {'C': 1, 'H': 4},
    'CO':  {'C': 1, 'O': 1},
    'CO2': {'C': 1, 'O': 2},
}
mixture_label = 'all'

# --- Conditions ------------------------------------------------------------------------
temperatures = [300.0, 700.0, 1000.0]            # K: temperature series, at reference_pressure
reference_pressure = 1.0e5                       # Pa (H2 partial pressure)
pressures = []                                   # no pressure series for this first mechanism
impurity_ratio = 0.25                            # p(impurity) / p(H2): the RMG input's 80:20
mixture_mode = 'same_each'
include_input_basis = False

species_labels = {
    'site(10)': 'empty site',
    'HX(11)': 'H*',
    '[H][H].[Pt](42)': 'H2* (molecular)',
    'OCX(13)': 'CO* (1 site)',
    'COXX(5137)': 'CO** (2 sites)',
    'OX(12)': 'O*',
    'HOX(15)': 'OH*',
    'CX(14)': 'C*',
    '[Pt].N(39)': 'NH3*',
    'N#[Pt](49)': 'N*',
    'C.[Pt](32)': 'CH4* (weakly bound)',
    '[Pt]#[C]OO[C]#[Pt](2071)': 'C2O2** (unusual)',
}

# --- What to run -----------------------------------------------------------------------
analyses = ['steady_state', 'composition', 'diagnostics', 'curves', 'sweep', 'orders',
            'rate_control', 'paths', 'transients']
plot_style = 'publication'

# --- Analysis settings (the mechanism is trusted at 300 K only, so most run there) ----
_H2, _S = 'H2(6)', 'site(10)'
diagnostics = {
    # Literature values are APPROXIMATE ranges for Ni(111) at low coverage: check them
    # against your own references before using them in a paper.
    'adsorption_checks': [
        ('H2 -> 2 H*',            {_H2: 1, _S: 2},          {'HX(11)': 2},           (-100, -90),  'Christmann et al. 1974/1979 (~ -96)'),
        ('H2 -> H2* (molecular)', {_H2: 1, _S: 1},          {'[H][H].[Pt](42)': 1},  (-10, 0),     'physisorption only; not stable on Ni(111) above ~50 K'),
        ('CO -> CO* (1 site)',    {'CO(7)': 1, _S: 1},      {'OCX(13)': 1},          (-135, -110), 'Stuckless et al. 1993 calorimetry (~ -130)'),
        ('CO -> CO** (2 sites)',  {'CO(7)': 1, _S: 2},      {'COXX(5137)': 1},       (-135, -110), 'same CO on Ni(111); bridge/hollow vs top differ by only a few kJ/mol'),
        ('1/2 O2 -> O*',          {'O2(2)': 0.5, _S: 1},    {'OX(12)': 1},           (-260, -220), 'approx., from O2 adsorption heat ~ -480 kJ/mol O2'),
        ('NH3 -> NH3*',           {'NH3(4)': 1, _S: 1},     {'[Pt].N(39)': 1},       (-80, -50),   'approx., TPD ~ 0.6-0.8 eV'),
        ('H2O -> H2O*',           {'H2O(5)': 1, _S: 1},     {'O.[Pt](40)': 1},       (-60, -40),   'approx., weak molecular adsorption'),
        ('CH4 -> CH4*',           {'CH4(1)': 1, _S: 1},     {'C.[Pt](32)': 1},       (-25, -10),   'physisorption'),
        ('CO2 -> CO2*',           {'CO2(3)': 1, _S: 1},     {'CO2X(16)': 1},         (-40, -10),   'physisorption / very weak'),
    ],
    'switch_off_cases': ['H2 only', 'CO'],           # molecular H2* switched off
    'top_reactions': [('CO2', 1000.0, '[H][H].[Pt](42)'), ('CO', 700.0, '[Pt]#[C]OO[C]#[Pt](2071)')],
}
curves = {'temperatures': [300.0]}
sweep = {'temperatures': [300.0 + 50 * i for i in range(15)]}     # 300-1000 K
orders = {'temperatures': [300.0]}
rate_control = {'temperatures': [300.0], 'cases': ['H2 only', 'CO', 'NH3']}
paths = {'conditions': [('H2 only', 300.0), ('CO', 300.0), ('NH3', 300.0)], 'elements': ['H']}
transients = {'temperatures': [300.0]}

# Optional (needs RMG-Py): flux-diagram videos of the RMG input's 27 batch reactors.
# Run with:  bash run.sh runs/ni_first_try --only flux
_RMG_RUN = '/Users/akinyemi.az/Desktop/PhD_Folder/research/MHI/Surface-Impurity-Effect/rmg_runs/Ni/first_try_saved'
flux = {
    'rmg_path': '/Users/akinyemi.az/CoMoCheng/WESTLAB/rms_rmg_new/RMG-Py',
    'rmg_input': _RMG_RUN + '/input.py',
    'seed_dir': _RMG_RUN + '/seed/seed',
}
