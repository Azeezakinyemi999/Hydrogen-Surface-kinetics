# H2 on Ni(111): first RMG mechanism (Danielle's run, first_try_saved).
# Anything not set here takes the package default (hydrogen_surface_kinetics/config.py).

# --- Files -----------------------------------------------------------------------------
# chem_annotated.yaml, NOT chem_annotated_Ni1.yaml (its site density is 10x too low)
mechanism = '../../../Surface-Impurity-Effect/rmg_runs/Ni/first_try_saved/chem_annotated.yaml'
output_dir = 'output'                            # relative to this file
expected_site_density = (3.148e-9, 'mol/cm2')    # surfaceSiteDensity in the RMG input.py

# --- Species roles ---------------------------------------------------------------------
probe_gas = 'H2(6)'
empty_site = 'site(10)'
probe_adsorbate = 'HX(11)'                       # H2 + 2 site <=> 2 HX
molecular_adsorbate = '[H][H].[Pt](42)'          # H2 + site <=> H2X
reference_label = 'pure H2'

impurities = {                                   # found by composition, so renumbering is harmless
    'H2O': {'H': 2, 'O': 1},
    'NH3': {'N': 1, 'H': 3},
    'CH4': {'C': 1, 'H': 4},
    'CO':  {'C': 1, 'O': 1},
    'CO2': {'C': 1, 'O': 2},
}
mixture_label = 'all'

# --- Conditions ------------------------------------------------------------------------
temperatures = [300.0, 700.0, 1000.0]            # K
pressure_mode = 'fixed'                          # same P_H2 at every T
reference_pressure = 1.0e5                       # Pa
input_pressures = {300.0: 5.0e5, 700.0: 3.0e5, 1000.0: 1.0e5}   # used if pressure_mode = 'input'
impurity_ratio = 0.25                            # p(impurity) / p(H2): the RMG input's 80:20
mixture_mode = 'same_each'
include_input_basis = False

# --- What to run -----------------------------------------------------------------------
analyses = ['steady_state']
plot_style = 'publication'
