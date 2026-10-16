# H2 on Ni(111): RMG run try_2/Ni (Danielle), mechanism and RMG files in this folder.
# The RMG input.py next to this file is RMG's; this settings.py is the package's.
# RMG built this mechanism for 650-1200 K and 0.01-100 bar total pressure (see input.py),
# so both series stay inside that range.
# Before running: `bash run.sh runs/try_2/Ni check` shows the probe reaction found and how many
# solves each analysis needs (nothing is solved by check).

# --- Files -----------------------------------------------------------------------------
mechanism = 'chem_annotated.yaml'
output_dir = 'output'
expected_site_density = (3.148e-09, 'mol/cm2')    # RMG metal database, Ni111

# --- Species roles (names read from chem_annotated.yaml) --------------------------------
probe_gas = 'H2(1)'
empty_site = 'vacantX(8)'
probe_adsorbate = 'HX(14)'                # H2 + 2 site <=> 2 H*
molecular_adsorbate = '[H][H].[Pt](40)'   # H2 + site <=> H2*
reference_label = 'H2 only'
probe_label = 'H2'

impurities = {                                   # the gases of the RMG reactors, by name
    'H2O':  'H2O(2)',
    'NH3':  'NH3(9)',
    'CH4':  'CH4(5)',
    'CO':   'CO(4)',
    'CO2':  'CO2(3)',
    'O2':   'O2(6)',
    'N2':   'N2(7)',
}
mixture_label = 'all'

species_labels = {
    'vacantX(8)': 'empty site',
    'HX(14)': 'H*',
    '[H][H].[Pt](40)': 'H2* (molecular)',
    'OX(15)': 'O*',
    '[Pt].N(66)': 'NH3*',
    'O.[Pt](44)': 'H2O*',
    'C.[Pt](57)': 'CH4* (weakly bound)',
    'CO2X(19)': 'CO2*',
    'COXX(7400)': 'CO** (2 sites)',
    'OCX(16)': 'CO* (OCX)',
    'COX(7423)': 'CO* (COX)',
}

# --- Conditions: two separate series (temperature and pressure effects never mix) --------
temperatures = [650.0, 800.0, 1000.0, 1200.0]   # K: temperature series, at reference_pressure
reference_pressure = 1.0e5                      # Pa: H2 partial pressure of the temperature series
pressures = [1.0e3, 1.0e5, 1.0e7]               # Pa: pressure series (0.01, 1, 100 bar, as in the RMG reactors)
pressure_temperatures = [650.0, 1200.0]         # K: a pressure series at each
# p(impurity) / p(H2): the proportions of the first RMG reactor (H2 = 190), for every case.
# Each impurity alone sits at its own level, and 'all' is that reactor's gas. The same in both series.
impurity_ratio = {'H2O': 10 / 190, 'CO2': 0.05 / 190, 'CO': 1 / 190, 'CH4': 0.8 / 190,
                  'N2': 0.85 / 190, 'O2': 0.15 / 190, 'NH3': 0.1 / 190}
mixture_mode = 'same_each'
include_input_basis = False

# --- What to run -----------------------------------------------------------------------
analyses = ['steady_state', 'composition', 'diagnostics', 'convergence', 'curves', 'sweep', 'arrhenius',
            'pressure_sweep', 'equilibrium', 'orders', 'rate_control', 'paths', 'transients']
plot_style = 'publication'

# --- Analysis settings -----------------------------------------------------------------
# Every point-based analysis runs on both series by default. rate_control is by far the
# most expensive (about 50 solves per case and point): `check` shows the counts. To narrow
# an analysis, give it its own series / temperatures / pressures / pressure_temperatures /
# cases, e.g.  rate_control = {'cases': ['H2 only', 'CO'], 'series': 'temperature'}
_H2, _S = 'H2(1)', 'vacantX(8)'
diagnostics = {
    # Literature values are APPROXIMATE ranges for Ni(111) at low coverage (the same as in
    # runs/ni_first_try): check them against your own references before using them in a paper.
    'adsorption_checks': [
        ('H2 -> 2 H*', {_H2: 1, _S: 2}, {'HX(14)': 2}, (-100, -90), 'Christmann et al. 1974/1979 (~ -96)'),
        ('H2 -> H2* (molecular)', {_H2: 1, _S: 1}, {'[H][H].[Pt](40)': 1}, (-10, 0), 'physisorption only; not stable on Ni(111) above ~50 K'),
        ('CO -> CO* (1 site)', {'CO(4)': 1, _S: 1}, {'OCX(16)': 1}, (-135, -110), 'Stuckless et al. 1993 calorimetry (~ -130)'),
        ('CO -> CO** (2 sites)', {'CO(4)': 1, _S: 2}, {'COXX(7400)': 1}, (-135, -110), 'same CO on Ni(111); bridge/hollow vs top differ by only a few kJ/mol'),
        ('1/2 O2 -> O*', {'O2(6)': 0.5, _S: 1}, {'OX(15)': 1}, (-260, -220), 'approx., from O2 adsorption heat ~ -480 kJ/mol O2'),
        ('NH3 -> NH3*', {'NH3(9)': 1, _S: 1}, {'[Pt].N(66)': 1}, (-80, -50), 'approx., TPD ~ 0.6-0.8 eV'),
        ('H2O -> H2O*', {'H2O(2)': 1, _S: 1}, {'O.[Pt](44)': 1}, (-60, -40), 'approx., weak molecular adsorption'),
        ('CH4 -> CH4*', {'CH4(5)': 1, _S: 1}, {'C.[Pt](57)': 1}, (-25, -10), 'physisorption'),
        ('CO2 -> CO2*', {'CO2(3)': 1, _S: 1}, {'CO2X(19)': 1}, (-40, -10), 'physisorption / very weak'),
    ],
    # switch_off: what-ifs that test how much the results depend on one species or route. Each
    # entry switches off reactions and solves the main conditions again (this costs solves): a
    # species name (all its surface reactions), a reaction number, or a list of those together.
    # Use it for a species the binding check flags TOO STRONG, or a route you suspect, e.g.
    # ['H(19)'] (H atoms desorbing to the gas). [] = no what-if.
    'switch_off': [],
}

# Optional (needs RMG-Py): flux-diagram videos of the RMG input's batch reactors, using the
# input.py and seed/seed in this folder. Run with:  bash run.sh runs/try_2/Ni --only flux
flux = {'rmg_path': '/Users/akinyemi.az/CoMoCheng/WESTLAB/rms_rmg_new/RMG-Py'}
