"""Small pieces: labels, condition bases, site-density units, plot styles, error messages."""
import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import conditions as C, plotting as P
from hydrogen_surface_kinetics.config import site_density_kmol_m2
from hydrogen_surface_kinetics.labels import formula, probe_label, species_label

from conftest import TOY_SETTINGS, write_settings

TWO_IMPURITIES = """
impurities = {'I': 'I', 'AI': 'AI'}
mixture_label = 'both'
"""


def test_formula_and_labels(toy_model):
    assert formula({'C': 1, 'O': 2, 'X': 1}) == 'CO2'
    assert formula({'H': 4, 'C': 1}) == 'CH4'
    assert formula({'N': 1, 'H': 3, 'Pt': 1}) == 'H3N'
    assert probe_label(toy_model) == 'Ax2'
    assert species_label(toy_model, 'S') == 'empty site'
    assert species_label(toy_model, 'I(ss)') == 'Ix**'                  # '(ss)' is not a species number


def test_species_labels_setting(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "species_labels = {'I(s)': 'poison*'}\nprobe_label = 'A2'\n"))
    m = hsk.Model(cfg)
    assert species_label(m, 'I(s)') == 'poison*' and probe_label(m) == 'A2'


def test_site_density_units():
    assert site_density_kmol_m2((3.148e-9, 'mol/cm2')) == pytest.approx(3.148e-8)
    assert site_density_kmol_m2((3.148e-5, 'mol/m2')) == pytest.approx(3.148e-8)


def test_mixture_modes(tmp_path):
    m = hsk.Model(hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, TWO_IMPURITIES)))
    df = C.main_conditions(m)
    assert list(df['case']) == ['probe only', 'I', 'AI', 'both']
    both = df[df['case'] == 'both'].iloc[0]
    assert both['P_impurity_bar'] == pytest.approx(2e-5)          # same_each: each at the ratio
    m.cfg.mixture_mode = 'split'
    both = C.main_conditions(m).query("case == 'both'").iloc[0]
    assert both['P_impurity_bar'] == pytest.approx(1e-5)          # split: shared


def test_input_basis(tmp_path):
    m = hsk.Model(hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, TWO_IMPURITIES + "include_input_basis = True\n")))
    df = C.main_conditions(m)
    assert set(df['basis']) == {'fixed_probe', 'input_X'}
    x = df[(df['basis'] == 'input_X') & (df['case'] == 'I')].iloc[0]
    assert x['P_total_bar'] == pytest.approx(1.0) and x['X']['A2'] == pytest.approx(0.8)
    mix = df[(df['basis'] == 'input_X') & (df['case'] == 'both')].iloc[0]
    assert sum(mix['X'].values()) == pytest.approx(1.0)           # 0.75 + 2 x 0.05, renormalized


def test_removed_pressure_mode(tmp_path):
    with pytest.raises(hsk.ConfigError, match="'pressure_mode' removed"):
        hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "pressure_mode = 'input'\ninput_pressures = {500: 2e5}\n"))


def test_series_points_are_separate(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS,
                                         "temperatures = [450, 500, 550]\npressures = [3e4, 3e5]\n"
                                         "pressure_temperatures = [500, 600]\n"))
    pts = {(p['T'], p['P']): p['series'] for p in C.series_points(cfg)}
    # temperature series at the reference pressure; pressure series at each pressure temperature
    assert pts[(450.0, 1e5)] == 'temperature' and pts[(550.0, 1e5)] == 'temperature'
    assert pts[(500.0, 1e5)] == 'both'                           # the shared point is solved once
    assert pts[(600.0, 1e5)] == 'pressure'                       # reference pressure joins every P series
    assert pts[(500.0, 3e4)] == pts[(600.0, 3e5)] == 'pressure'
    assert (450.0, 3e4) not in pts and (550.0, 3e5) not in pts   # no T x P grid
    assert len(pts) == 3 + 2 * 3 - 1
    only_T = C.series_points(cfg, series='temperature')
    assert {p['P'] for p in only_T} == {1e5}


def test_labels_do_not_change_the_cache_key(tmp_path):
    from hydrogen_surface_kinetics import runner as R
    from hydrogen_surface_kinetics.analyses.common import job
    row = {'basis': 'b', 'case': 'c', 'T_K': 500.0, 'X': {'A2': 1.0}, 'P_ref_bar': 1.0}
    a = job('k', {**row, 'series': 'temperature'})
    b = job('k', {**row, 'series': 'both'})
    assert R.cache_path(str(tmp_path), a, ()) == R.cache_path(str(tmp_path), b, ())


def test_unknown_role_species(tmp_path):
    with pytest.raises(hsk.ModelError, match="probe_adsorbate 'B\\(s\\)'"):
        hsk.Model(hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "probe_adsorbate = 'B(s)'\n")))


def test_duplicate_flag_without_duplicate(tmp_path):
    from conftest import TOY
    text = open(TOY).read()
    text = text.replace('# R3, impurity spreads onto a second site\n',
                        '# R3, impurity spreads onto a second site\n  duplicate: true\n')       # no partner
    r4 = '  rate-constant: {A: 4.0e+21, b: 0.0, Ea: 85.0}\n'
    text = text.replace(r4, '  duplicate: true\n' + r4 + '- equation: AI + 2 S <=> A(s) + I(s)\n  duplicate: true\n'
                        '  rate-constant: {A: 1.0e+20, b: 0.0, Ea: 120.0}\n')                 # a real pair (one reversed)
    (tmp_path / 'mech.yaml').write_text(text)
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "mechanism = 'mech.yaml'\n"))
    m = hsk.Model(cfg)
    fixed = tmp_path / 'output' / 'mech.duplicates_fixed.yaml'
    assert m.mechanism_file == str(fixed) and m.surf.n_reactions == 5
    assert fixed.read_text().count('duplicate: true') == 2                                    # only R3's flag removed
    assert any('1 reactions in mech.yaml' in w for w in m.warnings)
    assert m.new_surface().n_reactions == 5
    mtime = fixed.stat().st_mtime_ns
    assert hsk.Model(cfg).mechanism_file == str(fixed) and fixed.stat().st_mtime_ns == mtime  # reused


def test_unknown_section_key(tmp_path):
    with pytest.raises(hsk.ConfigError, match="unknown 'curves' settings"):
        hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "curves = {'ratio': [0.1]}\n"))


def test_subscripts():
    assert P.sub('H2O and CO2') == 'H$_{2}$O and CO$_{2}$'
    assert P.sub('300 K') == '300 K'
    assert P.sub('$E_{app}$') == '$E_{app}$'


def test_publication_size_and_cap():
    fig, ax = plt.subplots(figsize=(15, 4.6))
    P.apply_style(fig, 'publication', 17.8)
    assert fig.get_size_inches()[0] == pytest.approx(17.8 / 2.54)
    tall, _ = plt.subplots(2, 1, figsize=(3.4, 5.6))
    P.apply_style(tall, 'publication', 17.8)
    assert tall.get_size_inches()[1] == pytest.approx(P.PUB_MAX_HEIGHT_CM / 2.54)
    plt.close('all')


def test_panel_letters_past_z(tmp_path):
    assert [P.panel_letter(i) for i in (0, 25, 26, 27, 52)] == ['a', 'z', 'aa', 'ab', 'ba']
    fig, axes = plt.subplots(6, 5)
    for ax in axes.flat:
        ax.set_title('CO, 650 K, 1 bar', loc='left')
    P.save(fig, str(tmp_path / 'many.png'), style='publication')         # 30 panels: past (z)
    assert axes[-1][-1].get_title(loc='left').startswith(r'$\bf{(ad)}$')
    plt.close(fig)


def test_save_writes_pdf_in_publication(tmp_path):
    fig, ax = plt.subplots()
    ax.plot([1, 2], [1, 2])
    files = P.save(fig, str(tmp_path / 'f.png'), style='publication')
    assert sorted(f.rsplit('.', 1)[1] for f in files) == ['pdf', 'png']
    plt.close(fig)


def test_flux_needs_paths(toy_cfg):
    from hydrogen_surface_kinetics import analyses
    with pytest.raises(ValueError, match="rmg_input"):
        analyses.run(toy_cfg, only='flux', log=lambda *a: None)


def test_rmg_input_is_refused(tmp_path):
    rmg = tmp_path / 'input.py'
    rmg.write_text("database(\n    thermoLibraries=['primaryThermoLibrary'],\n)\ncatalystProperties(metal='Ni111')\n")
    with pytest.raises(hsk.ConfigError, match='looks like an RMG input file'):
        hsk.load_config(str(rmg))


def test_flux_finds_rmg_files_next_to_settings(tmp_path):
    from hydrogen_surface_kinetics.analyses import flux
    from hydrogen_surface_kinetics.analyses.common import Context
    (tmp_path / 'input.py').write_text("database()\n")
    (tmp_path / 'seed' / 'seed').mkdir(parents=True)
    (tmp_path / 'seed' / 'seed' / 'reactions.py').write_text('')
    cfg = hsk.load_config(write_settings(tmp_path))
    ctx = Context(hsk.Model(cfg), log=lambda *a: None)
    rmg_input, seed_dir = flux.default_rmg_files(ctx)
    assert rmg_input == str(tmp_path / 'input.py') and seed_dir == str(tmp_path / 'seed' / 'seed')


def test_flux_reactor_temperatures():
    """A reactor with one T keeps it; RMG leaves T unset for a range, which is run at both ends
    by default or at flux['temperatures']."""
    from types import SimpleNamespace as NS
    from hydrogen_surface_kinetics.analyses import flux
    fixed = NS(T=NS(value_si=600.0), Trange=None)
    ranged = NS(T=None, Trange=[NS(value_si=650.0), NS(value_si=1200.0)])
    assert flux.reactor_temperatures(fixed) == [600.0]
    assert flux.reactor_temperatures(fixed, [700, 800]) == [600.0]
    assert flux.reactor_temperatures(ranged) == [650.0, 1200.0]
    assert flux.reactor_temperatures(ranged, [650, 800, 1000]) == [650.0, 800.0, 1000.0]


def test_per_impurity_ratios(tmp_path):
    m = hsk.Model(hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, TWO_IMPURITIES
                                                 + "impurity_ratio = {'I': 1 / 190, 'AI': 10 / 190}\n")))
    df = C.main_conditions(m).set_index('case')
    assert df.loc['I', 'X']['I'] == pytest.approx(1 / 191)            # each impurity alone at its own level
    assert df.loc['AI', 'X']['AI'] == pytest.approx(10 / 200)
    x = df.loc['both', 'X']                                           # the mixture is the reactor's gas
    assert (x['A2'], x['I'], x['AI']) == pytest.approx((190 / 201, 1 / 201, 10 / 201))
    assert df.loc['both', 'P_probe_bar'] == pytest.approx(1.0)        # H2 still at the point's pressure


def test_per_impurity_ratio_errors(tmp_path):
    with pytest.raises(hsk.ConfigError, match=r"no ratio for impurities \['AI'\]"):
        hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, TWO_IMPURITIES + "impurity_ratio = {'I': 0.1}\n"))
    with pytest.raises(hsk.ConfigError, match="mixture_mode = 'split' needs one number"):
        hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, TWO_IMPURITIES
                                       + "impurity_ratio = {'I': 0.1, 'AI': 0.2}\nmixture_mode = 'split'\n"))


def test_path_diagram_styles(tmp_path):
    from PIL import Image
    from hydrogen_surface_kinetics.analyses import paths
    m = hsk.Model(hsk.load_config(write_settings(tmp_path)))
    Image.new('RGB', (20, 40), 'white').save(tmp_path / 'A(s).png')          # an image for A(s) only
    raw = ('digraph reaction_paths {\ncenter=1;\n'
           's1 -> s2[fontname="Helvetica", penwidth=6, arrowsize=3, color="0.7, 1.5, 0.9"\n, label = " 1"];\n'
           's2 -> s1[fontname="Helvetica", penwidth=2, arrowsize=1, color="0.7, 0.5, 0.9"\n, label = " 1e-07"];\n'
           's1 [ fontname="Helvetica", label="A2"];\ns2 [ fontname="Helvetica", label="A(s)"];\n'
           ' label = "Scale = 1\\l t (OneWayFlow, width ~ flux, threshold 1e-08 of max; labels = flux / max flux)";\n'
           ' fontname = "Helvetica";\n}\n')
    rmg = paths.style_dot(m, raw, str(tmp_path), image_names=False, orientation='TB', edge_labels=True,
                          style='rmg', amounts={'A2': 1.0, 'A(s)': 0.5}, threshold=1e-8)
    assert 'IMG SRC="' + str(tmp_path / 'A(s).png') in rmg                  # image node for A(s)
    assert 'label="A2"' in rmg                                              # no image: readable name
    assert rmg.count('color="black"') >= 2 and '0.7, 1.5, 0.9' not in rmg    # RMG look: black arrows
    assert 'penwidth=8.00' in rmg and 'penwidth=1.00' in rmg                # log widths: 1 -> 8 pt, 1e-7 -> 1 pt
    assert 'label=" 1e-07"' in rmg                                          # numbers kept
    clean = paths.style_dot(m, raw, None, orientation='LR', edge_labels=False)
    assert 'rankdir="LR"' in clean and 'label = " 1"' not in clean and 'arrowsize=1.40' in clean


def test_arrhenius_fit_and_intrinsic_constants(tmp_path):
    import cantera as ct
    from hydrogen_surface_kinetics.analyses import arrhenius as AR
    from hydrogen_surface_kinetics.analyses.common import Context
    T = np.array([400.0, 500.0, 600.0, 800.0])
    f = AR.fit(T, 3e7 * np.exp(-60e6 / (ct.gas_constant * T)))      # A = 3e7 1/s, Ea = 60 kJ/mol
    assert f['A'] == pytest.approx(3e7, rel=1e-9) and f['Ea_kJ_per_mol'] == pytest.approx(60.0, rel=1e-9)
    assert f['R2'] == pytest.approx(1.0) and f['n_points'] == 4
    assert np.allclose(AR.local_ea(T, 3e7 * np.exp(-60e6 / (ct.gas_constant * T))), 60.0, rtol=0.02)
    assert np.isnan(AR.fit([500.0], [1.0])['Ea_kJ_per_mol'])         # one point: no fit
    cfg = hsk.load_config(write_settings(tmp_path))
    ctx = Context(hsk.Model(cfg), log=lambda *a: None)
    intr = AR.intrinsic_table(ctx, [450.0, 500.0, 600.0])
    m = ctx.model
    for row in intr.itertuples():                                    # K = Kc P0 / RT, from Cantera itself
        surf = m.new_surface()
        m.gas.TP = row.T_K, cfg.reference_pressure
        surf.TP = row.T_K, cfg.reference_pressure
        Kc = surf.equilibrium_constants[m.i_probe]
        assert row.K_eq == pytest.approx(Kc * AR.P0 / (ct.gas_constant * row.T_K), rel=1e-9)
    vh = AR.fit(intr['T_K'], intr['K_eq'])                           # van 't Hoff dH vs the mechanism's thermo
    assert vh['Ea_kJ_per_mol'] == pytest.approx(intr['dH_thermo_kJ_per_mol'].mean(), abs=3.0)


def test_per_area_conversions(toy_model):
    """Format 2 from Format 1: J = TOF x Gamma, kd_eff = J/p, gamma_eff = J/Z, C_H = theta_H Gamma."""
    import pandas as pd
    from hydrogen_surface_kinetics import per_area
    G = per_area.site_density_mol_m2(toy_model)
    assert np.isclose(G, 2.5e-5)                                    # 2.5e-9 mol/cm2
    M = per_area.molar_mass_kg_mol(toy_model)
    Z = 1e5 / np.sqrt(2 * np.pi * M * 8.314462618 * 500.0)          # by hand
    assert np.isclose(per_area.impingement_flux(500.0, 1e5, M), Z)
    df = pd.DataFrame({'T_K': [500.0, 500.0], 'P_probe_bar': [1.0, 1.0], 'tof_fwd': [10.0, np.nan],
                       'tof_rev': [4.0, np.nan], 'tof_net': [6.0, np.nan], 'theta_probe': [0.5, np.nan]})
    out = per_area.add_columns(df, toy_model).iloc[0]
    assert np.isclose(out['J_diss_mol_m2_s'], 10 * G) and np.isclose(out['J_net_mol_m2_s'], 6 * G)
    assert np.isclose(out['J_exchange_mol_m2_s'], 4 * G)            # min(J_diss, J_rec)
    assert np.isclose(out['kd_eff_mol_m2_s_Pa'], 10 * G / 1e5) and np.isclose(out['gamma_eff'], 10 * G / Z)
    assert np.isclose(out['C_H_mol_m2'], 0.5 * G)
    assert per_area.add_columns(df, toy_model).iloc[1][list(per_area.COLUMNS)].isna().all()


def test_clean_surface_against_the_mechanism(toy_model):
    """On an empty surface the sticking probability is the mechanism's (A = 0.05, b = 0, Ea = 0);
    K agrees with the thermochemistry (converted from its reference pressure to P0 = 1 bar),
    the equilibrium coverage satisfies K, and
    J0 = kd p theta_v^2 = kr (theta_H Gamma)^2."""
    from hydrogen_surface_kinetics import per_area
    from hydrogen_surface_kinetics.analyses.diagnostics import reaction_thermo
    G = per_area.site_density_mol_m2(toy_model)
    for T, p in ((450.0, 1e5), (500.0, 3e4), (550.0, 1e6)):
        c = per_area.clean_surface(toy_model, T, p)
        assert np.isclose(c['gamma_clean'], 0.05, rtol=1e-9)
        _, dG = reaction_thermo(toy_model, {'A2': 1, 'S': 2}, {'A(s)': 2}, T)
        P_thermo = toy_model.gas.species('A2').thermo.reference_pressure      # 1 atm here; K uses 1 bar
        assert np.isclose(c['K_eq'], np.exp(-dG * 1e3 / (8.314462618 * T)) * 1e5 / P_thermo, rtol=1e-6)
        assert np.isclose(c['theta_H_eq'] ** 2 / (c['theta_v_eq'] ** 2 * p / 1e5), c['K_eq'])
        assert np.isclose(c['theta_H_eq'] + c['theta_v_eq'], 1.0)
        assert np.isclose(c['J0_clean_mol_m2_s'], c['kr_surface_m2_mol_s'] * (c['theta_H_eq'] * G) ** 2)
        assert np.isclose(c['kd_clean_mol_m2_s_Pa'] / c['kr_surface_m2_mol_s'], G ** 2 * c['K_eq'] / 1e5)


def test_switch_off_settings(tmp_path, toy_model):
    """switch_off: [] by default (no what-if); species, reaction numbers or groups; clear errors."""
    from hydrogen_surface_kinetics.analyses import diagnostics as D
    from hydrogen_surface_kinetics.analyses.common import Context
    from hydrogen_surface_kinetics.analyses.plan import counts
    cfg = toy_model.cfg
    assert cfg.diagnostics['switch_off'] == [] and D.switch_off_items(cfg) == []
    row = counts(Context(toy_model, log=lambda *a: None), ['diagnostics']).set_index('analysis').loc['diagnostics']
    assert row['to_solve'] == 0 and 'no switch-off' in row['note']
    cfg.diagnostics['switch_off'] = ['I(ss)', 3, ['I(s)', 3]]
    items = D.switch_off_items(cfg)
    assert [label for label, _ in items] == ['I(ss)', 'R3', 'I(s)+R3']
    assert D.switch_off_reactions(toy_model, ['I(ss)']) == [2]
    assert D.switch_off_reactions(toy_model, ['I(s)', 3]) == [1, 2, 3]
    assert D.file_label('H2_X(320)') == 'H2_X_320' and D.file_label('I(s)+R3') == 'I_s+R3'
    with pytest.raises(ValueError, match='does not exist'):
        D.switch_off_reactions(toy_model, [99])
    with pytest.raises(ValueError, match='not a species'):
        D.switch_off_reactions(toy_model, ['XYZ'])
    for bad in ("{'a': 1}", '2.5', '[[]]'):
        with pytest.raises(hsk.ConfigError, match='switch_off'):
            hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, f"diagnostics = {{'switch_off': {bad}}}\n"))
