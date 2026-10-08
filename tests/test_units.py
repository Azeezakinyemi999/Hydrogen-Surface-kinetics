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

from conftest import TOY_INPUT, write_input

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
    cfg = hsk.load_config(write_input(tmp_path, TOY_INPUT, "species_labels = {'I(s)': 'poison*'}\nprobe_label = 'A2'\n"))
    m = hsk.Model(cfg)
    assert species_label(m, 'I(s)') == 'poison*' and probe_label(m) == 'A2'


def test_site_density_units():
    assert site_density_kmol_m2((3.148e-9, 'mol/cm2')) == pytest.approx(3.148e-8)
    assert site_density_kmol_m2((3.148e-5, 'mol/m2')) == pytest.approx(3.148e-8)


def test_mixture_modes(tmp_path):
    m = hsk.Model(hsk.load_config(write_input(tmp_path, TOY_INPUT, TWO_IMPURITIES)))
    df = C.main_conditions(m)
    assert list(df['case']) == ['probe only', 'I', 'AI', 'both']
    both = df[df['case'] == 'both'].iloc[0]
    assert both['P_impurity_bar'] == pytest.approx(2e-5)          # same_each: each at the ratio
    m.cfg.mixture_mode = 'split'
    both = C.main_conditions(m).query("case == 'both'").iloc[0]
    assert both['P_impurity_bar'] == pytest.approx(1e-5)          # split: shared


def test_input_basis(tmp_path):
    m = hsk.Model(hsk.load_config(write_input(tmp_path, TOY_INPUT, TWO_IMPURITIES + "include_input_basis = True\n")))
    df = C.main_conditions(m)
    assert set(df['basis']) == {'fixed_probe', 'input_X'}
    x = df[(df['basis'] == 'input_X') & (df['case'] == 'I')].iloc[0]
    assert x['P_total_bar'] == pytest.approx(1.0) and x['X']['A2'] == pytest.approx(0.8)
    mix = df[(df['basis'] == 'input_X') & (df['case'] == 'both')].iloc[0]
    assert sum(mix['X'].values()) == pytest.approx(1.0)           # 0.75 + 2 x 0.05, renormalized


def test_input_pressure_mode(tmp_path):
    with pytest.raises(hsk.ConfigError, match='input_pressures'):
        hsk.load_config(write_input(tmp_path, TOY_INPUT, "pressure_mode = 'input'\ninput_pressures = {400: 2e5}\n"))
    cfg = hsk.load_config(write_input(tmp_path, TOY_INPUT, "pressure_mode = 'input'\ninput_pressures = {500: 2e5}\n"))
    df = C.main_conditions(hsk.Model(cfg))
    assert np.allclose(df['P_probe_bar'], 2.0)


def test_unknown_role_species(tmp_path):
    with pytest.raises(hsk.ModelError, match="probe_adsorbate 'B\\(s\\)'"):
        hsk.Model(hsk.load_config(write_input(tmp_path, TOY_INPUT, "probe_adsorbate = 'B(s)'\n")))


def test_unknown_section_key(tmp_path):
    with pytest.raises(hsk.ConfigError, match="unknown 'curves' settings"):
        hsk.load_config(write_input(tmp_path, TOY_INPUT, "curves = {'ratio': [0.1]}\n"))


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
