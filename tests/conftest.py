import os
import textwrap

import pytest

import hydrogen_surface_kinetics as hsk

TOY = os.path.join(os.path.dirname(__file__), 'data', 'toy_mechanism.yaml')

TOY_INPUT = f"""
mechanism = {TOY!r}
expected_site_density = (2.5e-9, 'mol/cm2')
probe_gas = 'A2'
empty_site = 'S'
probe_adsorbate = 'A(s)'
impurities = {{'I': 'I'}}
temperatures = [500.0]
reference_pressure = 1.0e5
impurity_ratio = 1e-5
run = {{'workers': 2}}
"""


def write_input(tmp_path, text=TOY_INPUT, extra=''):
    path = tmp_path / 'input.py'
    path.write_text(textwrap.dedent(text) + textwrap.dedent(extra))
    return str(path)


@pytest.fixture
def toy_cfg(tmp_path):
    return hsk.load_config(write_input(tmp_path))


@pytest.fixture
def toy_model(toy_cfg):
    return hsk.Model(toy_cfg)
