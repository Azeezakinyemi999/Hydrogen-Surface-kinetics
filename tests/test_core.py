"""Config, model, conditions, solver and runner on the toy mechanism."""
import os

import cantera as ct
import numpy as np
import pytest

import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import conditions as C, solver as S, runner as R

from conftest import TOY_SETTINGS, write_settings


# --- config ---------------------------------------------------------------------------
def test_defaults_and_merge(toy_cfg, tmp_path):
    assert toy_cfg.solver['abs_tol'] == 1e-8          # default kept
    assert toy_cfg.run['workers'] == 2                # override merged
    assert toy_cfg.run['job_timeout'] == 1800.0       # rest of the section kept
    assert toy_cfg.output_dir == os.path.join(str(tmp_path), 'output')


def test_config_lists_every_problem(tmp_path):
    path = write_settings(tmp_path, "mechanism = 'nope.yaml'\nsolver = {'abstol': 1}\nmixture_mode = 'x'\n")
    with pytest.raises(hsk.ConfigError) as e:
        hsk.load_config(path)
    msg = str(e.value)
    for text in ("'probe_gas' is required", 'mechanism file not found', 'abstol', 'mixture_mode'):
        assert text in msg


def test_keyword_overrides(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path), temperatures=[450, 550], run={'recompute': True})
    assert cfg.temperatures == [450.0, 550.0]
    assert cfg.run['recompute'] and cfg.run['workers'] == 2   # merged into the settings file's section


# --- model ----------------------------------------------------------------------------
def test_model_roles(toy_model):
    m = toy_model
    assert m.i_probe == 0 and m.sign_probe == 1 and m.probe_sites == 2
    assert m.impurities == {'I': 'I'}
    assert np.isclose(m.gamma, 2.5e-8)


def test_wrong_site_density(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "expected_site_density = (2.5e-8, 'mol/cm2')\n"))
    with pytest.raises(hsk.ModelError, match='site density'):
        hsk.Model(cfg)


def test_impurity_by_composition(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "impurities = {'imp': {'Ix': 1}}\n"))
    assert hsk.Model(cfg).impurities == {'imp': 'I'}


def test_enthalpy_shift(toy_model):
    m, dh = toy_model, 1e6
    surf0 = m.new_surface()
    h0 = surf0.partial_molar_enthalpies[surf0.species_index('A(s)')]
    surf1 = m.new_surface(h_shifts={'A(s)': dh})
    assert np.isclose(surf1.partial_molar_enthalpies[surf1.species_index('A(s)')] - h0, dh)


# --- conditions -----------------------------------------------------------------------
def test_conditions(toy_model):
    df = C.main_conditions(toy_model)
    assert list(df['case']) == ['probe only', 'I']          # one impurity: no mixture case
    row = df.iloc[1]
    assert np.isclose(row['P_probe_bar'], 1.0) and np.isclose(row['P_impurity_bar'], 1e-5)
    assert C.condition_key(row) == 'fixed_probe | I | 500 K'


# --- solver ---------------------------------------------------------------------------
def test_solve_and_residual(toy_model):
    m = toy_model
    row = C.main_conditions(m).iloc[1].to_dict()
    res = S.solve_condition(m, row)
    assert res['info_a']['converged'] and res['info_b']['converged']
    assert res['starts_agree']
    theta = res['theta_a']
    assert np.isclose(theta.sum(), 1.0)
    # Re-evaluate the solved state on a fresh object: every species at steady state
    S.set_state(m, row['T_K'], row['P_total_bar'] * 1e5, row['X'], S.coverages_to_dict(m, theta))
    assert np.abs(m.surf.net_production_rates[m.surf_kin_idx] * m.site_size / m.gamma).max() < 1e-6


def test_direct_and_reactor_agree(toy_model):
    m = toy_model
    row = C.main_conditions(m).iloc[1].to_dict()
    args = (row['T_K'], row['P_total_bar'] * 1e5, row['X'], S.clean_surface(m))
    th_d, info_d = S.solve_steady_state(m, *args, methods=('direct',))
    th_r, info_r = S.solve_steady_state(m, *args, methods=('reactor',))
    assert info_d['converged'] and info_r['converged']
    assert np.abs(th_d - th_r).max() < 1e-4


def test_multiplier_scales_both_directions(toy_model):
    m = toy_model
    row = C.main_conditions(m).iloc[0].to_dict()
    S.set_state(m, row['T_K'], 1e5, row['X'], {'S': 0.5, 'A(s)': 0.5})
    f0, r0, _ = m.probe_rates()
    S.set_state(m, row['T_K'], 1e5, row['X'], {'S': 0.5, 'A(s)': 0.5}, multipliers={0: 2.0})
    f1, r1, _ = m.probe_rates()
    assert np.isclose(f1, 2 * f0) and np.isclose(r1, 2 * r0)


# --- runner ---------------------------------------------------------------------------
def test_run_conditions_and_cache(toy_model):
    m = toy_model
    cond = C.main_conditions(m)
    first = S.run_conditions(m, cond)
    assert all(r['status'] == 'done' for r in first.values())
    files = os.listdir(os.path.join(m.cfg.output_dir, 'cache'))
    assert len(files) == len(cond)
    second = S.run_conditions(m, cond)                   # all from cache
    for k in first:
        assert np.allclose(first[k]['theta_a'], second[k]['theta_a'])
    status = S.status_table(second)
    assert status['converged_A'].all() and status['starts_agree'].all()


def test_cache_key_changes_with_settings(tmp_path):
    job = {'key': 'a | b', 'row': {'T_K': 500.0, 'X': {'A2': 1.0}}}
    p1 = R.cache_path(str(tmp_path), job, ('ctx', 1e-8))
    p2 = R.cache_path(str(tmp_path), job, ('ctx', 1e-9))
    p3 = R.cache_path(str(tmp_path), {**job, 'row': {'X': {'A2': 1.0}, 'T_K': 500.0}}, ('ctx', 1e-8))
    assert p1 != p2 and p1 == p3


def test_runner_timeout_and_error(tmp_path):
    import time

    def func(job):
        if job['kind'] == 'slow':
            time.sleep(30)
        if job['kind'] == 'bad':
            raise RuntimeError('boom')
        return {'value': 1}

    jobs = [{'key': k, 'kind': k} for k in ('ok', 'slow', 'bad')]
    res = R.run_jobs(func, jobs, str(tmp_path), workers=3, timeout=2, log=lambda *a: None)
    assert res['ok']['status'] == 'done' and res['ok']['value'] == 1
    assert res['slow']['status'] == 'timeout'
    assert res['bad']['status'] == 'error' and 'boom' in res['bad']['error']
