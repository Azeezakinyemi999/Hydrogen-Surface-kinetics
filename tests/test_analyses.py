"""Every analysis end to end on the toy mechanism (small settings)."""
import os

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')

import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import analyses
from hydrogen_surface_kinetics.cli import main as cli_main

from conftest import TOY_SETTINGS, write_settings

SMALL = """
pressures = [3e4, 3e5]                # pressure series at 500 K: 0.3, 1 (shared) and 3 bar
analyses = ['steady_state', 'composition', 'diagnostics', 'convergence', 'curves', 'sweep', 'arrhenius', 'pressure_sweep',
            'orders', 'rate_control', 'paths', 'transients']
plot_style = 'publication'
diagnostics = {'adsorption_checks': [('A2 -> 2 A*', {'A2': 1, 'S': 2}, {'A(s)': 2}, (-100, -60), 'made up')],
               'switch_off': 'I(ss)', 'top_reactions': [('I', 500.0, 'A(s)')]}
curves = {'ratios': [1e-6, 1e-5, 1e-4, 1.585e-4, 1e-3]}   # 1.585e-4: the toy has two steady states there
convergence = {'conditions': [('probe only', 500.0), ('I', 500.0, 3e5)]}
pressure_sweep = {'pressures': [3e4, 1e5, 3e5, 1e6]}
sweep = {'temperatures': [450.0, 500.0, 550.0]}
rate_control = {'drc_top_n': 4}
transients = {'ratios': [1e-5, 1e-3], 't_start': 1e-6, 't_stop': 1e2, 'per_decade': 2}
"""


@pytest.fixture(scope='module')
def ctx(tmp_path_factory):
    path = write_settings(tmp_path_factory.mktemp('toy'), TOY_SETTINGS, SMALL)
    return analyses.run(path, log=lambda *a: None)


def test_steady_state(ctx):
    s = ctx.outputs['steady_state']['summary']
    assert s['converged'].all()
    assert np.allclose(s['blocking_check'], 1.0, rtol=1e-6)
    ref = s[s['case'] == 'probe only']
    assert np.allclose(ref['poisoning_factor'], 1.0)              # every point has its own reference
    assert set(s['series']) == {'both', 'pressure'}               # 500 K / 1 bar is in both series
    assert sorted(s['P_ref_bar'].unique()) == [0.3, 1.0, 3.0]


def test_files_written(ctx):
    out = ctx.cfg.output_dir
    names = ['steady_state_results.csv', 'mechanism_binding_check.csv', 'mechanism_switch_off.csv',
             'top_reactions.csv', 'poisoning_curves.csv', 'temperature_sweep.csv', 'reaction_orders.csv',
             'sensitivity.csv', 'poisoning_over_time.csv', 'binding_check.pdf', 'sweep_coverages.pdf',
             'poisoning_tolerance.csv', 'pressure_sweep.pdf', 'pressure_sweep.csv', 'half_rate_pressure.csv',
             'pressure_sweep_coverages_500K.pdf', 'convergence.pdf', 'convergence.csv', 'site_blocking.pdf',
             'rate_control_500K_0.3bar.pdf', 'rate_control_500K_1bar.pdf', 'rate_control_500K_3bar.pdf',
             'poisoning_over_time_500K_1bar.pdf',          # temperature series: one point, one per point
             'poisoning_over_time_I_vs_P.pdf']             # pressure series: more points than impurities
    for stem in ('rates_summary', 'surface_composition', 'poisoning_curves', 'reaction_orders'):
        names += [f'{stem}_vs_T.pdf', f'{stem}_vs_P.pdf']                # one figure per series
    for name in names:
        assert os.path.exists(os.path.join(out, name)), name
    assert os.listdir(os.path.join(out, 'reaction_paths'))


def test_figures_split_at_nine_panels(tmp_path):
    import matplotlib.pyplot as plt
    from hydrogen_surface_kinetics import plotting as P
    assert [len(p) for p in P.pages(range(9), 4)] == [3, 3, 3] and P.pages([], 4) == []
    assert [P.grid_shape(n) for n in (3, 5, 7, 9)] == [(1, 3), (2, 3), (2, 4), (3, 3)]
    assert P.rows_per_page(3) == 3 and P.rows_per_page(4) == 2
    path = str(tmp_path / 'f.png')
    P.save_parts([plt.figure(), plt.figure()], path, style='screen')
    assert sorted(os.listdir(tmp_path)) == ['f_1.png', 'f_2.png']
    P.save_parts([plt.figure()], path, style='screen')                  # fewer parts: old parts removed
    assert sorted(os.listdir(tmp_path)) == ['f.png']


def test_curves_monotonic(ctx):
    c = ctx.outputs['curves']['curves']
    for _, part in c[c['case'] == 'I'].groupby('P_ref_bar'):
        part = part.sort_values('ratio')
        assert (np.diff(part['poisoning_factor']) < 0).all()      # more impurity, slower, at every pressure


def test_probe_order_near_one_on_clean_surface(ctx):
    o = ctx.outputs['orders']['orders']
    # p(I) tiny, surface mostly empty: dissociative sticking is first order in p(A2), minus blocking
    probe = o.loc[o['case'] == 'probe only', 'order_probe']
    assert len(probe) == 3 and ((probe > 0) & (probe <= 1.0 + 1e-6)).all()


def test_drc_probe_reaction(ctx):
    d = ctx.outputs['rate_control']['sensitivity']
    d = d[(d['kind'] == 'DRC') & (d['item'] == 0)]
    assert d['X_tof_fwd'].notna().all()
    # the forward rate of the probe reaction scales with its own rate constant
    assert (d['X_tof_fwd'] > 0.5).all()


def test_transient_starts_at_reference(ctx):
    t = ctx.outputs['transients']['summary']
    ref = ctx.outputs['steady_state']['summary']
    ref = ref[ref['case'] == 'probe only'].set_index('P_ref_bar')['tof_fwd']
    assert np.allclose(t['tof_t0'], t['P_ref_bar'].map(ref), rtol=1e-6)   # each starts from its own point
    assert (t['end_over_t0'] <= 1.0 + 1e-6).all()


def test_rerun_loads_cache(ctx):
    lines = []
    analyses.run(ctx.cfg, only=['steady_state', 'curves'], log=lines.append)
    assert sum('loaded from cache, 0 to solve' in l for l in lines) == 2


def test_cli(ctx, capsys):
    assert cli_main(['check', ctx.cfg.settings_file]) == 0
    assert cli_main(['status', ctx.cfg.settings_file]) == 0
    assert 'cached results' in capsys.readouterr().out
    assert cli_main(['run', ctx.cfg.settings_file, '--only', 'nonsense']) == 2


def test_unknown_analysis(tmp_path):
    with pytest.raises(hsk.ConfigError, match='unknown analyses'):
        hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "analyses = ['steady', 'sweep']\n"))


def test_species_colors_follow_the_species(ctx):
    """A species keeps its surface-composition color in the coverage-vs-temperature figure."""
    from hydrogen_surface_kinetics.analyses.composition import species_palette
    base = ctx.outputs['composition']['colors']
    cov = ctx.outputs['sweep']['coverages']
    named, colors = species_palette(cov, ctx.model.site, ctx.cfg.composition, base)
    assert all(colors[sp] == base[sp] for sp in named if sp in base)
    assert len(set(colors.values())) == len(colors)            # no two species share a color


def test_pressure_sweep(ctx):
    ps = ctx.outputs['pressure_sweep']['sweep']
    assert set(ps['P_ref_bar'].round(6)) == {0.3, 1.0, 3.0, 10.0} and set(ps['T_K']) == {500.0}
    ref = ps[ps['case'] == 'probe only']
    assert ref['converged'].all() and np.allclose(ref['poisoning_factor'], 1.0)
    slope = ref['slope_dlnTOF_dlnp'].dropna()
    assert len(slope) == 4 and ((slope > 0) & (slope <= 1.0 + 1e-6)).all()   # order in H2 of a 2-site step
    imp = ps[ps['case'] == 'I']                                    # fixed purity: impurity scales with H2
    assert np.allclose(imp['P_impurity_bar'] / imp['P_probe_bar'], 1e-5)


def test_convergence(ctx):
    out = ctx.outputs['convergence']
    for res in out.values():
        assert res['status'] == 'done' and res['accepted_s'] is not None
        traj = res['trajectory']
        at = traj[traj['t_s'] == res['accepted_s']].iloc[0]
        assert at['residual_ok'] and at['plateau_ok'] and at['norm_residual'] < 1
        assert traj['norm_residual'].iloc[0] > 1                    # far from steady state at 1e-8 s


def test_two_branches_and_site_blocking(ctx):
    c = ctx.outputs['curves']['curves']
    two = c[c['starts_agree'] == False]                              # noqa: E712
    assert len(two), 'the toy should have two steady states at ratio 1.585e-4'
    assert (np.abs(np.log(two['poisoning_factor_B'] / two['poisoning_factor'])) > 0.1).all()
    ok = c[c['converged'] & (c['case'] == 'I')]
    assert np.allclose(ok['poisoning_factor'], ok['blocking_full'], rtol=1e-6)   # elementary step: blocking only


def test_arrhenius(ctx):
    fits = ctx.outputs['arrhenius']['fits']
    assert set(fits['level']) == {'intrinsic', 'observed'}
    obs = fits[fits['level'] == 'observed']
    assert set(obs['case']) == {'probe only', 'I'} and (obs['n_points'] == 3).all()
    ref = obs[obs['case'] == 'probe only']
    assert np.allclose(ref['dEa_vs_reference_kJ_per_mol'], 0.0) and np.allclose(ref['A_over_reference'], 1.0)
    assert all(c for c in obs['conditions'])                        # every fit states its conditions
    out = ctx.cfg.output_dir
    for name in ('arrhenius_intrinsic.csv', 'arrhenius_fits.csv', 'arrhenius_local.csv', 'arrhenius.pdf'):
        assert os.path.exists(os.path.join(out, name)), name
