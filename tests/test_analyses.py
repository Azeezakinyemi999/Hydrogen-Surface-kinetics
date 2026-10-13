"""Every analysis end to end on the toy mechanism (small settings)."""
import os

import numpy as np
import pytest

import matplotlib
matplotlib.use('Agg')

import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import analyses
from hydrogen_surface_kinetics.cli import main as cli_main

from conftest import TOY_INPUT, write_input

SMALL = """
analyses = ['steady_state', 'composition', 'diagnostics', 'curves', 'sweep', 'orders',
            'rate_control', 'paths', 'transients']
plot_style = 'publication'
diagnostics = {'adsorption_checks': [('A2 -> 2 A*', {'A2': 1, 'S': 2}, {'A(s)': 2}, (-100, -60), 'made up')],
               'switch_off': 'I(ss)', 'top_reactions': [('I', 500.0, 'A(s)')]}
curves = {'ratios': [1e-6, 1e-5, 1e-4, 1e-3]}
sweep = {'temperatures': [450.0, 500.0, 550.0]}
rate_control = {'drc_top_n': 4}
transients = {'ratios': [1e-5, 1e-3], 't_start': 1e-6, 't_stop': 1e2, 'per_decade': 2}
"""


@pytest.fixture(scope='module')
def ctx(tmp_path_factory):
    path = write_input(tmp_path_factory.mktemp('toy'), TOY_INPUT, SMALL)
    return analyses.run(path, log=lambda *a: None)


def test_steady_state(ctx):
    s = ctx.outputs['steady_state']['summary']
    assert s['converged'].all()
    assert np.allclose(s['blocking_check'], 1.0, rtol=1e-6)
    ref = s[s['case'] == 'probe only'].iloc[0]
    assert ref['poisoning_factor'] == pytest.approx(1.0)


def test_files_written(ctx):
    out = ctx.cfg.output_dir
    for name in ('steady_state_results.csv', 'rates_summary.pdf', 'rates_summary.png', 'surface_composition.png',
                 'mechanism_binding_check.csv', 'mechanism_switch_off.csv', 'top_reactions.csv',
                 'poisoning_curves.csv', 'poisoning_curves.png', 'temperature_sweep.csv', 'reaction_orders.csv',
                 'sensitivity.csv', 'poisoning_over_time.csv', 'poisoning_over_time.png',
                 'binding_check.pdf', 'rate_control.pdf', 'sweep_coverages.pdf', 'poisoning_tolerance.csv'):
        assert os.path.exists(os.path.join(out, name)), name
    assert os.listdir(os.path.join(out, 'reaction_paths'))


def test_curves_monotonic(ctx):
    c = ctx.outputs['curves']['curves']
    c = c[c['case'] == 'I'].sort_values('ratio')
    assert (np.diff(c['poisoning_factor']) < 0).all()      # more impurity, slower


def test_probe_order_near_one_on_clean_surface(ctx):
    o = ctx.outputs['orders']['orders']
    # p(I) tiny, surface mostly empty: dissociative sticking is first order in p(A2), minus blocking
    assert 0 < o.loc[o['case'] == 'probe only', 'order_probe'].iloc[0] <= 1.0 + 1e-6


def test_drc_probe_reaction(ctx):
    d = ctx.outputs['rate_control']['sensitivity']
    d = d[(d['kind'] == 'DRC') & (d['item'] == 0)]
    assert d['X_tof_fwd'].notna().all()
    # the forward rate of the probe reaction scales with its own rate constant
    assert (d['X_tof_fwd'] > 0.5).all()


def test_transient_starts_at_reference(ctx):
    t = ctx.outputs['transients']['summary']
    ref = ctx.outputs['steady_state']['summary']
    tof_ref = ref.loc[ref['case'] == 'probe only', 'tof_fwd'].iloc[0]
    assert np.allclose(t['tof_t0'], tof_ref, rtol=1e-6)
    assert (t['end_over_t0'] <= 1.0 + 1e-6).all()


def test_rerun_loads_cache(ctx):
    lines = []
    analyses.run(ctx.cfg, only=['steady_state', 'curves'], log=lines.append)
    assert sum('loaded from cache, 0 to solve' in l for l in lines) == 2


def test_cli(ctx, capsys):
    assert cli_main(['check', ctx.cfg.input_file]) == 0
    assert cli_main(['status', ctx.cfg.input_file]) == 0
    assert 'cached results' in capsys.readouterr().out
    assert cli_main(['run', ctx.cfg.input_file, '--only', 'nonsense']) == 2


def test_unknown_analysis(tmp_path):
    with pytest.raises(hsk.ConfigError, match='unknown analyses'):
        hsk.load_config(write_input(tmp_path, TOY_INPUT, "analyses = ['steady', 'sweep']\n"))


def test_species_colors_follow_the_species(ctx):
    """A species keeps its surface-composition color in the coverage-vs-temperature figure."""
    from hydrogen_surface_kinetics.analyses.composition import species_palette
    base = ctx.outputs['composition']['colors']
    cov = ctx.outputs['sweep']['coverages']
    named, colors = species_palette(cov, ctx.model.site, ctx.cfg.composition, base)
    assert all(colors[sp] == base[sp] for sp in named if sp in base)
    assert len(set(colors.values())) == len(colors)            # no two species share a color
