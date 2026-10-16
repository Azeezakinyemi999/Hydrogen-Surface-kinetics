"""The PDF report's building blocks (no solving; the full report is tested in test_analyses)."""
import os
import subprocess

import numpy as np
import pandas as pd
import pytest

import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import report

from conftest import TOY_SETTINGS, write_settings

AWKWARD = ['CH2O(83)', 'HX_1#2', '50% H2 & CO', '{x}~^\\', 'θ* → Π ≈ 1 · 10⁻³', 'C#O', 'site(10)', 'X<-Y|Z>']


def test_tex_escape_compiles(tmp_path):
    """Awkward species names and labels typeset without errors."""
    exe = report.latexmk_path({})
    if exe is None:
        pytest.skip('no LaTeX')
    body = '\n\n'.join(report.tex_escape(s) + ' / ' + report.tex_formula(s) for s in AWKWARD)
    (tmp_path / 't.tex').write_text('\\documentclass{article}\\usepackage[T1]{fontenc}\\usepackage{amsmath,amssymb}'
                                    '\\begin{document}\n' + body + '\n\\end{document}\n')
    pdflatex = os.path.join(os.path.dirname(exe), 'pdflatex')
    done = subprocess.run([pdflatex, '-interaction=nonstopmode', '-halt-on-error', 't.tex'], cwd=tmp_path,
                          capture_output=True, text=True, timeout=120)
    assert done.returncode == 0, done.stdout[-2000:]


def test_tex_escape_specials():
    assert report.tex_escape('a_b#c%d&e$f') == r'a\_b\#c\%d\&e\$f'
    assert report.tex_formula('H2O*') == r'H\textsubscript{2}O$\ast$'
    assert report.tex_formula('CH2O(83)') == r'CH\textsubscript{2}O(83)'          # the mechanism number stays
    assert report.tex_escape('é') == '?'                                       # unknown non-ASCII


def test_numbers():
    f = report.tex_num
    assert [f(v) for v in (500.0, 0.0999959, 0.99999, 1234.5, 0.00526)] == ['$500$', '$0.100$', '$1.00$', '$1230$', '$0.00526$']
    assert f(3.47e5) == '$3.47{\\times}10^{5}$' and f(-2.08e-10) == '$-2.08{\\times}10^{-10}$'
    assert f(np.nan) == '--' and f(None) == '--' and f(True) == 'yes' and f(0) == '$0$'
    assert report.tex_coverage(1e-16) == '$0$'
    assert report.tex_arrhenius(4.21e6, -2.06) == '$4.21{\\times}10^{6}\\,\\exp(+2.1/RT)$'
    assert report.tex_arrhenius(1.0, 1e-3) == '$1$ (no $T$ dependence)'
    assert report.fig_refs(['fig:3', 'fig:4', 'fig:5', 'fig:8']) == r'Figures \ref{fig:3}--\ref{fig:5}, \ref{fig:8}'
    assert report.fig_refs(['fig:2']) == r'Figure \ref{fig:2}'
    assert report.tex_units('mol H2 m-2 s-1 Pa-1') == r'mol H\textsubscript{2} m$^{-2}$ s$^{-1}$ Pa$^{-1}$'
    assert report.tex_units('-') == '--' and report.tex_units('m2 mol-1 s-1') == r'm$^{2}$ mol$^{-1}$ s$^{-1}$'


def test_table_borders_and_row_limit(toy_cfg):
    rep = report.Report(toy_cfg, log=lambda *a: None)
    df = pd.DataFrame({'case': ['a'] * 35 + ['b'] * 35, 'x': np.arange(70.0), 'note': [np.nan] * 70})
    tex = rep.table(df, [('case', 'case', 'l', report.tex_escape), ('x', '$x$', 'r', None), ('note', 'note', 'P3', str)],
                    'A table.', 'data.csv', group_by='case')
    assert '{|l|r|' in tex                                           # a rule between all columns
    rows = [line for line in tex.splitlines() if line.endswith((r'\\ \hline', r'\\ \thickhline')) and '&' in line]
    assert len(rows) == 60 + 2                                      # 60 shown (max_rows) + the header twice
    assert 'First 60 of 70 rows' in tex and r'data.csv' in tex
    assert ' nan' not in tex and '-- \\\\' in tex                      # missing values print as --
    full = rep.table(df, [('x', '$x$', 'r', None)], 'All.', full=True)
    assert full.count(r'\\ \hline') == 70


def test_wide_tables(toy_cfg):
    """A wrapping column is not narrowed below its header; a landscape table keeps its
    heading, and the landscape tables of one series share pages."""
    rep = report.Report(toy_cfg, log=lambda *a: None)
    df = pd.DataFrame({**{f'c{i}': ['x' * 17] for i in range(5)}, 'note': ['some words']})
    cols = [(f'c{i}', f'c{i}', 'l', str) for i in range(5)] + [('note', 'a very long header line of text', 'P5', str)]
    tex = rep.table(df, cols, 'Narrow.')
    assert tex.startswith(r'\begingroup\tiny') and 'p{5.00cm}' in tex      # not \scriptsize at 3.7 cm
    wide = pd.DataFrame([{'case': 'probe only', 'series': 'pressure', 'T_K': T, 'P_ref_bar': Pb,
                          **{f'v{i}': 'y' * 20 for i in range(14)}} for T in (400.0, 500.0) for Pb in (0.3, 3.0)])
    tables = rep.series_tables(wide, [(f'v{i}', f'v{i}', 'l', str) for i in range(14)], 'Wide', 'wide.csv')
    assert tables.count(r'\begin{landscape}') == 1 and tables.count(r'\end{landscape}') == 1
    assert tables.count(r'\begin{longtable}') == 2
    sec = rep.section('diagnostics', [('Wide tables', tables)])
    assert r'\begin{landscape}' + '\n' + r'\subsection*{Wide tables}' in sec
    assert sec.count(r'\subsection*{Wide tables}') == 1


def test_figure_parts_and_stems(toy_cfg):
    rep = report.Report(toy_cfg, log=lambda *a: None)
    os.makedirs(rep.out, exist_ok=True)
    for name in ('a.pdf', 'a.png', 'b_1.pdf', 'b_2.pdf', 'b_1.png', 'b_2.png', 'c_650K_DRC_1.png', 'c_650K_DRC_2.png',
                 'a_per_area.pdf', 'b_per_area_1.pdf', 'b_per_area_2.pdf'):
        open(os.path.join(rep.out, name), 'w').close()
    assert rep.figure_files('a') == ['a.pdf']                         # PDF preferred; not the per-area copy
    assert rep.figure_files('b') == ['b_1.pdf', 'b_2.pdf']            # parts, in order
    assert rep.figure_files('b_per_area') == ['b_per_area_1.pdf', 'b_per_area_2.pdf']
    assert rep.stems(r'.+') == ['a', 'a_per_area', 'b', 'b_per_area', 'c_650K_DRC']   # parts merged into their figure
    tex, labels = rep.figure('b', 'Caption.')
    assert labels == ['fig:1', 'fig:2'] and 'part 2 of 2' in tex


def test_path_diagram_names(toy_cfg):
    rep = report.Report(toy_cfg, log=lambda *a: None)
    folder = os.path.join(rep.out, 'reaction_paths')
    os.makedirs(folder)
    for case in ('I', 'probeonly'):
        for P in ('0.3', '1'):
            open(os.path.join(folder, f'path_H_{case}_500K_{P}bar_OneWayFlow.png'), 'w').close()
    found = report.path_diagrams(rep)
    assert [(d[1], d[3]) for d in found] == [('probe only', 1.0), ('I', 1.0), ('probe only', 0.3), ('I', 0.3)]
    assert report.paths_reference_point(rep, found) == (500.0, 1.0)  # at the reference pressure


def test_report_settings_checked(tmp_path):
    cfg = hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, "report = {'max_rows': 20}\n"))
    assert cfg.report['max_rows'] == 20 and cfg.report['enabled']     # merged with the defaults
    for bad, msg in (("report = {'path_diagrams': 'some'}\n", 'path_diagrams'),
                     ("report = {'max_rows': 0}\n", 'max_rows'), ("report = {'pages': 3}\n", 'unknown report')):
        with pytest.raises(hsk.ConfigError, match=msg):
            hsk.load_config(write_settings(tmp_path, TOY_SETTINGS, bad))


def test_pending_switch_off_note(toy_cfg):
    """A what-if listed in the settings but not run yet gets a diagnostics section with a note
    and the command that runs it."""
    toy_cfg.diagnostics['switch_off'] = ['I(ss)']
    rep = report.Report(toy_cfg, log=lambda *a: None)
    os.makedirs(rep.out, exist_ok=True)
    with open(report.write_tex(rep)) as f:
        tex = f.read()
    assert r'\label{sec:diagnostics}' in tex
    assert 'not run yet' in tex and '--only diagnostics' in tex
