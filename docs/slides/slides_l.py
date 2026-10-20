"""Section L (slides 59-61): running it and reading it: run.sh and its four commands, the output
folder, and the PDF report. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PL = 'L · Running it and reading it'
RUN = 'runs/try_2/Cu'


def s59_commands(d):
    sl = d.slide('Running it: one script, four commands',
                 'The run folders come first, then the command (run if none is given), then the options.',
                 part=PL, tag='package · run.sh, cli.py')
    rows = (('bash run.sh runs/try_2/Cu', 'every analysis listed in settings.py, in the package’s order, then the PDF '
                                          'report'),
            ('bash run.sh runs/try_2/Cu check', 'loads the settings and the mechanism, lists the 72 conditions and the '
                                                'solves per analysis; solves nothing'),
            ('bash run.sh runs/try_2/Cu status', 'counts the saved results in each cache and those not converged (none '
                                                 'for try_2/Cu); lists the tables and figures'),
            ('bash run.sh runs/try_2/Cu report', 'rebuilds the PDF from the saved results; solves nothing '
                                                 '(--tex-only writes the .tex only)'))
    for i, (cmd, what) in enumerate(rows):
        y = 1.7 + i * 0.69
        shape(sl, 'rrect', 0.55, y, 4.4, 0.6, fill='codebg', line='codeline', lw=0.75, radius=0.06, content=cmd,
              size=13, font=CODE, align='l', margin=0.12)
        text(sl, 5.15, y, 7.6, 0.6, what, size=14, anchor='m')
    card(sl, 0.55, 4.55, 7.25, 1.33, 'Options of a run')
    text(sl, 0.75, 4.95, 6.95, 0.9, ['`--only curves,sweep`: only those, in the package’s order; no report',
                                     '`--workers 8` and `--timeout 1800`: processes; time limit per job (s)',
                                     '`--recompute`: ignore the cache;  `--no-report`: skip the PDF'], size=13)
    card(sl, 8.0, 4.55, 4.78, 1.33, 'Several runs in one call')
    text(sl, 8.2, 4.95, 4.45, 0.3, '`bash run.sh runs/try_2/*`', size=13)
    text(sl, 8.2, 5.27, 4.45, 0.58, 'One after another; a failed run does not stop the rest. Each appends to its '
         'own run.log.', size=13)
    actions(sl, ['Check first: the conditions and the solve counts.',
                 'Run: every analysis, then the report.',
                 'Status: what is cached and what did not converge.',
                 'Report: rebuild the PDF without solving.'])
    notes(sl, ['run.sh: usage lines 2-19; run folders or settings files first (26-36: a folder needs settings.py, so an RMG',
               'input.py is refused), then check, status or report (38-41, default run), then the options for every run.',
               'PYTHON defaults to the rmg4.0 env (43-44); one thread per worker (47); caffeinate -dims for run only (53-56,',
               'slide 38); Ctrl-C stops every run (58-60); each run is teed to <run folder>/run.log (64-80); "N of M runs OK"',
               'with several runs (82-87); the exit status is that of the last failed run.',
               'cli.py: cmd_check (lines 26-59: the conditions and analyses.plan counts; nothing solved), cmd_run (62-69:',
               'unknown analysis names refused), cmd_report (72-74), cmd_status (77-97), main (100-125: --only, --workers,',
               '--timeout, --recompute, --no-report; report --tex-only). status on try_2/Cu (2026-10-20, read-only, 8 min',
               'to open 4,679 cached files): 0 not converged in each of the 5 caches.',
               'analyses/__init__.py run (lines 14-40): analyses always run in ANALYSIS_ORDER; the report only after the',
               'full list (line 37), so --only skips it; build_report (43-53) logs a failed report instead of raising.',
               'Docs: "Running an analysis" (line 1159).'])


def s60_output(d):
    sl = d.slide('The output folder: tables, figures and caches',
                 'Everything one run writes, here for try_2/Cu: 94 MB. Each solve is saved under its fingerprint (slide 38).',
                 part=PL, tag=f'{RUN}/output/')
    card(sl, 0.55, 1.66, 6.3, 4.22, 'What is in it')
    items = (('*.csv', '21 tables, e.g. steady_state_results.csv'),
             ('*.png, *.pdf', '61 figures, each in both formats'),
             ('reaction_paths/', '72 path diagrams, .dot and .png (23 MB)'),
             ('report/', 'the PDF, its .tex and the settings listing (21 MB)'),
             ('cache*/', '5 folders of saved solves (24 MB)'))
    for i, (name, what) in enumerate(items):
        y = 2.12 + i * 0.6
        shape(sl, 'rrect', 0.75, y, 2.1, 0.42, fill='white', line='faint', radius=0.06, content=name, size=13,
              font=CODE, bold=True, align='l', margin=0.1)
        text(sl, 3.0, y, 3.75, 0.42, what, size=13, anchor='m')
    text(sl, 0.75, 5.18, 5.95, 0.6, 'Next to output/: run.log, everything printed by every command, appended.',
         size=13, color='muted')
    card(sl, 7.05, 1.66, 5.73, 4.22, 'The caches: one file per solve')
    text(sl, 7.2, 2.08, 1.9, 0.3, 'folder', size=12, color='muted', bold=True)
    text(sl, 9.05, 2.08, 0.8, 0.3, 'files', size=12, color='muted', bold=True, align='r')
    text(sl, 10.05, 2.08, 2.6, 0.3, 'from', size=12, color='muted', bold=True)
    caches = (('cache/', '1,326', 'main, curves, sweeps, orders'), ('cache_sensitivity/', '3,056', 'rate control'),
              ('cache_transients/', '224', 'poisoning over time'), ('cache_convergence/', '1', 'convergence traces'),
              ('cache_diagnostics/', '72', 'switch-off what-ifs'))
    for i, (name, n, what) in enumerate(caches):
        y = 2.42 + i * 0.4
        line(sl, 7.2, y - 0.02, 12.6, y - 0.02, color='rule', lw=0.75)
        text(sl, 7.2, y, 1.95, 0.36, name, size=12, font=CODE, anchor='m', color='amber' if i == 4 else 'ink')
        text(sl, 9.05, y, 0.8, 0.36, n, size=13, bold=True, align='r', anchor='m')
        text(sl, 10.05, y, 2.65, 0.36, what, size=13, anchor='m')
    caveat(sl, 7.2, 4.5, 5.43, 1.25, 'cache_diagnostics/ holds 72 results of an H₂∗ switch-off from Oct 13, the old '
           'default. Nothing reads them while switch_off = []; they would be reused if that what-if came back.',
           size=13)
    actions(sl, ['Read the numbers in the CSV tables.',
                 'Look at the figures, each as PNG and PDF.',
                 'Open the report for the whole run.',
                 'Keep the caches: they make every re-run fast.'])
    notes(sl, [f'Counts and sizes of {RUN}/output (ls and du, 2026-10-20): 21 .csv, 61 .png and 61 .pdf in the folder itself;',
               'reaction_paths/ 144 files (72 diagrams x .dot, .png), 23 MB; report/ try_2_Cu_report.pdf, .tex and',
               'settings_listing.py, 21 MB; caches 9.1 MB + 12 MB + 2.6 MB + 0.3 MB + 24 kB; 94 MB in all.',
               'cache/ 1,326 = 72 main + 680 curves + 108 sweep + 162 pressure sweep + 304 orders; cache_sensitivity/',
               '3,056 = 64 base + 2,992 perturbed; cache_transients/ 224; cache_convergence/ 1 trace.',
               'cache_diagnostics/: 72 files switch_off_H_H_._Pt_34_<case>_<T>_K_<P>_bar_<hash>.pkl, all dated 2026-10-13',
               '20:34, from the old default that switched off the molecular adsorbate [H][H].[Pt](34); settings.py now has',
               "'switch_off': [] (line 75), so diagnostics solves nothing and the files are unused (a cache hit if the",
               "what-if ['[H][H].[Pt](34)'] is set again with the same mechanism and solver settings).",
               'File names: runner.py cache_path (lines 42-45). Docs: Appendix C, Output files (line 1600).'])


def s61_report(d):
    sl = d.slide('The PDF report: one section per analysis, written to be read',
                 'Built from the saved results at the end of a full run, or by report. Nothing is solved to make it.',
                 part=PL, tag='package · report.py, report_text.py')
    card(sl, 0.55, 1.66, 5.6, 4.22, 'try_2_Cu_report.pdf: 133 pages')
    text(sl, 0.75, 2.08, 5.3, 0.3, 'Title page: the run, the mechanism’s SHA-256, the package commit', size=12,
         color='muted')
    toc = ['Run summary', 'Main steady states: dissociation rate and poisoning', 'What covers the surface',
           'Did the solver reach a true steady state?', 'Poisoning curves: how much impurity is tolerated?',
           'Temperature sweep', 'Hydrogen rates in Arrhenius form', 'Pressure sweep', 'Equilibrium coverage of hydrogen',
           'Reaction orders', 'Rate control: which species and steps set the rate', 'Reaction paths of hydrogen',
           'Poisoning over time', 'Reaction-path diagrams at the other points', 'Output files',
           'The settings of this run']
    nums = [str(i + 1) for i in range(13)] + ['A', 'B', 'C']
    text(sl, 0.7, 2.42, 0.42, 3.4, nums, size=12, color='muted', align='r')
    text(sl, 1.2, 2.42, 4.85, 3.4, toc, size=12)
    card(sl, 6.35, 1.66, 6.43, 2.62, 'Inside every section')
    parts = (('Goal', 'the question the analysis answers'), ('Key numbers', 'the results to quote, in words'),
             ('Tables', 'split by series; long ones as a summary'), ('Figures', 'each with a caption'),
             ('How to read the tables', 'what each column means'),
             ('How to read the figures', 'each panel, and what each trend means'))
    for i, (h, t) in enumerate(parts):
        y = 2.1 + i * 0.35
        shape(sl, 'rrect', 6.55, y, 2.6, 0.3, fill='white', line='faint', radius=0.06, content=h, size=12, bold=True)
        text(sl, 9.3, y - 0.01, 3.4, 0.32, t, size=13, anchor='m')
    card(sl, 6.35, 4.42, 6.43, 1.46, 'Settings: the report section')
    text(sl, 6.55, 4.82, 6.1, 1.0, ['`max_rows` 60: a longer table becomes a summary and its file name',
                                    '`path_diagrams` reference+appendix: the rest in Appendix A',
                                    'LaTeX through latexmk; without it, only the .tex is written'], size=13)
    actions(sl, ['Start with the run summary: system, conditions, warnings.',
                 'In each section: the goal, then the key numbers.',
                 'Read each table and figure with its guide.',
                 'Rebuild it any time with report: no solving.'])
    notes(sl, ['report.py: write_tex (lines 1621-1639: title page, contents, run summary, one section per analysis with',
               'outputs, in ANALYSIS_ORDER, and the appendices), Report.section (442-475: Goal, the parts, "How to read the',
               'tables" from column_list 438-440, "How to read the figures" from how_to_read 424-436: per panel what it',
               'shows, how to read it, and trend => meaning), title_page (574-603: run, settings, output folder, mechanism',
               'and its SHA-256, the catalyst from input.py, build time, package commit, Cantera version), run_summary',
               '(604-726: the system, conditions, analyses, settings that differ from the defaults, mechanism warnings,',
               'convergence), build (1641-1670: latexmk -pdf, 900 s limit; ReportError with the LaTeX errors).',
               'The texts of every section: report_text.py SECTIONS (title, goal, columns, figure specs).',
               'Settings: config.py REPORT_DEFAULTS (lines 42-48): enabled True, max_rows 60, path_diagrams',
               "'reference+appendix', latexmk None = on PATH or /Library/TeX/texbin.",
               f'{RUN}/output/report/try_2_Cu_report.pdf: 133 pages, "H2 on Cu111: surface kinetics report", built',
               '2026-10-15 23:49, after the last report.py commit (7d4ac3f, 22:09). Its sections are listed from the .tex',
               '(no diagnostics section: switch_off = []). Docs: "The run report" (line 1333).'])


SLIDES = [s59_commands, s60_output, s61_report]
