"""Section D (slides 19-23): the run folder and settings.py, the checks at load (Stage 0), the
roles (Stage 1), find_role_reaction(), and the check command. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PD = 'D · Settings and loading · Stages 0–1'
RUN = 'runs/try_2/Cu'


def s19_settings(d):
    sl = d.slide('A run is one folder: RMG’s files plus one settings.py the package reads',
                 'RMG’s files stay as RMG wrote them. settings.py says what the package should do with them.',
                 part=PD, tag='settings.py')
    text(sl, 0.55, 1.62, 5.9, 0.34, f'`{RUN}/`', size=15, bold=True)
    files = [('input.py', "RMG's input: catalyst, molecules, 15 reactors", 'panel2'),
             ('chem_annotated.yaml', 'the mechanism Cantera loads', 'panel2'),
             ('seed/   species/', "RMG's copy of the core; structure pictures", 'panel2'),
             ('settings.py', "the package's: roles, series, analyses", 'white'),
             ('output/', 'tables, figures, the PDF report (written by the package)', 'white')]
    for i, (fn, what, fill) in enumerate(files):
        y = 2.02 + i * 0.66
        shape(sl, 'doc', 0.75, y, 2.45, 0.56, fill=fill, line='faint', content=f'`{fn}`', size=13, align='l',
              anchor='m', margin=0.1)
        text(sl, 3.35, y + 0.08, 3.2, 0.5, what, size=13, color='muted' if fill == 'panel2' else 'ink')
    text(sl, 0.55, 5.4, 5.9, 0.5, 'Anything settings.py leaves out takes its default (docs, Appendix B).', size=13,
         color='muted')
    lines = [(9, "mechanism = 'chem_annotated.yaml'"), (11, "expected_site_density = (2.943e-09, 'mol/cm2')"),
             (14, "probe_gas = 'H2(1)'"), (15, "empty_site = 'vacantX(8)'"), (16, "probe_adsorbate = '[PtH](33)'"),
             (17, "molecular_adsorbate = '[H][H].[Pt](34)'"), (21, "impurities = {'H2O': 'H2O(2)', …}   # 7, by name"),
             (45, 'temperatures = [650.0, 800.0, 1000.0, 1200.0]'), (46, 'reference_pressure = 1.0e5'),
             (47, 'pressures = [1.0e3, 1.0e5, 1.0e7]'), (48, 'pressure_temperatures = [650.0, 1200.0]'),
             (51, "impurity_ratio = {'H2O': 10 / 190, …}"), (57, "analyses = ['steady_state', …]   # 13"),
             (66, "diagnostics = {…, 'switch_off': []}")]
    code(sl, 6.7, 1.66, 6.08, 3.38, [f'{n:>3}  {s}' for n, s in lines], size=12)
    ref(sl, 6.7, 5.08, 6.08, f'{RUN}/settings.py (excerpt)', align='l')
    caveat(sl, 6.7, 5.38, 6.08, 0.5, 'The loader refuses a file that looks like an RMG input (config.py:246–248).',
           size=13)
    actions(sl, ['Keep RMG’s files as RMG wrote them.',
                 'Write one settings.py: roles, series and analyses.',
                 'Leave everything else to the defaults.',
                 'Every result goes to output/.'])
    notes(sl, [f'{RUN}/settings.py sections: files (lines 8-11), species roles (13-30), species_labels (32-42),',
               'conditions (44-54: the temperature series at reference_pressure, the pressure series at',
               'pressure_temperatures, impurity_ratio per impurity from RMG reactor 1), what to run (56-59), analysis',
               'settings (61-80: diagnostics, flux). The excerpt shows real lines; "…" shortens two dictionaries.',
               'Defaults: config.py DEFAULTS (docs Appendix B). load_config (config.py:239) reads the file with runpy;',
               'lines 246-248 stop it if it contains database(, surfaceReactor(, catalystProperties( or simulator(.',
               'Docs: "Writing a settings.py" (simulation_workflow.md, line 1194).'])


def s20_stage0(d):
    sl = d.slide('Stage 0: one number catches a mis-converted file; broken duplicate flags are repaired',
                 'Both checks run every time the package loads a mechanism, before anything is solved.',
                 part=PD, tag='Stage 0 · model.py')
    # left: the site density
    card(sl, 0.55, 1.66, 5.95, 3.3, 'Site density: one number, two units')
    text(sl, 0.75, 2.12, 5.6, 0.3, "`expected_site_density = (2.943e-09, 'mol/cm2')`", size=12)
    text(sl, 0.75, 2.42, 5.6, 0.3, 'in settings.py, copied from RMG’s metal library', size=12, color='muted')
    step(sl, 0.75, 2.8, 2.6, 0.8, ['expected Γ (× 10)', ('**2.943 × 10^{−8} kmol/m²**', {'size': 15})], size=13)
    step(sl, 3.7, 2.8, 2.6, 0.8, ['the YAML’s Γ (mol, m)', ('**2.943 × 10^{−8} kmol/m²**', {'size': 15})], size=13)
    text(sl, 3.7, 3.62, 2.6, 0.3, '`site-density: 2.943e-05`', size=12)
    tick(sl, 3.52, 3.2, ok=True, d=0.3)
    text(sl, 0.75, 4.02, 5.6, 0.85, 'They must agree within 0.1 %. If not, the load stops: "wrong or mis-converted '
         'file?" (model.py:96–102).', size=13)
    # right: duplicate flags
    card(sl, 6.83, 1.66, 5.95, 3.3, 'Duplicate flags without a partner')
    flow = ['Cantera refuses the whole file: “No duplicate found for declared duplicate reaction”',
            'The package copies it without those flags: output/<stem>.duplicates_fixed.yaml',
            'and loads the copy. The flag changes no rate; a warning says how many.']
    for i, t in enumerate(flow):
        step(sl, 7.03, 2.12 + i * 0.78, 5.55, 0.62, t, size=13, align='l')
        if i < 2:
            arrow(sl, 9.8, 2.74 + i * 0.78, 9.8, 2.9 + i * 0.78)
    text(sl, 7.03, 4.45, 5.6, 0.4, 'try_2/Cu: none needed.  try_2/Fe: 17 flags.  try_2/Ni: 1.', size=13, bold=True)
    caveat(sl, 0.55, 5.08, W - 2 * LM, 0.8, 'The trap it catches: an export in kmol units that copies the mol/cm² number '
           'would read 2.943 × 10^{−9} kmol/m² here, 10× too small. The load stops before anything is solved.', size=13)
    actions(sl, ['Convert the expected site density to kmol/m² (× 10).',
                 'Compare it with the YAML’s own value, within 0.1 %.',
                 'If Cantera refuses only duplicate flags, load a copy without them.',
                 'Stop on anything else: the file is wrong.'])
    notes(sl, ['Site density: model.py lines 96-102 (np.isclose(surf.site_density, expected, rtol=1e-3));',
               'site_density_kmol_m2 in config.py converts (value, unit); 1 mol/cm2 = 10 kmol/m2.',
               f'{RUN}/settings.py line 11: (2.943e-09, "mol/cm2"), the Cu111 value of RMG-database metal.py;',
               'chem_annotated.yaml line 42: site-density 2.943e-05 in mol/m (units line 5). run.log: 2.9430e-08 kmol/m2.',
               'Duplicate flags: model.py _load_mechanism (137-163) and remove_unmatched_duplicates (34-66): the copy\'s',
               'first line records the source\'s sha256, so it is reused while the mechanism is unchanged.',
               'Counts from the first lines of runs/try_2/Fe and try_2/Ni output/chem_annotated.duplicates_fixed.yaml.',
               'Docs Stage 0, procedure step 2: "A writer that copies the mol cm-2 number into a kmol-based file gives a',
               'value 10x too small." RMG may write several exports (cantera1/, cantera2/): use the one that matches.',
               'Docs: Stage 0, procedure step 2 and the checks table.'])


def s21_roles(d):
    sl = d.slide('Stage 1: name the roles, so the code knows what to watch',
                 'Every role is given by its name in the mechanism. The model checks that each one exists.',
                 part=PD, tag='settings.py · model.py')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'Roles, given by name')
    rows = [('H2(1)', 'probe gas'), ('vacantX(8)', 'empty site'), ('[PtH](33)', 'dissociated adsorbate, H∗'),
            ('[H][H].[Pt](34)', 'molecular adsorbate, H₂∗ (optional)')]
    for i, (nm, role) in enumerate(rows):
        y = 2.15 + i * 0.9
        cx, cy = 1.05, y + 0.36
        if i == 0:
            h2(sl, cx, cy, d=0.3)
        elif i == 1:
            site(sl, cx, cy, d=0.38)
        elif i == 2:
            adsorbate(sl, cx, cy, 'H', 'H', d=0.38)
        else:
            site(sl, cx, cy + 0.12, d=0.3)
            h2(sl, cx, cy - 0.18, d=0.24)
        text(sl, 1.55, y, 2.8, 0.36, nm, size=15, font=CODE, bold=True)
        text(sl, 1.55, y + 0.36, 4.8, 0.36, role, size=14, color='muted')
    card(sl, 6.83, 1.66, 5.95, 4.22, 'Impurities, given by name')
    imps = [('H₂O', 'H2O', 'H2O(2)'), ('NH₃', 'NH3', 'NH3(9)'), ('CH₄', 'CH4', 'CH4(5)'), ('CO', 'CO', 'CO(4)'),
            ('CO₂', 'CO2', 'CO2(3)'), ('O₂', 'O2', 'O2(6)'), ('N₂', 'N2', 'N2(7)')]
    for i, (lab, key, name) in enumerate(imps):
        y = 2.12 + i * 0.42
        chip(sl, 7.0, y, lab, key=key, w=0.85, h=0.34, size=13)
        arrow(sl, 7.95, y + 0.17, 8.4, y + 0.17)
        text(sl, 8.5, y, 1.8, 0.36, name, size=14, font=CODE, bold=True)
        text(sl, 10.25, y, 1.7, 0.36, 'a gas species', size=12, color='muted')
        tick(sl, 12.4, y + 0.17, d=0.28)
    text(sl, 7.0, 5.1, 5.6, 0.75, "Or by composition, e.g. {'C': 1, 'O': 1}: exactly one gas species must match, "
         'or the load stops: “give the species name instead”.', size=13, color='muted')
    actions(sl, ['Name the probe gas, the empty site and the adsorbates.',
                 'Check that each name exists in the mechanism.',
                 'Resolve each impurity: a name, or a composition with one match.',
                 'Stop with a clear message if anything is missing.'])
    notes(sl, [f'{RUN}/settings.py lines 14-29 (roles and the 7 impurities, by name); reference_label "H2 only".',
               'model.py lines 104-112 check the role names; line 128 resolves the impurities with _resolve_gas',
               '(lines 169-179: a name must be a gas species; a composition must match exactly one, else ModelError).',
               "run.log: impurities: {'H2O': 'H2O(2)', 'NH3': 'NH3(9)', 'CH4': 'CH4(5)', 'CO': 'CO(4)', 'CO2': 'CO2(3)',",
               "'O2': 'O2(6)', 'N2': 'N2(7)'}. Names are exact; compositions survive renumbering between RMG runs,",
               'but isomers sharing a composition must be named (doc Stage 1).'])


def s22_find_reaction(d):
    sl = d.slide('find_role_reaction(): find R0 by what it does, not where it sits',
                 'Reaction numbers change between exports and RMG runs. Species and stoichiometry do not.',
                 part=PD, tag='package · model.py:69–82')
    text(sl, 0.55, 1.62, 7.3, 0.34, 'looking for:  `{H2(1), vacantX(8)}`  ⇌  `{[PtH](33)}`', size=15)
    rx = [('#0', 'H2(1) + 2 vacantX(8) <=> 2 [PtH](33)', True, 'hit: sign +1, n = 2 sites'),
          ('#1', 'H2(1) + vacantX(8) <=> [H][H].[Pt](34)', False, 'products differ'),
          ('#2', 'H2O(2) + vacantX(8) <=> O.[Pt](39)', False, 'no match'),
          ('#3', 'CO(4) + vacantX(8) <=> OCX(15)', False, 'no match'),
          ('⋮', '… 1,391 more, up to #1394', None, 'no match')]
    for i, (num, eq, ok, why) in enumerate(rx):
        y = 2.05 + i * 0.5
        if i % 2 == 0:
            shape(sl, 'rect', 0.55, y - 0.03, 7.3, 0.47, fill='panel', line=None)
        text(sl, 0.62, y, 0.6, 0.4, num, size=14, font=CODE, bold=True)
        text(sl, 1.2, y, 4.65, 0.4, eq, size=13, font=CODE)
        if ok is not None:
            tick(sl, 5.98, y + 0.2, ok=ok, d=0.28)
        text(sl, 6.22, y, 1.63, 0.45, why, size=12, color='ok' if ok else 'muted', bold=bool(ok))
    ref(sl, 0.55, 4.55, 7.3, 'loop over every surface reaction: lines 74–78, either direction')
    arrow(sl, 7.95, 3.2, 8.4, 3.2, color='ink', lw=1.75)
    step(sl, 8.45, 2.75, 1.65, 0.9, ['hits', ('`[(0, +1, 2)]`', {'size': 13})])
    arrow(sl, 10.1, 3.2, 10.3, 3.2)
    gate(sl, 10.3, 2.55, 2.45, 1.3, 'exactly one hit? (line 79)', size=13)
    arrow(sl, 11.52, 3.85, 11.52, 4.22)
    yes_no(sl, 11.6, 3.88, 'yes')
    step(sl, 10.3, 4.22, 2.45, 0.6, 'return (0, +1, 2)', fill='okbg', line='ok', size=14)
    text(sl, 8.3, 4.3, 1.9, 0.6, 'R0 = reaction 0, forward, 2 empty sites', size=12, color='muted')
    text(sl, 10.35, 1.8, 2.45, 0.6, 'no: ModelError (lines 80–81)', size=12, color='ink', bold=True, align='c')
    arrow(sl, 11.52, 2.55, 11.52, 2.25, color='dark')
    cases = ['Written the other way round: sign −1. probe_rates() then swaps forward and reverse (lines 216–217).',
             'Two matches, e.g. a duplicate reaction: ModelError. RMG can make the same reaction twice.',
             'Other reactions using H₂ directly: a warning. None here; try_3–5/Cu each have one (CX + H₂ ⇌ CH₂X).']
    for i, t in enumerate(cases):
        x = 0.55 + i * 4.13
        shape(sl, 'rrect', x, 4.95, 3.97, 0.95, fill='white', line='rule', radius=0.08)
        text(sl, x + 0.1, 4.98, 3.77, 0.9, t, size=13, anchor='m')
    actions(sl, ['Compare each reaction’s species sets with the wanted ones.',
                 'Accept either direction, and remember which (the sign).',
                 'Read the number of empty sites from the stoichiometry.',
                 'Insist on exactly one match.'])
    notes(sl, ['Hydrogen-Surface-kinetics model.py find_role_reaction, lines 69-82. Called in Model.__init__ lines',
               '115-120 for the probe reaction and again for the molecular adsorption (finds #1).',
               'Other reactions using the probe gas: lines 121-126; in try_2/Cu only #0 and #1 contain H2(1).',
               'try_3, try_4 and try_5 Cu run.log: "1 other surface reactions use the probe gas directly ... CX(296) + H2(1)',
               '<=> CH2X(297)". Equations are written as Cantera prints them; the YAML writes R0 as',
               '"vacantX(8) + vacantX(8) + H2(1) <=> [PtH](33) + [PtH](33)": comparing species SETS makes order and',
               'repetition irrelevant, and n_sites = r.reactants[site] = 2 (line 76).',
               'The last surface reaction, #1394 (2 vacantX(8) + [Pt]OC=N(1282) <=> [PtH](33) + [Pt]OC(=N)[Pt](1244)), also',
               'makes H*, but its reactants differ: no match.'])


def s23_check(d):
    sl = d.slide('check: see the roles, both series and every solve before running anything',
                 '`bash run.sh runs/try_2/Cu check` loads the settings and the mechanism, and solves nothing.',
                 part=PD, tag='package · cli.py, analyses/plan.py')
    out = ['$ bash run.sh runs/try_2/Cu check',
           '  probe reaction #0: H2(1) + 2 vacantX(8) <=> 2 [PtH](33)  (+1, 2 sites)',
           '  molecular adsorption #1: H2(1) + vacantX(8) <=> [H][H].[Pt](34)',
           '  6 sticking coefficients exceed 1 at some temperature',
           'temperature series: 650, 800, 1000, 1200 K at 1 bar',
           'pressure series:    0.01, 1, 100 bar at 650 K',
           'pressure series:    0.01, 1, 100 bar at 1200 K',
           '72 main conditions (8 points x cases)',
           'solves per analysis (nothing is solved by check):',
           '      analysis  points  solves',
           '  steady_state       8      72',
           '        curves       8     680',
           '         sweep      12     108',
           'pressure_sweep      18     162',
           '        orders       8     304',
           '  rate_control       8    3392',
           '    transients       8     224']
    code(sl, 0.55, 1.62, 7.75, 3.7, out, size=12)
    line_y = [1.62 + 0.09 + 0.2 * i + 0.1 for i in range(len(out))]
    ann = [(1, 1.86, 'roles and reactions found (Stage 1)', 'muted'),       # (code line, label top, text)
           (5, 2.62, 'both series (Stage 3)', 'muted'),
           (7, 3.08, 'every case at every point', 'muted'),
           (9, 3.54, 'solves needed, and how many are cached', 'muted'),
           (15, 4.62, 'rate control dominates: narrow it first', 'amber')]
    for i, yt, t, col in ann:
        arrow(sl, 8.95, yt + 0.15, 8.38, line_y[i], color='faint', lw=1.0)
        text(sl, 9.0, yt, 3.8, 0.32, t, size=13, color=col)
    actions(sl, ['Load the settings and the mechanism.',
                 'Show the roles, the series and every main condition.',
                 'Count the solves each analysis needs, and how many are cached.',
                 'Narrow an expensive analysis before running it.'])
    notes(sl, ['cli.py cmd_check (from line 26) prints Model.summary, the series and the main conditions table; analyses/',
               'plan.py counts (from line 37) gives the solves per analysis, matching cache files by the same fingerprints as',
               'the run (Stage 6), so "cached" and "to_solve" are exact for the current settings.',
               f'Output excerpt from {RUN}/run.log, check of 2026-10-15 22:32 (probe line shortened; the conditions table',
               'and the columns cached / to_solve / note are left out; total to solve was 0 then: every solve cached).',
               'rate_control: about 53 solves per case and point (one start each). Narrow an analysis with its own',
               "series / temperatures / pressures / cases, e.g. rate_control = {'cases': ['H2 only', 'CO'], 'series': 'temperature'}.",
               'Docs: "Running an analysis" (line 1159).'])


SLIDES = [s19_settings, s20_stage0, s21_roles, s22_find_reaction, s23_check]
