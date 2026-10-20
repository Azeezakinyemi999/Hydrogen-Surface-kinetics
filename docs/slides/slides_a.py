"""Section A (slides 1-5): the question, the four questions, the pipeline, and the two ways to
simulate. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PA = 'A · The question'
RUN = 'runs/try_2/Cu'
COMMIT = '51f05c2'


def s01_title(d):
    sl = d.prs.slides.add_slide(d.layout)
    d.count += 1
    bg = sl.background.fill
    bg.solid()
    bg.fore_color.rgb = rgb('white')
    text(sl, 0.9, 1.15, 11.6, 0.35, 'HYDROGEN-SURFACE-KINETICS · HOW THE WORKFLOW WORKS', size=12, color='faint',
         bold=True)
    text(sl, 0.9, 1.5, 11.8, 0.95, 'From an RMG mechanism to H₂ kinetics on Cu(111)', size=34, font=TITLE)
    text(sl, 0.9, 2.45, 11.6, 0.5, 'What each step of the workflow does, and why. Running example: try_2/Cu',
         size=22, color='muted')
    xs = slab(sl, 1.05, 5.0, 18, occ={3: ('H', 'H', 1), 4: ('H', 'H', 1), 10: ('CO', 'CO', 1),
                                      14: ('H', 'H', 1), 16: ('H2O', 'H₂O', 1)})
    cx = (xs[7] + xs[8]) / 2
    h2(sl, cx, 3.55)
    arrow(sl, cx, 3.8, cx, 4.62, color='H', lw=2)
    h2(sl, xs[1], 3.85)
    gas(sl, xs[12], 3.6, 'CO')
    gas(sl, xs[16] + 0.15, 3.75, 'H2O', 'H₂O')
    text(sl, 0.9, 6.25, 11.6, 0.35, 'Azeez Akinyemi · October 2026', size=16, color='muted')
    text(sl, 0.9, 6.62, 11.6, 0.35,
         f'Built from docs/simulation_workflow.md, {RUN} (RMG input.py, chem_annotated.yaml, settings.py) and the '
         f'package Hydrogen-Surface-kinetics (commit {COMMIT})', size=12, color='faint')
    notes(sl, ['Audience: me, the researcher. Goal: understand what each step does and why; results are left out.',
               f'Running example: {RUN}, the RMG run try_2/Cu (Danielle). Its folder holds RMG\'s input.py and',
               'chem_annotated.yaml, and the package\'s settings.py.',
               'Picture: schematic side view of Cu(111) sites (open circles = empty sites).',
               f'Sources: docs/simulation_workflow.md; {RUN}/input.py, chem_annotated.yaml and settings.py;',
               f'Hydrogen-Surface-kinetics at commit {COMMIT}.'])


def s02_question(d):
    sl = d.slide('The question: do gas impurities slow H₂ splitting on copper?',
                 'H₂ needs two empty neighbouring sites to split. Anything sitting on the surface can take them.',
                 part=PA, tag='concept')
    # left: H2 only
    card(sl, 0.55, 1.66, 5.95, 3.72, 'H₂ only')
    xs = slab(sl, 1.05, 4.2, 8, occ={5: ('H', 'H', 1), 6: ('H', 'H', 1)})
    cx = (xs[2] + xs[3]) / 2
    h2(sl, cx, 2.62)
    arrow(sl, cx, 2.85, cx, 3.86, color='H', lw=2)
    label(sl, cx + 0.25, 2.95, 2.9, 'lands on two empty neighbours and splits')
    label(sl, xs[4] - 0.2, 4.86, 2.6, 'H∗ + H∗ from an earlier split', size=14, color='muted')
    # right: with an impurity
    card(sl, 6.83, 1.66, 5.95, 3.72, 'H₂ with an impurity, here CO')
    xs2 = slab(sl, 7.33, 4.2, 8, occ={0: ('CO', 'CO', 1), 2: ('CO', 'CO', 1), 4: ('CO', 'CO', 1),
                                      6: ('CO', 'CO', 1)})
    cx2 = (xs2[3] + xs2[4]) / 2
    h2(sl, cx2, 2.62)
    arrow(sl, cx2, 2.85, cx2, 3.42, color='H', lw=2, end=False)
    line(sl, cx2 - 0.22, 3.45, cx2 + 0.22, 3.45, color='dark', lw=3)
    label(sl, cx2 + 0.3, 3.12, 2.4, 'no two empty neighbours')
    gas(sl, xs2[1], 2.45, 'CO')
    gas(sl, xs2[6] + 0.1, 2.75, 'CO')
    label(sl, xs2[0] - 0.2, 4.86, 5.0, 'CO∗ holds every other site', size=14, color='muted')
    text(sl, 0.55, 5.45, W - 2 * LM, 0.42,
         'R0 in the mechanism:  H₂ + 2 ∗ ⇌ 2 H∗      written there as  `H2(1) + 2 vacantX(8) <=> 2 [PtH](33)`',
         size=16, align='c', anchor='m')
    actions(sl, ['An H₂ molecule from the gas lands on the surface.',
                 'To split, it needs two empty neighbouring sites.',
                 'An impurity such as CO sits on sites.',
                 'Fewer empty pairs means fewer splitting events.'])
    notes(sl, [f'R0 = surface reaction #0 of {RUN}/chem_annotated.yaml (line 3494):',
               'vacantX(8) + vacantX(8) + H2(1) <=> [PtH](33) + [PtH](33).',
               'Rate: sticking-Arrhenius, A = 0.03, b = 0, Ea = 5000 J/mol (5 kJ/mol). Note: "Template reaction:',
               'Surface_Adsorption_Dissociative", i.e. estimated by RMG, not taken from a library.',
               'The package finds it at load (run.log): "probe reaction #0: H2(1) + 2 vacantX(8) <=> 2 [PtH](33)',
               '(direction +1, 2.0 empty sites)".',
               '∗ = empty site, vacantX(8); H∗ = [PtH](33); CO∗ = OCX(15), "CO* (1 site)" in settings.py species_labels.',
               'One row of sites stands for the two-dimensional surface (schematic).',
               'The question: docs/simulation_workflow.md section 1.1.'])


def s03_four_questions(d):
    sl = d.slide('The workflow answers four questions about the surface',
                 'Each later stage exists to answer one of them.', part=PA, tag='concept')
    xs0 = [0.55 + i * (2.87 + 0.25) for i in range(4)]
    heads = ['How fast does H₂ split?', 'What covers the surface?', 'What do impurities do?',
             'Can it be trusted?']
    body = ['The forward, reverse and net turnover frequency (TOF) of R0, per site, and the same per area.',
            'The steady-state coverage θ of every surface species, and θ at equilibrium with the gas.',
            'The same numbers with each impurity, and whether lost empty sites explain the change.',
            'Where each number in the mechanism came from, and whether conclusions survive a what-if.']
    foot = ['Stages 4–7', 'Stage 8 and 11.2c', 'Stages 7 and 11', 'Stage 10']
    for i, x in enumerate(xs0):
        card(sl, x, 1.66, 2.87, 5.0, heads[i], size=16)
        text(sl, x + 0.15, 4.36, 2.57, 1.5, body[i], size=14)
        text(sl, x + 0.15, 6.18, 2.57, 0.3, foot[i], size=14, color='muted', bold=True)
    # 1: one splitting event on a pair of sites, counted two ways
    x = xs0[0]
    xs = slab(sl, x + 0.44, 3.08, 4, pitch=0.5, d=0.36, metal_label='')
    h2(sl, (xs[1] + xs[2]) / 2, 2.3, d=0.3)
    arrow(sl, (xs[1] + xs[2]) / 2, 2.5, (xs[1] + xs[2]) / 2, 2.83, color='H', lw=2)
    text(sl, x + 0.15, 3.68, 2.57, 0.56, ['per site: events per second', 'per area: mol H₂ m⁻² s⁻¹'], size=12,
         color='muted', align='c')
    # 2: a covered surface
    x = xs0[1]
    grid(sl, x + 0.42, 2.3, 5, 3, pitch=0.42, d=0.32,
         occ={(0, 0): ('H', 'H', 1), (0, 2): ('H', 'H', 1), (0, 3): ('H', 'H', 1), (1, 1): ('H', 'H', 1),
              (1, 4): ('H', 'H', 1), (2, 0): ('H', 'H', 1), (2, 2): ('H', 'H', 1)})
    text(sl, x + 0.15, 3.68, 2.57, 0.3, 'θ = share of sites', size=12, color='muted', align='c')
    # 3: reference vs impurity
    x = xs0[2]
    grid(sl, x + 0.18, 2.4, 3, 3, pitch=0.36, d=0.3, plate=False,
         occ={(0, 0): ('H', 'H', 1), (1, 1): ('H', 'H', 1), (2, 2): ('H', 'H', 1)})
    arrow(sl, x + 1.41, 2.862, x + 1.75, 2.862, color='muted', lw=1.5)     # between the offset middle rows
    grid(sl, x + 1.62, 2.4, 3, 3, pitch=0.36, d=0.3, plate=False,
         occ={(0, 0): ('CO', 'CO', 1), (0, 1): ('CO', 'CO', 1), (1, 0): ('CO', 'CO', 1), (1, 2): ('CO', 'CO', 1),
              (2, 1): ('CO', 'CO', 1), (2, 2): ('H', 'H', 1)})
    text(sl, x + 0.1, 3.68, 1.3, 0.3, 'H₂ only', size=12, color='muted', align='c')
    text(sl, x + 1.5, 3.68, 1.3, 0.3, 'with CO', size=12, color='muted', align='c')
    # 4: two real provenance notes of this mechanism
    x = xs0[3]
    shape(sl, 'rrect', x + 0.2, 2.2, 2.47, 0.74, fill='white', line='ok', lw=1.25, radius=0.06,
          content=['Library reaction:', 'O2_dissociative_', 'adsorption'], size=12, font=CODE)
    shape(sl, 'rrect', x + 0.2, 3.04, 2.47, 0.74, fill='white', line='amber', lw=1.25, radius=0.06,
          content=['Template reaction:', 'Surface_Adsorption_', 'Dissociative'], size=12, font=CODE)
    notes(sl, ['The four questions of docs/simulation_workflow.md section 1.1.',
               'TOF: section 2.5 and Stage 7. Per area (Format 2): Stage 7, per_area.py, docs/rates_units_equilibrium.md.',
               'Coverages: Stage 8. Equilibrium coverage: Stage 11.2c (analyses/equilibrium.py).',
               'Impurities and site blocking: Stages 7 and 11.',
               'Trust: Stage 10 (binding check, switch-off what-ifs, sources and sinks, provenance).',
               f'Card 4 shows two real note types of {RUN}/chem_annotated.yaml: "Library reaction:',
               'O2_dissociative_adsorption" (2 reactions, e.g. line 3531) and R0\'s "Template reaction:',
               'Surface_Adsorption_Dissociative". The other pictures are schematic.'])


def s04_pipeline(d):
    sl = d.slide('The pipeline: RMG writes the mechanism, Cantera runs it, a report collects it',
                 'Two branches share the first stages and meet again in the diagnostics.', part=PA, tag='concept')
    bh = 0.8
    yA, yB, yC = 1.7, 3.12, 4.55
    # row A: inputs, Stage 0, Stage 1, settings.py
    boxes_a = [(0.55, 2.25, ['**RMG job**', '`input.py`: 15 reactors'], 'panel2'),
               (3.2, 2.75, ['**Exports**', '`chem_annotated.yaml`, seed, species'], 'panel2'),
               (6.35, 2.15, ['**Stage 0**', 'check the files belong together'], 'white'),
               (8.9, 2.15, ['**Stage 1**', 'load in Cantera, name the roles'], 'white')]
    for x, w, t, f in boxes_a:
        step(sl, x, yA, w, bh, t, fill=f, size=13)
    for (x1, w1, _, _), (x2, _, _, _) in zip(boxes_a[:-1], boxes_a[1:]):
        arrow(sl, x1 + w1, yA + bh / 2, x2, yA + bh / 2)
    shape(sl, 'doc', 11.4, yA, 1.38, bh + 0.05, fill='white', line='faint', content=['`settings.py`', 'roles, series,', 'analyses'],
          size=12)
    arrow(sl, 11.4, yA + bh / 2, 11.05, yA + bh / 2)
    # row B: branch B, then diagnostics, further analyses and the report
    bw, gap = 1.79, 0.3
    xb = [0.55 + i * (bw + gap) for i in range(6)]
    shape(sl, 'rrect', 0.47, yB - 0.12, 3 * bw + 2 * gap + 0.16, 1.27, fill='panel', line=None, radius=0.08)
    text(sl, 0.55, yB + 0.86, 3 * bw + 2 * gap, 0.28, 'Branch B · steady state at fixed gas', size=12, color='muted')
    texts_b = [['**Stage 3**', 'two series: T and p'], ['**Stages 4–6**', 'solve, two starts, parallel, cached'],
               ['**Stages 7–9**', 'rates per site and per area, coverages'], ['**Stage 10**', 'can it be trusted?'],
               ['**Stage 11**', 'further analyses'], ['**Report**', 'one PDF: every table and figure']]
    for i, (x, t) in enumerate(zip(xb, texts_b)):
        step(sl, x, yB, bw, bh, t, fill='panel2' if i == 5 else 'white', size=13)
    for x1, x2 in zip(xb[:-1], xb[1:]):
        arrow(sl, x1 + bw, yB + bh / 2, x2, yB + bh / 2)
    # Stage 1 -> Stage 3 and Stage 2
    s1 = 8.9 + 2.15 / 2
    poly(sl, [(s1, yA + bh), (s1, 2.75), (0.33, 2.75), (0.33, yC + bh / 2), (0.62, yC + bh / 2)], end=True)
    arrow(sl, xb[0] + bw / 2, 2.75, xb[0] + bw / 2, yB)
    # row C: branch A
    step(sl, 0.62, yC, 4.08, bh, ['**Branch A · Stage 2** (optional)', "replay RMG's 15 reactors, draw flux diagrams"],
         size=13)
    s10 = xb[3] + bw / 2
    poly(sl, [(4.7, yC + bh / 2), (s10, yC + bh / 2), (s10, yB + bh)], end=True)
    # Stage 10 -> RMG: regenerate
    poly(sl, [(xb[3] + bw - 0.3, yB + bh), (xb[3] + bw - 0.3, 5.58), (0.2, 5.58), (0.2, yA + bh / 2), (0.55, yA + bh / 2)],
         end=True, dash=True, color='amber')
    text(sl, 2.0, 5.62, 5.0, 0.3, 'if not trustworthy: regenerate the mechanism', size=14, color='amber', italic=True)
    actions(sl, ['RMG builds the mechanism from its input.py.',
                 'Stages 0–1 check the files and load them, with the roles from settings.py.',
                 "Branch B solves the surface at fixed gas; branch A replays RMG's reactors.",
                 'Stage 10 judges, Stage 11 extends, and the report collects every result.'])
    notes(sl, ['Diagram 1 and the stage table of docs/simulation_workflow.md section 1.3.',
               f'15 = surfaceReactor blocks in {RUN}/input.py: 5 gases (the full mixture, H2:CO 8:2, H2:CO2 8:2,',
               'H2:CH4 9:1, H2:NH3 9:1) at 3 pressures (1, 100, 0.01 bar), each over 650-1200 K.',
               'settings.py is the package\'s file (roles, conditions, which analyses); input.py is RMG\'s.',
               'Stage 3 defines a temperature series and a pressure series (conditions.py). Stages 4-6: solver.py,',
               'runner.py. Stages 7-9: analyses/steady_state.py, composition.py; per area: per_area.py.',
               'Branch A: analyses/flux.py, optional, run with  bash run.sh runs/try_2/Cu --only flux  (needs RMG-Py).',
               'Report: report.py and report_text.py; output/report/<run>_report.pdf, built at the end of a run of every',
               'analysis in settings.py, or rebuilt with  bash run.sh runs/try_2/Cu report.'])


def _frame(sl, x, y, w, h, n_h2, n_co, label_txt):
    shape(sl, 'rrect', x, y, w, h, fill='white', line='faint', lw=1.25, radius=0.08)
    shape(sl, 'rect', x + 0.12, y + h - 0.3, w - 0.24, 0.18, fill='ni', line=None)
    spots = [(0.42, 0.35), (1.2, 0.5), (0.42, 0.8), (1.2, 0.95), (0.42, 1.25), (1.2, 1.4),
             (0.42, 1.7), (1.2, 1.85)]
    k = 0
    for _ in range(n_h2):
        gas(sl, x + spots[k][0], y + spots[k][1], 'H', 'H₂', w=0.52, h=0.3, size=12)
        k += 1
    for _ in range(n_co):
        gas(sl, x + spots[k][0], y + spots[k][1], 'CO', 'CO', w=0.52, h=0.3, size=12)
        k += 1
    text(sl, x, y + h + 0.04, w, 0.3, label_txt, size=14, color='muted', align='c')


def s05_two_ways(d):
    sl = d.slide('Two ways to simulate: a closed box, or a gas that never changes',
                 'In a closed box the surface uses up the gas, so "the rate at 20 % CO" has no single value.',
                 part=PA, tag='concept')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'Closed batch reactor · Stage 2')
    for i, (nh, lab) in enumerate([(6, 't = 0'), (4, 'later'), (2, 'much later')]):
        _frame(sl, 0.85 + i * 1.9, 2.25, 1.62, 2.45, nh, 2, lab)
        if i < 2:
            arrow(sl, 0.85 + i * 1.9 + 1.66, 3.45, 0.85 + (i + 1) * 1.9 - 0.04, 3.45)
    text(sl, 0.75, 5.18, 5.6, 0.6, 'The gas and the surface both change: H₂ is used up.', size=14, color='ink')
    card(sl, 6.83, 1.66, 5.95, 4.22, 'Fixed gas · Stages 3–11')
    shape(sl, 'rrect', 7.1, 2.25, 2.85, 2.45, fill='white', line='faint', lw=1.25, radius=0.08)
    pos = [(0.45, 0.4), (1.25, 0.32), (2.05, 0.45), (0.7, 0.95), (1.55, 0.88), (2.35, 1.05), (0.4, 1.5),
           (1.2, 1.45), (2.05, 1.6), (0.75, 2.0), (1.6, 2.02)]
    keys = ['H', 'H', 'CO', 'H', 'H', 'H', 'CO', 'H', 'H', 'H', 'CO']
    for (px, py), k in zip(pos, keys):
        gas(sl, 7.1 + px, 2.25 + py, k, 'H₂' if k == 'H' else 'CO', w=0.52, h=0.3, size=12)
    shape(sl, 'rect', 8.33, 4.52, 0.4, 0.1, fill='dark', line=None)
    text(sl, 7.1, 4.75, 2.85, 0.3, 'gas held at T, P, x', size=14, color='muted', align='c')
    arrow(sl, 10.05, 3.45, 10.45, 3.45)
    grid(sl, 10.65, 2.75, 4, 3, pitch=0.44, d=0.34,
         occ={(0, 1): ('H', 'H', 1), (0, 3): ('H', 'H', 1), (1, 0): ('CO', 'CO', 1), (1, 2): ('H', 'H', 1),
              (2, 1): ('H', 'H', 1), (2, 3): ('CO', 'CO', 1)})
    text(sl, 10.3, 4.05, 2.4, 0.62, 'only the surface evolves, to its steady state', size=14, color='muted', align='c')
    text(sl, 7.03, 5.18, 5.6, 0.6, 'Steady-state rates and coverages at exactly that gas.', size=14, color='ink')
    actions(sl, ['In a closed box, adsorption uses up H₂: the gas drifts.',
                 'So "the rate at this impurity level" has no single value.',
                 'Holding T, P and the gas composition fixed leaves only the surface to change.',
                 'Its steady state gives the rate at exactly that gas.'])
    notes(sl, ['docs/simulation_workflow.md section 1.2 and Figure 1.',
               'Closed batch (Stage 2): RMG\'s surfaceReactor, closed, isothermal, constant volume. V = RT/P for 1 mol',
               f'of gas, A = V x 1e5 1/m (surfaceVolumeRatio in {RUN}/input.py); analyses/flux.py simulate, lines 198-222,',
               'sampled from 1e-8 s to terminationTime (1000 s).',
               '"20 % CO": RMG reactor 2 of try_2/Cu is H2:CO = 8:2.',
               'Fixed gas (Stages 4-5): 1 m3 of gas over 1e-10 m2 of catalyst (reactor_area_per_volume, config.py:22),',
               'energy off, gas reactions off (solver.py fixed_gas_reactor, lines 98-107).',
               'Molecule counts in the frames are schematic.'])


SLIDES = [s01_title, s02_question, s03_four_questions, s04_pipeline, s05_two_ways]
