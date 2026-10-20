"""Section B (slides 6-10): how RMG builds the mechanism, for runs/try_2/Cu: its input, the
core/edge growth, where the numbers come from, what is exported, and how to read the names."""
from deck_lib import *  # noqa: F401,F403

PB = 'B · RMG builds the mechanism'
RUN = 'runs/try_2/Cu'


def s06_input(d):
    sl = d.slide('RMG starts from what input.py asks for',
                 'The catalyst, the starting molecules, the reactors to simulate, and when to stop.',
                 part=PB, tag='RMG input.py')
    xs0 = [0.55 + i * (2.87 + 0.25) for i in range(4)]
    heads = ['Catalyst: Cu(111)', 'Starting molecules', '15 reactors', 'Limits and stopping']
    for x, hd in zip(xs0, heads):
        card(sl, x, 1.66, 2.87, 4.24, hd)
    # catalyst: RMG looks the metal up
    x = xs0[0]
    text(sl, x + 0.15, 2.1, 2.6, 0.5, ['`catalystProperties(`', '`    metal="Cu111")`'], size=12)
    text(sl, x + 0.15, 2.66, 2.6, 0.3, 'site density', size=14, color='muted')
    text(sl, x + 0.15, 2.93, 2.6, 0.34, '2.943 × 10^{−9} mol/cm²', size=16, bold=True)
    text(sl, x + 0.15, 3.36, 2.6, 0.3, 'binding energies, eV', size=14, color='muted')
    for i, (el, v) in enumerate([('H', '−2.584'), ('C', '−4.960'), ('N', '−3.584'), ('O', '−4.208')]):
        yy = 3.66 + i * 0.36
        circle(sl, x + 0.35, yy + 0.15, 0.3, fill='H' if el == 'H' else 'dark', line=None, content=el, size=12,
               color='white', bold=True)
        text(sl, x + 0.62, yy, 1.6, 0.32, v, size=16)
    text(sl, x + 0.15, 5.12, 2.6, 0.7, "From RMG's metal library; used to shift Pt(111) values to Cu(111).",
         size=13, color='muted')
    # starting molecules
    x = xs0[1]
    text(sl, x + 0.15, 2.08, 2.6, 0.28, 'fed by the reactors', size=12, color='muted')
    fed = [('H₂', 'H'), ('H₂O', 'H2O'), ('CO₂', 'CO2'), ('CO', 'CO'), ('CH₄', 'CH4'), ('O₂', 'O2'), ('N₂', 'N2'),
           ('NH₃', 'NH3')]
    for i, (t, k) in enumerate(fed):
        chip(sl, x + 0.2 + (i % 2) * 1.25, 2.38 + (i // 2) * 0.42, t, key=k, w=1.1, h=0.34)
    text(sl, x + 0.15, 4.08, 2.6, 0.28, 'listed, but no reactor feeds them', size=12, color='muted')
    for i, t in enumerate(['N₂O', 'C₂H₆', 'C₂H₄', 'C₂H₂']):
        chip(sl, x + 0.2 + (i % 2) * 1.25, 4.38 + (i // 2) * 0.42, t, w=1.1, h=0.34)
    site(sl, x + 0.38, 5.42, d=0.32)
    text(sl, x + 0.62, 5.27, 2.2, 0.32, '`vacantX`: the empty site', size=13)
    # reactors: 5 gases x 3 pressures
    x = xs0[2]
    mixes = [('full mix', 'dark'), ('H₂ : CO = 8 : 2', 'CO'), ('H₂ : CO₂ = 8 : 2', 'CO2'), ('H₂ : CH₄ = 9 : 1', 'CH4'),
             ('H₂ : NH₃ = 9 : 1', 'NH3')]
    cols = [(0, '0.01'), (1, '1'), (2, '100')]               # pressure in bar; reactors 1-5 at 1, 6-10 at 100, 11-15 at 0.01
    first = {'0.01': 11, '1': 1, '100': 6}
    text(sl, x + 0.12, 2.08, 1.36, 0.28, 'p (bar):', size=12, color='muted')
    for j, p in cols:
        text(sl, x + 1.42 + j * 0.42, 2.08, 0.5, 0.28, p, size=12, color='muted', align='c')
    for i, (m, k) in enumerate(mixes):
        yy = 2.4 + i * 0.4
        text(sl, x + 0.12, yy, 1.36, 0.34, m, size=13, anchor='m')
        for j, p in cols:
            shape(sl, 'rect', x + 1.5 + j * 0.42, yy, 0.34, 0.34, fill=k, line=None, content=str(first[p] + i),
                  size=12, color='white', bold=True, margin=0.0)
    text(sl, x + 0.12, 4.47, 2.65, 0.6, 'Every reactor spans 650–1200 K, with 5 simulations each.', size=13)
    text(sl, x + 0.12, 5.06, 2.65, 0.8, 'full mix = reactor 1: H₂ 190, H₂O 10, CO 1, N₂ 0.85, CH₄ 0.8, O₂ 0.15, '
         'NH₃ 0.1, CO₂ 0.05', size=12, color='muted')
    # limits
    x = xs0[3]
    lim = ['≤ 2 C, O and N atoms, and ≤ 2 sites, per species', 'no surface singlet carbenes (forbidden)',
           'stop at 1000 s, or once 99 % of the H₂ is converted', 'reaction families: surface and default',
           'thermo: the Pt(111) surface library and 3 gas libraries', 'one reaction library: O₂ dissociative adsorption']
    for i, t in enumerate(lim):
        text(sl, x + 0.15, 2.1 + i * 0.62, 2.6, 0.6, t, size=14)
    actions(sl, ['Name the catalyst; RMG looks up its site density and binding energies.',
                 'List the molecules RMG starts from.',
                 'Define the 15 reactors RMG must simulate.',
                 'Limit species size and say when each simulation stops.'])
    notes(sl, [f'{RUN}/input.py: database() line 4 (thermoLibraries surfaceThermoPt111, primaryThermoLibrary,',
               "thermo_DFT_CCSDTF12_BAC, DFT_QCI_thermo and 'auto'; reactionLibraries O2_dissociative_adsorption;",
               "kineticsFamilies surface, default; kineticsEstimator 'rate rules'); catalystProperties line 13;",
               'forbidden line 17 (surface_singlet_carbene); 13 species() blocks, lines 27-103;',
               '15 surfaceReactor() blocks, lines 108-494: 1-5 at 1 bar, 6-10 at 100 bar, 11-15 at 0.01 bar; each',
               'temperature 650-1200 K, nSims = 5, surfaceVolumeRatio 1e5 1/m, stop at 99 % H2 or 1000 s;',
               'simulator line 497; model line 502; generatedSpeciesConstraints line 521.',
               'Site density and binding energies: RMG-database input/surface/libraries/metal.py, entry "Cu111"',
               '(lines 127-134; the same in all four local copies). The run\'s own database copy is not saved, but these',
               'values reproduce the Pt(111) -> Cu(111) shifts recorded in chem_annotated.yaml (slide 8).',
               "'auto' = libraries chosen by RMG (rmgpy/rmg/input.py:92, main.py:620 in the local RMG-Py)."])


def s07_growth(d):
    sl = d.slide('RMG grows the mechanism: react what it has, keep what carries flux',
                 'The core is the mechanism. The edge holds candidates that have not earned a place yet.',
                 part=PB, tag='RMG algorithm')
    # left: core and edge
    shape(sl, 'oval', 1.0, 1.75, 4.3, 3.45, fill=None, line='faint', lw=1.25, dash=True)
    shape(sl, 'oval', 1.95, 2.5, 2.4, 1.95, fill='panel2', line=None)
    text(sl, 1.95, 2.58, 2.4, 0.3, 'core', size=14, color='muted', align='c', bold=True)
    for (cx, cy, t, k) in [(2.55, 3.12, 'H₂', 'H'), (3.65, 3.12, 'CO', 'CO'), (2.45, 3.62, 'H₂O', 'H2O'),
                           (3.62, 3.62, 'H∗', 'H'), (3.1, 4.1, 'CO∗', 'CO')]:
        gas(sl, cx, cy, k, t, w=0.68, h=0.34)
    text(sl, 4.75, 1.62, 1.2, 0.3, 'edge', size=14, color='muted', bold=True)
    for (cx, cy) in [(1.55, 2.55), (1.3, 3.75), (2.4, 4.8), (4.8, 2.5), (4.25, 4.8)]:
        shape(sl, 'rrect', cx - 0.5, cy - 0.17, 1.0, 0.34, fill='white', line='faint', dash=True, radius='pill',
              content='new', size=13, color='muted', margin=0.0)
    shape(sl, 'rrect', 4.95 - 0.5, 3.45 - 0.17, 1.0, 0.34, fill='white', line='ink', lw=1.5, radius='pill',
          content='new', size=13, color='ink', bold=True, margin=0.0)
    arrow(sl, 4.42, 3.45, 4.12, 3.45, color='ink', lw=2.5)
    text(sl, 4.38, 3.02, 0.9, 0.3, 'moves in', size=13, color='ink')            # between core and edge oval
    shape(sl, 'rrect', 0.75, 5.32, 5.65, 0.52, fill='white', line='rule', radius=0.06,
          content='characteristic rate  R = √( Σ r_{k}^{2} )  over the core species k', size=14)
    # right: the loop
    s = [(6.9, 1.72, '1 · React the core species with each other, using the reaction families'),
         (6.9, 2.47, '2 · New products go to the edge, not to the core'),
         (6.9, 3.22, ['3 · Simulate each reactor with the core reactions only,', 'at 5 temperatures in 650–1200 K'])]
    for x, y, t in s:
        step(sl, x, y, 4.6, 0.58, t, align='l')
    for y in (2.3, 3.05):
        arrow(sl, 9.2, y, 9.2, y + 0.17)
    gate(sl, 7.35, 3.98, 3.7, 1.08, 'rate of an edge species > 0.5 R ?', size=14)
    arrow(sl, 9.2, 3.8, 9.2, 3.98)
    step(sl, 11.25, 4.22, 1.5, 0.6, 'move it into the core', size=13)
    arrow(sl, 11.05, 4.52, 11.25, 4.52)
    yes_no(sl, 10.9, 4.2, 'yes', w=0.4)
    poly(sl, [(12.0, 4.22), (12.0, 2.01), (11.5, 2.01)], end=True)
    step(sl, 6.95, 5.3, 4.5, 0.58, ['no: done. The core is the mechanism:', '111 species and 1,634 reactions.'],
         size=13, fill='okbg', line='ok')
    arrow(sl, 9.2, 5.06, 9.2, 5.3)
    yes_no(sl, 9.3, 5.02, 'no', w=0.4)
    actions(sl, ['Apply reaction families to the core species.',
                 'New products land in the edge, not in the core.',
                 'Simulate the 15 reactors and compare each edge rate with R.',
                 'Move species above 0.5 R into the core; repeat until none.'])
    notes(sl, [f'{RUN}/input.py model() (line 502): toleranceMoveToCore = 0.5, toleranceInterruptSimulation = 0.5,',
               'toleranceKeepInEdge = 0.0 (the branching options are commented out). simulator() line 497: atol 1e-18,',
               'rtol 1e-12. Families: surface, default (database() line 4).',
               'RMG-Py, local checkout d1f99b7 (2026-07-07; the run used RMG-Py commit b342f9b, YAML line 2):',
               'R = char_rate = sqrt(sum(core_species_rates^2)) (rmgpy/solver/base.pyx:812); edge ratio =',
               '|edge rate| / char_rate (824); moved to the core if above toleranceMoveToCore (1045-1047).',
               'nSims = 5 simulations per reactor and iteration, each at a condition from get_cond()',
               '(rmgpy/rmg/main.py:1001-1004), chosen far from earlier conditions (2513-2516).',
               'The exported core is chem_annotated.yaml: 40 gas + 71 surface species, 239 + 1,395 reactions.',
               'Chips labelled "new" stand for edge species (schematic). The run\'s RMG log is not in the run folder,',
               'so its run time is not known here.'])


def s08_estimates(d):
    sl = d.slide('RMG estimates almost every number: only 2 of 1,634 reactions come from a library',
                 'Rates come from reaction templates; adsorbate thermochemistry is borrowed from Pt(111) and shifted to Cu.',
                 part=PB, tag='chem_annotated.yaml')
    # rate constants
    card(sl, 0.55, 1.66, 5.95, 3.58, 'Rate constants of 1,634 reactions')
    text(sl, 0.75, 2.12, 1.0, 0.7, '2', size=40, font=TITLE, bold=True)
    text(sl, 1.55, 2.2, 4.8, 0.6, 'from the library O2_dissociative_adsorption: CO and O₂ adsorption', size=14)
    text(sl, 0.75, 2.88, 1.6, 0.7, '1,632', size=40, font=TITLE, bold=True)
    text(sl, 2.45, 2.96, 3.9, 0.6, 'estimated from reaction templates and rate rules: 239 gas, 1,393 surface', size=14)
    text(sl, 0.75, 3.68, 5.6, 0.28, 'R0, H₂ + 2 ∗ ⇌ 2 H∗:', size=13, color='muted')
    text(sl, 0.75, 3.94, 5.6, 0.3, '`Template reaction: Surface_Adsorption_Dissociative`', size=13)
    text(sl, 0.75, 4.34, 5.6, 0.28, 'R1, H₂ + ∗ ⇌ H₂∗ (physisorbed H₂):', size=13, color='muted')
    text(sl, 0.75, 4.6, 5.6, 0.3, '`Template reaction: Surface_Adsorption_vdW`', size=13)
    # thermochemistry
    card(sl, 6.83, 1.66, 5.95, 3.58, 'Thermochemistry of 71 surface species')
    x0, y0, wbar = 7.08, 2.15, 5.45
    xx = x0
    for n, col in [(61, 'faint'), (10, 'rule')]:
        w = wbar * n / 71
        shape(sl, 'rect', xx, y0, w, 0.46, fill=col, line='white', lw=1.0, content=str(n), size=14,
              color='white' if col == 'faint' else 'ink', bold=True)
        xx += w
    text(sl, x0, 2.7, 5.5, 0.3, '[[faint|■]]  61  from the Pt(111) surface library', size=14)
    text(sl, x0, 3.0, 5.5, 0.3, '[[rule|■]]  10  gas-phase thermo plus an adsorption estimate', size=14)
    text(sl, x0, 3.38, 5.5, 0.3, 'All shifted to Cu(111), in kJ/mol per bond to the surface:', size=14)
    for i, (el, v) in enumerate([('H', '+16'), ('C', '+199'), ('N', '+101'), ('O', '−38')]):
        cx = x0 + 0.18 + i * 1.38
        circle(sl, cx, 3.88, 0.3, fill='H' if el == 'H' else 'dark', line=None, content=el, size=12, color='white',
               bold=True)
        text(sl, cx + 0.22, 3.72, 1.1, 0.32, v, size=16)
    text(sl, x0, 4.16, 5.5, 0.6, 'Cu binds C and N far more weakly than Pt, H a little more weakly, O more strongly.',
         size=13, color='muted')
    text(sl, x0, 4.7, 5.5, 0.32, 'No shift for the empty site and 17 physisorbed species.', size=13, color='muted')
    caveat(sl, 0.55, 5.34, W - 2 * LM, 0.6,
           'R0 itself carries a nickel value: its rate matched RMG training reaction 1, "Deutschmann Ni" '
           '(A = 0.03, Ea = 5 kJ/mol, metal = Ni). Stage 10 checks the numbers that matter.', size=14)
    actions(sl, ['Use a library value when one exists.',
                 'Otherwise estimate the rate from templates and rate rules.',
                 'Take adsorbate thermo from Pt(111), or from the gas plus an adsorption estimate.',
                 'Shift it to Cu(111) with the binding energies.'])
    notes(sl, [f'Counted with a YAML parser from {RUN}/chem_annotated.yaml. "Template reaction:" 1,632 times',
               '(239 gas, 1,393 surface). "Library reaction: O2_dissociative_adsorption" 2 times:',
               'vacantX(8) + CO(4) <=> OCX(15) (line 3526, sticking A = 0.84) and 2 vacantX(8) + O2(6) <=> 2 OX(14)',
               '(line 4204, A = 0.07). Largest surface families: Surface_Abstraction_vdW 343, Surface_Abstraction 218,',
               'Surface_Abstraction_Beta_double_vdW 212; gas: Disproportionation 98, H_Abstraction 91.',
               'Surface thermo: 61 "Thermo library: surfaceThermoPt111 Binding energy corrected by LSR (...)"; 10 "Gas',
               'phase thermo for X from Thermo library: ... Adsorption correction: + Thermo group additivity estimation:',
               'adsorptionPt111(...)". Gas thermo (40): primaryThermoLibrary 14, thermo_DFT_CCSDTF12_BAC 23, DFT_QCI_thermo 3.',
               'Shifts per bond (LSR scaling 1.00) from the notes: H +16 ([PtH](33)), O -38 (OX(14)), C +199',
               '([Pt]=[CH][CH]=[Pt](73); OCX(15) has 0.50C -> +100), N +101 (N#[Pt](86)). They equal the Cu111 - Pt111',
               'binding energies of RMG-database metal.py: H -2.584 vs -2.754, C -4.960 vs -7.025, N -3.584 vs -4.632,',
               'O -4.208 vs -3.812 eV (1 eV = 96.485 kJ/mol). No shift (LSR ()): vacantX(8) and 17 vdW species.',
               'R0 note: "Matched reaction 1 H2 + Ni_3 + Ni_4 <=> HX_3 + HX_4 in Surface_Adsorption_Dissociative/training".',
               'That entry (RMG-database input/kinetics/families/Surface_Adsorption_Dissociative/training/reactions.py,',
               'index 1): StickingCoefficient A = 0.03, n = 0, Ea = 5 kJ/mol, shortDesc "Deutschmann Ni", metal = "Ni";',
               'its longDesc cites Delgado et al., Catalysts 2015, 5, 871-904 (as given in the database entry).',
               'Doc sections 1.4 and 10.4 (provenance and order of confidence).'])


def s09_yaml(d):
    sl = d.slide('What comes out: one YAML file with 111 species and 1,634 reactions',
                 'chem_annotated.yaml holds everything Cantera needs: phases, species and reactions.',
                 part=PB, tag='chem_annotated.yaml')
    xs0 = [0.55 + i * (2.87 + 0.25) for i in range(4)]
    card(sl, xs0[0], 1.66, 2.87, 1.82, 'Gas phase')
    text(sl, xs0[0] + 0.15, 2.08, 2.6, 0.4, '**40** species', size=18)
    text(sl, xs0[0] + 0.15, 2.5, 2.6, 0.4, '**239** gas reactions', size=18)
    text(sl, xs0[0] + 0.15, 2.95, 2.6, 0.4, 'ideal gas', size=14, color='muted')
    card(sl, xs0[1], 1.66, 2.87, 1.82, 'Surface phase')
    text(sl, xs0[1] + 0.15, 2.06, 2.6, 0.6, ['**71** species:', '53 one-site, 18 two-site'], size=15)
    text(sl, xs0[1] + 0.15, 2.7, 2.6, 0.6, ['**1,395** reactions:', '76 sticking, 1,319 Arrhenius'], size=15)
    card(sl, xs0[2], 1.66, 2.87, 1.82, 'Surface sites')
    text(sl, xs0[2] + 0.15, 2.08, 2.6, 0.3, 'site density Γ', size=14, color='muted')
    text(sl, xs0[2] + 0.15, 2.38, 2.6, 0.4, '2.943 × 10^{−5} mol/m²', size=17, bold=True)
    text(sl, xs0[2] + 0.15, 2.82, 2.6, 0.4, '= 2.943 × 10^{−8} kmol/m²', size=15)
    card(sl, xs0[3], 1.66, 2.87, 1.82, 'Also exported')
    text(sl, xs0[3] + 0.15, 2.08, 2.6, 0.62, "`seed/`: RMG's own copy of the core (Stage 2)", size=14)
    text(sl, xs0[3] + 0.15, 2.7, 2.6, 0.62, '`species/`: 111 pictures (path diagrams)', size=14)
    code(sl, 0.55, 3.62, 7.75, 2.05, ['- equation: vacantX(8) + vacantX(8) + H2(1) <=> [PtH](33) + [PtH](33)',
                                       '  type: sticking-Arrhenius',
                                       '  sticking-coefficient: {A: 0.03, b: 0.0, Ea: 5000.0}',
                                       '  note: |',
                                       '    Reaction index: Chemkin #240; RMG #22',
                                       '    Template reaction: Surface_Adsorption_Dissociative',
                                       '    Flux pairs: H2(1), [PtH](33); vacantX(8), [PtH](33); …',
                                       '    Matched reaction 1 H2 + Ni_3 + Ni_4 <=> HX_3 + HX_4 in …/training',
                                       '    This reaction matched rate rule [H2;VacantSite1;VacantSite2]'], size=12)
    line_y = [3.62 + 0.09 + 0.2 * i + 0.1 for i in range(9)]       # centre of each code line
    ann = [(0, 3.66, 'the reaction, as RMG wrote it (R0)', 'muted'),
           (2, 4.04, 'its rate: a sticking probability (J/mol)', 'muted'),
           (5, 4.42, 'where the number came from: a template', 'muted'),
           (6, 4.8, 'species that swap atoms (used in Stage 2)', 'muted'),
           (7, 5.18, 'matched to training reaction 1, on Ni', 'amber')]
    for i, yt, t, col in ann:
        arrow(sl, 8.95, yt + 0.15, 8.38, line_y[i], color='faint', lw=1.0)
        text(sl, 9.0, yt, 3.8, 0.3, t, size=14, color=col)
    actions(sl, ['The gas phase lists 40 species and 239 gas reactions.',
                 'The surface phase lists 71 adsorbates and the site density Γ.',
                 'Each species carries NASA polynomials for H, S and c_{p}.',
                 'Each reaction carries its rate parameters and where they came from.'])
    notes(sl, [f'{RUN}/chem_annotated.yaml: generator line 2 (RMG-Py CanteraWriter2, RMG-Py commit b342f9b); units',
               'line 5 {length: m, time: s, quantity: mol, activation-energy: J/mol}; phases from line 8 (gas line 9,',
               'surface line 21, site-density 2.943e-05 mol/m2 at line 42).',
               'Counted with a YAML parser: gas 40 species, 239 reactions; surface 71 species (53 with 1 site, 18 with 2),',
               '1,395 reactions (76 sticking-Arrhenius, 1,319 Arrhenius). Cantera says the same at load (run.log: "gas: 40',
               'species, 239 reactions; surface: 71 species, 1395 reactions; site density: 2.9430e-08 kmol/m2").',
               'Excerpt: lines 3494-3502, two lines shortened (…). Cantera writes R0 as H2(1) + 2 vacantX(8) <=> 2 [PtH](33).',
               'seed/: seed/seed (reactions.py, dictionary.txt), seed/filters and seed/seed_edge (seed_edge is kept on',
               'disk only, not in git). species/: 111 PNG pictures, one per species.'])


def s10_names(d):
    sl = d.slide("Reading RMG's names: X and [Pt] mean a surface site, even on copper",
                 'A name is an RMG label plus its RMG index. SMILES-style labels write the site atom as Pt, whatever the metal.',
                 part=PB, tag='chem_annotated.yaml')
    rows = [('H2(1)', 'H₂ in the gas. (1) is its RMG index.'),
            ('vacantX(8)', 'An empty surface site; vacantX is the label input.py gave it.'),
            ('[PtH](33)', 'H on one site: H∗. The site atom is written Pt, even on Cu.'),
            ('[H][H].[Pt](34)', 'H₂ resting on a site, not bonded to it (the dot): physisorbed H₂∗.'),
            ('OCX(15)', 'CO bonded to one site. X is RMG’s symbol for a site.'),
            ('[Pt]N=N[Pt](57)', 'N₂ bonded to two sites (two site atoms). keteneXX(434): XX means the same.')]
    y0, rh = 1.72, 0.7
    for i, (nm, mean) in enumerate(rows):
        y = y0 + i * rh
        if i % 2 == 0:
            shape(sl, 'rect', 0.55, y, W - 2 * LM, rh, fill='panel', line=None)
        text(sl, 0.75, y, 4.0, rh, nm, size=20, font=CODE, bold=True, anchor='m')
        cy = y + rh / 2
        cx = 5.6
        if i == 0:
            h2(sl, cx, cy, d=0.32)
        elif i == 1:
            shape(sl, 'rect', cx - 0.45, cy + 0.17, 0.9, 0.12, fill='ni', line=None)
            site(sl, cx, cy - 0.02, d=0.32)
        elif i == 2:
            shape(sl, 'rect', cx - 0.45, cy + 0.17, 0.9, 0.12, fill='ni', line=None)
            adsorbate(sl, cx, cy - 0.02, 'H', 'H', d=0.32)
        elif i == 3:
            shape(sl, 'rect', cx - 0.45, cy + 0.2, 0.9, 0.1, fill='ni', line=None)
            site(sl, cx, cy + 0.06, d=0.26)
            h2(sl, cx, cy - 0.2, d=0.22)
        elif i == 4:
            shape(sl, 'rect', cx - 0.45, cy + 0.17, 0.9, 0.12, fill='ni', line=None)
            adsorbate(sl, cx, cy - 0.02, 'CO', 'CO', d=0.32, size=12)
        else:
            shape(sl, 'rect', cx - 0.6, cy + 0.17, 1.2, 0.12, fill='ni', line=None)
            adsorbate(sl, cx, cy - 0.02, 'N2', 'N₂', d=0.32, n_sites=2, pitch=0.5, size=12)
        text(sl, 6.6, y, 6.1, rh, mean, size=16, anchor='m')
    actions(sl, ['Drop the trailing (index) to get the base label.',
                 'X, or Pt in a SMILES label, marks a bond to a site.',
                 'A dot means "on the site, but not bonded".',
                 'Two site marks mean the species takes two sites.'])
    notes(sl, ['docs/simulation_workflow.md Stage 0 (line 348): Cantera names are label(index); RMG writes its site atom X',
               f'as Pt in SMILES, whatever the metal. {RUN}/chem_annotated.yaml line 7 even gives element X the atomic',
               'weight 195.083 (Pt).',
               'All names are species of the mechanism; [Pt]N=N[Pt](57) and keteneXX(434) have "sites: 2".',
               'vacantX is the label of the empty site in input.py (line 70). Readable names come from settings.py',
               "species_labels (lines 32-42): 'empty site', 'H*', 'H2* (molecular)', 'CO* (1 site)'."])


SLIDES = [s06_input, s07_growth, s08_estimates, s09_yaml, s10_names]
