"""Section J (slides 44-47): Stage 10, can the mechanism be trusted? The binding check, the
switch-off what-ifs, the sources and sinks of a species, and where each number came from.
Running example: runs/try_2/Cu (none of Stage 10 is configured for it yet)."""
from deck_lib import *  # noqa: F401,F403

PJ = 'J · Can it be trusted? · Stage 10'
RUN = 'runs/try_2/Cu'


def s44_binding(d):
    sl = d.slide('Stage 10: check the mechanism’s binding energies against experiments',
                 'A binding energy 10 kJ/mol too strong makes that species’ coverage about 6 times too high at 650 K.',
                 part=PJ, tag='package · diagnostics.py:26–91')
    card(sl, 0.55, 1.66, 8.6, 3.55, 'One adsorption, three possible verdicts')
    X0, X1, Y = 1.0, 8.6, 2.85
    arrow(sl, X0, Y + 0.55, X1, Y + 0.55, color='muted', lw=1.25)
    text(sl, X0, Y + 0.62, X1 - X0, 0.3, 'ΔH_{ads} at 298 K (kJ/mol)   ← stronger binding', size=13, color='muted',
         align='c')
    shape(sl, 'rrect', 4.0, Y + 0.08, 2.0, 0.3, fill='rule', line=None, radius='pill')
    text(sl, 4.0, Y - 0.32, 2.0, 0.3, 'literature range', size=13, color='muted', align='c')
    for x, lab in [(2.3, 'TOO STRONG'), (5.0, 'OK'), (7.4, 'TOO WEAK')]:
        if lab == 'OK':
            circle(sl, x, Y + 0.23, 0.26, fill='ink', line=None)
        else:
            shape(sl, 'tri', x - 0.15, Y + 0.08, 0.3, 0.3, fill='dark', line=None, rotation=270 if lab == 'TOO STRONG' else 90)
        text(sl, x - 0.8, Y + 1.0, 1.6, 0.3, lab, size=13, bold=True, align='c', color='ink')
    schematic(sl, 7.6, 1.98)
    text(sl, 0.75, 4.42, 8.2, 0.3, 'One entry of diagnostics[\'adsorption_checks\'] in settings.py:', size=13,
         color='muted')
    text(sl, 0.75, 4.72, 8.3, 0.4, "`('H2 -> 2 H*', {'H2(1)': 1, 'vacantX(8)': 2}, {'[PtH](33)': 2}, (lo, hi), 'source')`",
         size=12)
    card(sl, 9.4, 1.66, 3.38, 3.55, 'Why it matters')
    text(sl, 9.55, 2.1, 3.15, 0.45, 'ΔH_{ads} = Σ_{k} ν_{k} H_{k}°(298 K)', size=15, font=MATH)
    text(sl, 9.55, 2.6, 3.15, 0.8, 'θ_{k}/θ_{∗} = exp(−ΔG_{ads}/RT) p_{k}/P°', size=15, font=MATH)
    text(sl, 9.55, 3.3, 3.15, 1.8, 'Too strong by δ raises θ_{k}/θ_{∗} by e^{δ/RT}: δ\u00a0=\u00a010\u00a0kJ/mol gives × 6.4 at '
         '650 K, × 2.7 at 1200 K.', size=14)
    caveat(sl, 0.55, 5.32, W - 2 * LM, 0.58, 'Not set for Cu yet: adsorption_checks = [] in try_2/Cu/settings.py. The '
           'ranges must be for Cu(111); Ni(111) values do not apply.', size=14)
    actions(sl, ['Compute each adsorption enthalpy from the mechanism’s own thermo.',
                 'Compare it with a literature range for Cu(111).',
                 'Verdict: OK, too strong or too weak.',
                 'Too strong by δ: coverage too high by e^{δ/RT}.'])
    notes(sl, ['analyses/diagnostics.py: reaction_thermo() (lines 26-36: ΔH and ΔG from each species\' NASA polynomials',
               'at T, reference states), binding_check() (39-58: OK if lo <= ΔH298 <= hi, TOO STRONG below, TOO WEAK above;',
               'ΔG at every main temperature), plot_binding_check() (61-91) -> mechanism_binding_check.csv, binding_check.png.',
               f'{RUN}/settings.py lines 66-70: adsorption_checks = [] with the comment that Cu(111) ranges are needed.',
               'The entry format follows config.py diagnostics (line 69): (label, {reactants}, {products}, (lo, hi) kJ/mol,',
               'note); the names in the example are real try_2/Cu species, lo and hi are placeholders.',
               'exp(10e3 / (8.314 x 650)) = 6.4; at 1200 K, 2.7. Doc Stage 10.1 (line 853). The plot is schematic.'])


def s45_whatif(d):
    sl = d.slide('What-if tests: switch reactions off with rate multipliers, never by editing energies',
                 'A multiplier scales forward and reverse together, so every equilibrium constant stays as it was.',
                 part=PJ, tag='package · diagnostics.py:94–195')
    card(sl, 0.55, 1.66, 6.25, 4.22, 'Example for try_2/Cu: switch off H atoms, H(16)')
    shape(sl, 'rrect', 1.0, 2.25, 1.3, 0.5, fill='white', line='H', lw=2, dash=True, radius='pill',
          content='H(16)', size=16, color='H', bold=True)
    arrow(sl, 2.3, 2.5, 3.1, 2.5, color='faint', lw=1.5, dash=True)
    shape(sl, 'rrect', 3.15, 2.27, 0.62, 0.46, fill='dark', line=None, radius='pill', content='× 0', size=13,
          color='white', bold=True)
    text(sl, 3.95, 2.3, 2.8, 0.4, '#26  H + ∗ ⇌ H∗', size=15, bold=True)
    text(sl, 1.0, 2.95, 5.6, 0.62, 'H(16) is in 1 surface reaction: H atoms adsorbing from, and desorbing to, the gas.',
         size=13, color='muted')
    text(sl, 0.75, 3.62, 5.9, 0.3, 'in settings.py, each entry is one what-if:', size=13, color='muted')
    forms = [("`diagnostics = {'switch_off': ['H(16)']}`", 'a species: all its surface reactions'),
             ("`[26]`", 'a reaction number'),
             ("`[['H(16)', 26]]`", 'a list, switched off together')]
    for i, (code_txt, what) in enumerate(forms):
        y = 3.95 + i * 0.5
        text(sl, 0.75, y, 3.6, 0.4, code_txt, size=12)
        text(sl, 4.45, y, 2.3, 0.45, what, size=12, color='muted')
    text(sl, 0.75, 5.45, 5.9, 0.35, 'Outputs: mechanism_switch_off_H_16.csv and switch_off_H_16.png', size=12,
         color='muted')
    card(sl, 7.05, 1.66, 5.73, 2.25, 'Why not edit the energy instead?')
    text(sl, 7.2, 2.1, 5.45, 1.75, ['It changes reverse rates by e^{δH/RT}; forward rates keep RMG’s old estimates.',
                                    'Barriers can turn negative, and the equations get stiffer.',
                                    'The real fix belongs in RMG: regenerate with corrected thermo.'], size=14,
         space_after=4)
    cases = ['Each what-if is solved again from a clean start; the direct solve fails on zero Jacobian rows, so the '
             'reactor method runs.',
             'Switch off one species at a time: the network can re-route. A sensitivity test, not a correction.']
    for i, t in enumerate(cases):
        y = 4.05 + i * 0.66
        shape(sl, 'rrect', 7.05, y, 5.73, 0.58, fill='white', line='rule', radius=0.08)
        text(sl, 7.15, y + 0.01, 5.55, 0.56, t, size=13, anchor='m')
    caveat(sl, 7.05, 5.38, 5.73, 0.5, 'Off by default (switch_off = []); not run for try_2/Cu.', size=13)
    actions(sl, ['Pick a suspect species or reaction.',
                 'Set every reaction that makes or uses it to zero.',
                 'Solve the main conditions again.',
                 'If the conclusions survive, they don’t hinge on it.'])
    notes(sl, ['analyses/diagnostics.py: switch_off_items() (lines 94-106: an entry is a species name, a reaction number, or',
               'a list of those), file_label() (109-112: H(16) -> H_16), switch_off_reactions() (115-131: clear errors for an',
               'unknown species or number), switch_off() (134-195: multipliers 0, clean start, the main cases at the',
               'diagnostics points), plot_switch_off() (217-263: (a) rate ratio, (b) slowdown normal dot / switched-off ring,',
               '(c) the equilibrium check q). Multipliers: Model.new_surface(multipliers) -> surf.set_multiplier (model.py',
               f'207-208). {RUN}: H(16) appears in exactly one surface reaction, #26 vacantX(8) + H(16) <=> [PtH](33).',
               'Defaults: config.py diagnostics (lines 68-76): switch_off [], switch_off_cases None (reference + every impurity).',
               'Docs: Stage 10.2 (line 867) and 10.5 (why thermochemistry is not edited in place).'])


def s46_sources(d):
    sl = d.slide('Trace who makes and uses a species, reaction by reaction',
                 'At a steady state each reaction adds ν q / Γ to a species. Sorting the terms shows which reactions matter.',
                 part=PJ, tag='package · diagnostics.py:266–289')
    card(sl, 0.55, 1.66, 6.0, 4.22, 'Sources and sinks of one species, e.g. H∗')
    cx, cy = 3.55, 3.75
    shape(sl, 'rrect', cx - 0.6, cy - 0.35, 1.2, 0.7, fill='H', line=None, radius='pill', content='H∗',
          size=18, color='white', bold=True)
    ins = [(1.0, 2.8, 3.5), (0.85, 3.75, 2.0), (1.0, 4.85, 1.0)]          # first chip clear of the label
    for x, y, w in ins:
        arrow(sl, x + 0.75, y, cx - 0.65, cy + (y - cy) * 0.25, color='ink', lw=w, head='sm')
        chip(sl, x - 0.35, y - 0.18, 'template', w=1.05, h=0.32, size=12)
    outs = [(6.0, 2.85, 2.5), (6.0, 4.65, 1.25)]
    for x, y, w in outs:
        arrow(sl, cx + 0.65, cy + (y - cy) * 0.25, x - 0.75, y, color='ink', lw=w, head='sm')
        chip(sl, x - 0.7, y - 0.18, 'template', w=1.05, h=0.32, size=12)
    text(sl, 0.75, 2.16, 2.0, 0.3, 'producers', size=13, color='muted', bold=True)
    text(sl, 5.0, 2.16, 1.5, 0.3, 'consumers', size=13, color='muted', bold=True, align='r')
    schematic(sl, 0.75, 5.45, w=5.5, label='schematic widths · H∗ is in 126 surface reactions here')
    card(sl, 6.8, 1.66, 5.98, 3.0, 'What it lists')
    text(sl, 6.98, 2.1, 5.65, 2.5, ['the top_n (5) largest makers (+) and users (−) of the species, per site (s^{−1}),',
                                     'with the reaction number, its equation and its provenance (the YAML note).',
                                     'Setting: diagnostics[\'top_reactions\'] = [(case, T, species)] or (case, T, P, species).',
                                     'Output: top_reactions.csv.'], size=13, space_after=5)
    caveat(sl, 6.8, 4.8, 5.98, 1.08, "Not set for try_2/Cu (top_reactions = []). For example ('CO', 650.0, '[PtH](33)') "
           'would trace H∗ in the CO case at 650 K.', size=13)
    actions(sl, ['Pick a species and a steady state.',
                 'List every reaction’s contribution to the species.',
                 'Sort producers and consumers by size.',
                 'Read each top reaction’s provenance note.'])
    notes(sl, ['analyses/diagnostics.py top_reactions() (lines 266-289): at a main steady state (P None = reference_pressure),',
               'q = net rates of progress / Γ; for each reaction with the species, the contribution ν q; sorted by |ν q|;',
               'the source is the second line of the reaction\'s YAML note. Columns of top_reactions.csv: case, T_K,',
               'P_ref_bar, species, label, coverage, rate_1_per_s, reaction, equation, source.',
               f'{RUN}/settings.py has no top_reactions entry; the example entry is illustrative (real species name).',
               '[PtH](33) appears in 126 surface reactions (YAML parse). Docs: Stage 10.3. The arrows are schematic.'])


def s47_provenance(d):
    sl = d.slide('Where each number came from: in try_2/Cu, none is from a Cu(111) library',
                 'The order of confidence, filled in for this mechanism. Trust library values most, scaled estimates least.',
                 part=PJ, tag='chem_annotated.yaml')
    ladder = [('amber', '1', 'a library value for this metal and facet',
               '2 reactions, from O2_dissociative_adsorption (CO and O₂ adsorption)',
               'that library is not in the local RMG-database: its metal is not recorded here', 2, 1634, 'reactions'),
              ('faint', '2', 'training reactions and rate rules',
               '1,632 reactions estimated from templates, R0 among them',
               'R0 matched a Ni training reaction, “Deutschmann Ni” (A = 0.03, Ea = 5 kJ/mol)', 1632, 1634, 'reactions'),
              ('amber', '3', 'values scaled from another metal',
               'the thermochemistry of all 71 surface species: Pt(111), shifted to Cu(111)',
               'by linear scaling with the Cu(111) binding energies (slide 8)', 71, 71, 'surface species')]
    for i, (c, n, h, sub, flag, k, tot, unit) in enumerate(ladder):
        y = 1.72 + i * 1.32
        card(sl, 0.55, y, W - 2 * LM, 1.2, None, fill='panel')
        circle(sl, 1.0, y + 0.6, 0.48, fill=c, line=None, content=n, size=17, color='white', bold=True)
        text(sl, 1.45, y + 0.1, 6.6, 0.35, h, size=16, bold=True)
        text(sl, 1.45, y + 0.46, 6.9, 0.35, sub, size=14)
        text(sl, 1.45, y + 0.8, 6.9, 0.35, flag, size=13, color='amber')
        bx, bw = 8.75, 3.7
        shape(sl, 'rect', bx, y + 0.42, bw, 0.32, fill='white', line='rule', lw=0.75)
        shape(sl, 'rect', bx, y + 0.42, max(bw * k / tot, 0.04), 0.32, fill='dark' if i else 'amber', line=None)
        text(sl, bx, y + 0.78, bw, 0.3, f'{k:,} of {tot:,} {unit}', size=12, color='muted', align='r')
    text(sl, 0.55, 5.72, W - 2 * LM, 0.3, 'Look hardest at R0 (a Ni value on Cu) and at the species that hold the most '
         'sites (Stage 8): compare them with Cu(111) literature (slide 44).', size=14, color='amber', bold=True)
    actions(sl, ['Read each number’s provenance note in the YAML.',
                 'Library values for this metal come first.',
                 'Then rate rules and training reactions.',
                 'Values scaled from another metal come last.'])
    notes(sl, ['Docs Stage 10.4 (provenance; order of confidence) and section 1.4.',
               f'Counts from {RUN}/chem_annotated.yaml (YAML parser): "Library reaction: O2_dissociative_adsorption" 2 times',
               '(vacantX(8) + CO(4) <=> OCX(15), line 3526; 2 vacantX(8) + O2(6) <=> 2 OX(14), line 4204); "Template reaction:"',
               '1,632 times. Surface thermo: 61 surfaceThermoPt111 + LSR, 10 gas-phase thermo + adsorptionPt111 + LSR; none',
               'from a Cu library (input.py thermoLibraries: surfaceThermoPt111, primaryThermoLibrary, thermo_DFT_CCSDTF12_BAC,',
               "DFT_QCI_thermo, 'auto').",
               'O2_dissociative_adsorption: searched in all four local RMG-database copies (input/kinetics/libraries): not',
               'found; it was probably a library on the cluster that ran try_2/Cu. Ask for its source before trusting it.',
               'R0: training entry 1 of Surface_Adsorption_Dissociative, shortDesc "Deutschmann Ni", metal "Ni" (slide 8).'])


SLIDES = [s44_binding, s45_whatif, s46_sources, s47_provenance]
