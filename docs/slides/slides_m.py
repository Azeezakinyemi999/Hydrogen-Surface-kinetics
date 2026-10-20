"""Section M (slides 62-65): limits and next: the model's assumptions, what the workflow does not do
yet, where the twelve Cu, Fe and Ni runs stand, and the whole workflow on one page.
Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PM = 'M · Limits and next'
RUN = 'runs/try_2/Cu'


def s62_assumptions(d):
    sl = d.slide('What the model assumes, and what that means for try_2/Cu',
                 'Every number in this deck holds only inside these seven assumptions.',
                 part=PM, tag='docs · simulation_workflow.md §1.4')
    text(sl, 0.55, 1.62, 2.4, 0.3, 'assumption', size=12, color='muted', bold=True)
    text(sl, 3.05, 1.62, 4.6, 0.3, 'what it says', size=12, color='muted', bold=True)
    text(sl, 7.85, 1.62, 4.9, 0.3, 'for try_2/Cu', size=12, color='muted', bold=True)
    rows = (('Mean field', 'rates depend only on average coverages: no islands, no ordering',
             'no coverage-dependent parameters in any of the 12 YAMLs, so no lateral interactions either'),
            ('One site type', 'an ideal surface with a constant site density Γ',
             'Γ = 2.943 × 10^{−5} mol/m²; every species sits on the same Cu(111) site'),
            ('Isothermal', 'heats of reaction do not change the temperature',
             'every reactor runs with its energy equation off'),
            ('Intrinsic kinetics', 'no mass-transfer limits; the gas is uniform',
             'the rates are those of the surface alone; a real bed can be slower'),
            ('Fixed gas', 'the surface does not change the gas: negligible conversion',
             '1 m³ of gas over 10^{−10} m² of catalyst, gas reactions off'),
            ('Parameters as given', 'RMG’s estimates: libraries, rate rules, scaling from other metals',
             'R0 is a Ni training value; all 71 surface species scaled from Pt(111) (slide 47)'),
            ('Cantera’s site counting', 'species on several sites follow Cantera’s convention',
             'the site-count factor f enters Π (slides 13 and 49)'))
    for i, (h, says, here) in enumerate(rows):
        y = 1.98 + i * 0.68
        line(sl, 0.55, y - 0.04, W - LM, y - 0.04, color='rule', lw=0.75)
        text(sl, 0.55, y, 2.4, 0.6, h, size=14, bold=True, anchor='m')
        text(sl, 3.05, y, 4.6, 0.6, says, size=14, anchor='m')
        text(sl, 7.85, y, 4.93, 0.6, here, size=13, anchor='m', color='ink')
    line(sl, 0.55, 1.98 + 7 * 0.68 - 0.04, W - LM, 1.98 + 7 * 0.68 - 0.04, color='rule', lw=0.75)
    notes(sl, ['Docs section 1.4 "Model assumptions and limits" (simulation_workflow.md line 173): mean field, ideal surface,',
               'isothermal, intrinsic kinetics, fixed gas means negligible conversion, parameters as given, Cantera\'s site',
               'counting (section 2.1).',
               'For try_2/Cu: grep "coverage-dependencies" finds 0 in each of the 12 runs/try_*/*/chem_annotated.yaml.',
               f'Site density: {RUN}/chem_annotated.yaml line 42 (site-density 2.943e-05, mol/m2 as on slide 9). Energy off:',
               'solver.py fixed_gas_reactor (lines 98-107, IdealGasReactor energy="off") and analyses/flux.py simulate.',
               'Fixed gas: config.py reactor_area_per_volume 1e-10 m2 per m3 (line 22), gas multiplier 0 (slide 5).',
               'Provenance: slide 47 (chem_annotated.yaml notes). Site counting: slide 13; Π with f: slide 49.'])


def s63_not_yet(d):
    sl = d.slide('Not in the workflow yet: four things it does not do',
                 'Each would answer a question the steady states cannot, and each can start from a piece that exists.',
                 part=PM, tag='concept')
    cards = (('A flow reactor', 'Along a real bed the gas changes as reactions convert it. The fixed gas holds only at '
                                'the inlet, or at low conversion.',
              'Start from: fixed-gas steps in a chain, the gas updated between them.'),
             ('TPD', 'Temperature-programmed desorption: heat a covered surface and record when each species leaves. '
                     'The peak temperatures test the binding energies directly.',
              'Start from: the transient reactor of slide 58, with a temperature ramp and no gas.'),
             ('Uncertainty', 'RMG’s estimates have errors, and nothing carries them through to the rates and '
                             'coverages.',
              'Start from: TRC and DRC (slides 55–56), which rank the estimates that matter most.'),
             ('Coverage effects', 'No lateral interactions: a crowded adsorbate binds as strongly as a lone one, so '
                                  'high coverages can come out too high.',
              'Start from: coverage-dependent parameters in the YAML, which Cantera already applies.'))
    for i, (h, what, start) in enumerate(cards):
        x, y = (0.55, 6.78)[i % 2], (1.66, 4.27)[i // 2]
        card(sl, x, y, 6.0, 2.45, h, size=17)
        text(sl, x + 0.2, y + 0.55, 5.6, 1.0, what, size=14)
        shape(sl, 'rrect', x + 0.17, y + 1.62, 5.66, 0.66, fill='white', line='rule', radius=0.08)
        text(sl, x + 0.3, y + 1.64, 5.45, 0.62, start, size=13, anchor='m')
    notes(sl, ['Not implemented in the package (no flow-reactor, ramp, uncertainty or lateral-interaction code in',
               'hydrogen_surface_kinetics/).',
               'Flow reactor: docs section 1.2 (line 108): the fixed gas represents a differential reactor or the inlet of',
               'a flow reactor; section 1.4: fixed gas means negligible conversion.',
               'TPD: would reuse transients.run_transient (lines 30-56: set_state, fixed_gas_reactor, net.advance) with T',
               'stepped in time and the gas emptied; compare peak temperatures with the binding check of Stage 10.1.',
               'Uncertainty: Stage 10 judges parameters but propagates nothing; sensitivity.csv (Stage 11.4) ranks reactions',
               '(DRC) and species (TRC).',
               'Coverage effects: docs section 1.4 (lateral interactions only through coverage-dependent parameters) and',
               'section 2.2 (line 249: Cantera applies them if the YAML has them); none of the 12 YAMLs has any.'])


def s64_runs(d):
    sl = d.slide('The same workflow on Cu, Fe and Ni: where each run stands',
                 'As of Oct 20, from each run.log. The Fe and Ni mechanisms have 5 to 17 times the surface reactions of try_2/Cu.',
                 part=PM, tag='runs/try_*/')
    tries = ('try_2', 'try_3', 'try_4', 'try_5')
    report = ('okbg', 'ok')
    bad = ('amberbg', 'amber')
    run_ = ('white', 'H')
    plog = 'Cantera rejects the PLOG rate of C₂H₃ + CH₂O ⇌ C₂H₄ + HCO'
    grid_ = {
        'Cu': [(report, 'report, Oct 15', '133 pages', '1,395'), (report, 'report, Oct 16', '134 pages', '1,197'),
               (report, 'report, Oct 15', '139 pages', '4,604'), (report, 'report, Oct 16', '134 pages', '1,369')],
        'Fe': [(run_, 'running', 'poisoning curves: 636 of 680 done', '23,334'),
               (bad, 'does not load', plog, '13,756'), (bad, 'does not load', plog, '15,332'),
               (bad, 'does not load', plog, '6,339')],
        'Ni': [(report, 'report, Oct 18', '136 pages', '17,037'),
               (bad, 'killed at load', 'exit 137 after 3.5 min; 17,430 gas reactions', '18,771'),
               (bad, 'killed at load', 'exit 137 after 2.7 min; 28,532 gas reactions', '19,741'),
               (bad, 'does not load', plog, '10,892')]}
    X0, CW, GAP = 1.5, 2.7, 0.12
    for j, t in enumerate(tries):
        text(sl, X0 + j * (CW + GAP), 1.62, CW, 0.3, t, size=14, color='muted', bold=True, align='c')
    for i, metal in enumerate(('Cu', 'Fe', 'Ni')):
        y = 1.98 + i * 1.33
        text(sl, 0.55, y, 0.85, 1.22, metal, size=20, bold=True, anchor='m', font=TITLE)
        for j, ((fill, edge), status, detail, size) in enumerate(grid_[metal]):
            x = X0 + j * (CW + GAP)
            shape(sl, 'rrect', x, y, CW, 1.22, fill=fill, line=edge, lw=1.25, radius=0.08)
            text(sl, x + 0.12, y + 0.06, CW - 0.24, 0.32, status, size=14, bold=True,
                 color='H' if fill == 'white' else 'ink')
            text(sl, x + 0.12, y + 0.38, CW - 0.24, 0.55, detail, size=12)
            text(sl, x + 0.12, y + 0.9, CW - 0.24, 0.28, f'{size} surface reactions', size=12, color='muted')
    actions(sl, ['Resume the Cu batch: paused since Oct 18 at 23:02.',
                 'Repair that PLOG rate so the four failing runs load.',
                 'Find why try_3/Ni and try_4/Ni are killed at load.',
                 'Then compare the metals at the same 72 conditions.'])
    notes(sl, ['Status from runs/try_*/*/run.log and output/report/ on 2026-10-20:',
               'Cu: reports built 2026-10-15 23:49 (try_2, 133 pages), 10-16 00:05 (try_3, 134), 10-15 22:07 (try_4, 139),',
               '10-16 00:21 (try_5, 134). The update batch "bash run.sh runs/try_2/Cu runs/try_3/Cu runs/try_4/Cu',
               'runs/try_5/Cu" started 2026-10-18 23:02:07 is stopped (ps state T, 0.1 s CPU in 37 h): resume it with fg in',
               'its terminal. try_2/Ni: exit 0 on 2026-10-18 21:54, report 136 pages.',
               'try_2/Fe: running since 2026-10-18 23:02:20; at curves 636/680 on 2026-10-20 (one TIMEOUT after 1800 s).',
               'try_3/Fe, try_4/Fe, try_5/Fe (2026-10-15) and try_5/Ni (2026-10-18): exit 1 at load, CanteraError from',
               'PlogRate::validate: "Invalid rate coefficient for reaction C2H3(145) + CH2O(83) <=> C2H4(12) + HCO(89)" (try_3/Fe,',
               'chem_annotated.yaml line 4659; try_5/Ni: C2H3(156) + CH2O(94) <=> C2H4(12) + HCO(100), line 3825). RMG\'s own',
               'Cantera check flagged the same reaction (Surface-Impurity-Effect commits 19b3ebb and 76d0cf1).',
               'try_3/Ni and try_4/Ni (2026-10-18): exit 137 (killed) 3.5 and 2.7 min after starting, before the model summary;',
               'the system log shows no memory kill. Sizes: surface reactions from each YAML (surface-reactions list); gas',
               'reactions 17,430 and 28,532. All 12 settings.py use the same 8 points (temperatures, pressures,',
               'pressure_temperatures).'])


def s65_summary(d):
    sl = d.slide('One picture: from an RMG mechanism to a report',
                 'The whole workflow on one page, with try_2/Cu’s numbers.',
                 part=PM, tag='summary')
    steps = (('RMG', 'RMG builds the mechanism', '15 reactors, 650–1200 K; 1,395 surface reactions'),
             ('0–1', 'Load it, name the roles', 'checks at load; R0 and the 7 impurities found by name'),
             ('2', 'Replay RMG’s reactors', 'optional: flux-diagram videos at 650 and 1200 K'),
             ('3', 'Design the gas', '8 points × 9 cases = 72 conditions, p(H₂) held fixed'),
             ('4–6', 'Solve the steady states', 'two starts each; saved under a fingerprint'),
             ('7–9', 'Measure', 'TOF and Π, per site and per area; who covers the surface'),
             ('10', 'Judge the mechanism', 'binding check, what-ifs, provenance; none set for Cu yet'),
             ('11', 'Go further', '9 analyses, 4,534 more solves'),
             ('PDF', 'Read the report', 'every table and figure, each with how to read it'))
    for i, (tag_, h, what) in enumerate(steps):
        y = 1.72 + i * 0.57
        optional = tag_ == '2'
        shape(sl, 'rrect', 0.55, y + 0.05, 0.8, 0.36, fill='white' if optional else 'dark',
              line='dark' if optional else None, lw=1.0, dash=optional, radius='pill', content=tag_, size=12,
              color='dark' if optional else 'white', bold=True)
        if i < len(steps) - 1:
            line(sl, 0.95, y + 0.41, 0.95, y + 0.62, color='faint', lw=1.25)
        text(sl, 1.5, y, 3.3, 0.46, h, size=14, bold=True, anchor='m')
        text(sl, 4.85, y, 4.35, 0.46, what, size=13, anchor='m')
    card(sl, 9.4, 1.66, 3.38, 5.1, 'try_2/Cu in numbers')
    nums = (('72', 'main conditions'), ('4,606', 'saved solves, all converged'), ('21 + 61', 'tables and figures'),
            ('133', 'pages of report'))
    for i, (n, lab) in enumerate(nums):
        y = 2.2 + i * 1.1
        text(sl, 9.6, y, 3.0, 0.5, n, size=26, bold=True, font=TITLE)
        text(sl, 9.6, y + 0.52, 3.0, 0.4, lab, size=13, color='muted')
    notes(sl, [f'Counts for {RUN}: input.py 15 reactors at 650-1200 K (slide 6); chem_annotated.yaml 71 surface species and',
               '1,395 surface reactions (slide 9); 7 impurities (slide 21); 72 main conditions (slide 27); two starts and the',
               'fingerprint cache (slides 36-38); Format 2 (slide 18); Stage 10 not set for Cu (slides 44-46).',
               'Stage 11: curves 680 + sweep 108 + pressure sweep 162 + orders 304 + rate control 3,056 + transients 224 =',
               '4,534 solves (arrhenius, equilibrium and paths solve nothing); with the 72 main ones, 4,606 saved solves =',
               'cache/ 1,326 + cache_sensitivity/ 3,056 + cache_transients/ 224; status: 0 not converged.',
               'Stage 2 is optional and not run for any try_* mechanism (slide 24). 21 tables and 61 figures in output/,',
               'report 133 pages (slides 60-61).'])


SLIDES = [s62_assumptions, s63_not_yet, s64_runs, s65_summary]
