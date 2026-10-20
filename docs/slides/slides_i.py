"""Section I (slides 39-43): Stages 7-9, what is measured: forward and net rates, the poisoning
factor and the site-blocking check, the results table, the surface composition and the rate
plots. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PI = 'I · What is measured · Stages 7–9'
RUN = 'runs/try_2/Cu'


def s39_forward_net(d):
    sl = d.slide('Stage 7: the forward rate of R0 measures splitting; its net rate is nearly zero',
                 'In pure H₂, splitting and recombining balance: the step runs fast in both directions.',
                 part=PI, tag='package · model.py:212–218')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'Pure H₂ at steady state')
    slab(sl, 0.95, 4.6, 8, occ={i: ('H', 'H', 1) for i in (0, 1, 2, 4, 5, 6)})
    h2(sl, 2.2, 2.45)
    h2(sl, 4.6, 2.45)
    arrow(sl, 2.2, 2.7, 2.2, 4.2, color='H', lw=4)
    arrow(sl, 4.6, 4.2, 4.6, 2.7, color='H', lw=4)
    text(sl, 0.75, 3.2, 1.35, 0.6, 'splitting (forward)', size=13, color='H', align='r')
    text(sl, 4.75, 3.2, 1.7, 0.6, 'recombining (reverse)', size=13, color='H')
    text(sl, 0.75, 5.3, 5.6, 0.4, 'Equal widths: R0 is at equilibrium.', size=14, color='muted', align='c')
    card(sl, 6.83, 1.66, 5.95, 4.22, 'What the workflow reports')
    defs = [('TOF_{f} = q_{f} / Γ', 'splitting events per site per second: tof_fwd'),
            ('TOF_{r} = q_{r} / Γ', 'recombinations per site per second: tof_rev'),
            ('net TOF = TOF_{f} − TOF_{r}', 'close to zero in pure H₂: tof_net'),
            ('−ṡ_{H₂} / Γ', 'all H₂ use per site, R1 (H₂∗) included: tof_probe_consumed')]
    for i, (e, t) in enumerate(defs):
        y = 2.15 + i * 0.85
        text(sl, 7.0, y, 5.6, 0.4, e, size=18, font=MATH)
        text(sl, 7.0, y + 0.4, 5.6, 0.35, t, size=13, color='muted')
    actions(sl, ['Count splitting events per site per second: TOF_{f}.',
                 'Count recombinations the same way: TOF_{r}.',
                 'In pure H₂ they nearly cancel: R0 is at equilibrium.',
                 'So the forward TOF is reported as the splitting rate.'])
    notes(sl, ['docs/simulation_workflow.md section 2.5 (forward vs net) and Stage 7 (line 733).',
               'Model.probe_rates() (model.py 212-218) swaps forward and reverse when the file writes R0 backwards',
               '(sign -1), so "forward" always means splitting. The columns are computed by analyses/common.py evaluate()',
               '(lines 136-161): tof_fwd, tof_rev, tof_net, tof_probe_consumed, and tof_molec_fwd / tof_molec_net for R1,',
               f'H2(1) + vacantX(8) <=> [H][H].[Pt](34) in {RUN}. The picture is schematic.'])


def s40_poisoning(d):
    sl = d.slide('The poisoning factor: the forward rate with the impurity over the rate without',
                 'If the impurity acts only by taking sites, the empty sites alone predict the same ratio.',
                 part=PI, tag='package · common.py:217–233')
    occR = {(0, 0): ('H', 'H', 1), (0, 1): ('H', 'H', 1), (0, 3): ('H', 'H', 1), (1, 0): ('H', 'H', 1),
            (1, 2): ('H', 'H', 1), (1, 3): ('H', 'H', 1), (2, 1): ('H', 'H', 1), (2, 2): ('H', 'H', 1)}
    occI = {(0, 0): ('CO', 'CO', 1), (0, 1): ('H', 'H', 1), (0, 2): ('CO', 'CO', 1), (0, 3): ('CO', 'CO', 1),
            (1, 0): ('CO', 'CO', 1), (1, 1): ('CO', 'CO', 1), (1, 3): ('H', 'H', 1), (2, 0): ('CO', 'CO', 1),
            (2, 1): ('CO', 'CO', 1), (2, 2): ('CO', 'CO', 1), (2, 3): ('H', 'H', 1)}
    text(sl, 0.55, 1.68, 2.6, 0.32, 'reference: H₂ only', size=14, bold=True)
    grid(sl, 0.75, 2.15, 4, 3, pitch=0.5, d=0.4, occ=occR)
    text(sl, 0.55, 3.75, 2.8, 0.3, 'θ_{∗,ref}: 4 of 12 empty', size=14, font=MATH)
    text(sl, 3.35, 1.68, 2.6, 0.32, 'with CO', size=14, bold=True)
    grid(sl, 3.55, 2.15, 4, 3, pitch=0.5, d=0.4, occ=occI)
    text(sl, 3.35, 3.75, 2.8, 0.3, 'θ_{∗}: 1 of 12 empty', size=14, font=MATH)
    text(sl, 0.55, 4.2, 5.6, 0.3, 'example surfaces (not results)', size=12, color='faint', italic=True)
    card(sl, 6.45, 1.66, 6.33, 2.75, None)
    text(sl, 6.6, 1.8, 6.1, 0.45, 'Π = TOF_{f} / TOF_{f,ref}        slowdown = 1 / Π', size=20, font=MATH)
    text(sl, 6.6, 2.42, 6.1, 0.35, 'site-blocking prediction, for a step needing n empty sites:', size=14,
         color='muted')
    text(sl, 6.6, 2.8, 6.1, 0.5, 'Π = (p_{H₂} / p_{H₂,ref}) (θ_{∗} / θ_{∗,ref})^{n} (f / f_{ref})^{n}', size=19,
         font=MATH)
    text(sl, 6.6, 3.45, 6.1, 0.8, 'R0: n = 2, read from the mechanism. p(H₂) is the same in both (Stage 3), so the '
         'first factor is 1.', size=14, color='muted')
    text(sl, 0.55, 4.62, 5.7, 0.5, 'blocking_check = Π / prediction  →  should be 1', size=18, font=MATH)
    caveat(sl, 6.45, 4.62, 6.33, 1.26, 'A blocking_check other than 1 points to coverage-dependent parameters, an '
           'indexing error, or a probe step that is not elementary. The run lists any above 10^{−6}.', size=14)
    actions(sl, ['Solve the H₂-only reference at the same point.',
                 'Divide the forward TOF with the impurity by the reference.',
                 'Predict the same ratio from empty sites and f alone.',
                 'If they match, the impurity works only by taking sites.'])
    notes(sl, ['analyses/common.py add_poisoning() (lines 217-233): the reference is matched by (basis, T_K, P_ref_bar);',
               'poisoning_factor = tof_fwd / tof_fwd_ref, slowdown_x = 1 / poisoning_factor, poisoning_factor_B for start B;',
               'blocking_from_theta = (θ_∗/θ_∗,ref)^n (p/p_ref); blocking_full = blocking_from_theta (f/f_ref)^n;',
               'blocking_check = poisoning_factor / blocking_full. n = Model.probe_sites (2, found by find_role_reaction).',
               'steady_state.py run(): lists conditions not converged, with starts that disagree, or |blocking_check - 1| > 1e-6',
               '(line 38). Doc Stage 7: the site-blocking derivation, TOF_f = k_f C_A2 (f θ_∗ Γ)^n / Γ, and Figure 6.',
               'The two grids are illustrations (4/12 and 1/12 empty), not results.'])


def s41_table(d):
    sl = d.slide('steady_state_results.csv: one row per case and point, everything Stage 7 measures',
                 f'The real {len(_COLS_FLAT)} columns of {RUN}/output/steady_state_results.csv, grouped by what they are for.',
                 part=PI, tag='package · steady_state.py')
    for i, (head, cols, purpose) in enumerate(_GROUPS):
        row, col = divmod(i, 4)
        x = 0.55 + col * 3.1
        y = 1.66 + row * 2.45
        w = 2.95
        card(sl, x, y, w, 2.32, head, size=15)
        text(sl, x + 0.12, y + 0.45, w - 0.24, 1.25, ', '.join(cols), size=12, font=CODE)
        text(sl, x + 0.12, y + 1.8, w - 0.24, 0.5, purpose, size=13, color='muted')
    text(sl, 9.95, 4.2, 2.85, 2.1, 'Written by steady_state.py from common.results_frame; the same columns, '
         'per point and level, fill the tables of Stage 11.', size=13, color='muted')
    notes(sl, [f'Column list: the header of {RUN}/output/steady_state_results.csv ({len(_COLS_FLAT)} columns).',
               'analyses/steady_state.py build_summary() (lines 19-24) and run() (27-60); analyses/common.py results_frame()',
               '(lines 192-214, which also adds the Format 2 columns through per_area.add_columns) and add_poisoning()',
               '(217-233). Docs: Appendix C (output files and their columns) and Table R3.'])


_GROUPS = [
    ('the condition', ['key', 'basis', 'case', 'T_K', 'P_ref_bar', 'series', 'P_total_bar', 'P_probe_bar',
                       'P_impurity_bar'], 'which gas, at which point'),
    ('rates per site (s⁻¹)', ['tof_fwd', 'tof_rev', 'tof_net', 'tof_probe_consumed', 'tof_molec_fwd', 'tof_molec_net'],
     'R0 forward, reverse, net; all H₂ use; R1'),
    ('coverages', ['theta_free', 'theta_probe', 'theta_molec', 'theta_impurity', 'site_count_factor', 'top_species'],
     'θ_{∗}, θ_{H}, θ_{H₂∗}, θ_{imp}, f, the top three'),
    ('start B', ['tof_fwd_B', 'theta_free_B', 'theta_impurity_B', 'poisoning_factor_B', 'starts_agree'],
     'the second steady state, if any'),
    ('the solver', ['converged', 'start_A_solved_by', 'start_B_solved_by', 'residual_1_per_s'],
     'did it converge, and how'),
    ('per area (Format 2)', ['J_diss_mol_m2_s', 'J_rec_mol_m2_s', 'J_net_mol_m2_s', 'J_exchange_mol_m2_s',
                             'kd_eff_mol_m2_s_Pa', 'gamma_eff', 'C_H_mol_m2'], 'the same rates per m²'),
    ('poisoning', ['poisoning_factor', 'slowdown_x', 'blocking_from_theta', 'blocking_full', 'blocking_check'],
     'Π, slowdown, the site-blocking check'),
]
_COLS_FLAT = [c for _, cols, _ in _GROUPS for c in cols]


def s42_composition(d):
    sl = d.slide('Stage 8: who holds the sites, with one colour per species in every figure',
                 'Each bar is one case at one point; its segments are the coverages, and they add up to 1.',
                 part=PI, tag='package · composition.py')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'surface_composition_vs_T: one panel per point')
    base, top = 4.95, 2.3
    bars = [('H₂ only', [('panel2', 0.15), ('H', 0.8), ('faint', 0.05)]),
            ('H₂O', [('panel2', 0.1), ('H', 0.45), ('H2O', 0.4), ('faint', 0.05)]),
            ('CO', [('panel2', 0.05), ('H', 0.2), ('CO', 0.7), ('faint', 0.05)]),
            ('all', [('panel2', 0.05), ('H', 0.15), ('H2O', 0.3), ('CO', 0.45), ('faint', 0.05)])]
    for i, (lab, segs) in enumerate(bars):
        x = 0.95 + i * 0.95
        y = base
        for c, v in segs:
            h = v * (base - top)
            shape(sl, 'rect', x, y - h, 0.6, h, fill=c, line='white', lw=0.75)
            y -= h
        text(sl, x - 0.2, base + 0.05, 1.0, 0.3, lab, size=12, align='c')
    leg = [('panel2', 'empty site'), ('H', 'H*'), ('H2O', 'H2O*'), ('CO', 'CO* (1 site)'), ('faint', 'other')]
    for i, (c, t) in enumerate(leg):
        shape(sl, 'rect', 4.82, 2.35 + i * 0.42, 0.24, 0.24, fill=c, line=None)
        text(sl, 5.12, 2.28 + i * 0.42, 1.35, 0.36, t, size=12)
    schematic(sl, 0.7, 5.42, w=3.0, label='layout only, not results')
    card(sl, 6.83, 1.66, 5.95, 4.22, 'The rules (settings: composition)')
    rules = [('listed', 'every species above 10^{−3}: steady_state_surface_species.csv'),
             ('own segment', 'if it reaches 2 % anywhere; at most 7 named, the rest “other”'),
             ('names', 'readable names from species_labels in settings.py, e.g. “CO* (1 site)”'),
             ('colours', 'one per species, the same in every figure of the run'),
             ('figures', 'surface_composition_vs_T and _vs_P: one panel per point')]
    for i, (h, t) in enumerate(rules):
        y = 2.12 + i * 0.72
        text(sl, 7.0, y, 1.55, 0.35, h, size=14, bold=True)
        text(sl, 8.55, y, 4.1, 0.66, t, size=13)
    actions(sl, ['Collect every coverage at every main condition.',
                 'List the species above 10⁻³, with readable names.',
                 'Stack the coverages per case; lump the small ones as “other”.',
                 'Keep each species’ colour across all figures.'])
    notes(sl, ['analyses/composition.py: coverage_table() (line 20), surface_species() (30), species_palette() (39: a fixed',
               'colour per species), run() (63), plot_surface_composition() (88). Settings: config.py composition (lines',
               '63-67): min_listed 1e-3, min_plotted 0.02, max_named 7.',
               f'Legend names are the {RUN}/settings.py species_labels (lines 32-42): vacantX(8) "empty site",',
               '[PtH](33) "H*", O.[Pt](39) "H2O*", OCX(15) "CO* (1 site)". Docs: Stage 8 (line 789).',
               'The bar heights are a layout sketch, not results.'])


def s43_rate_plots(d):
    sl = d.slide('Stage 9: compare rates on log axes, one figure per series',
                 'Rates are dots, not bars, because they span decades. One marker per temperature (or pressure).',
                 part=PI, tag='package · steady_state.py:63–131')
    marks = [('oval', 'ink', '650 K'), ('rect', 'muted', '800 K'), ('tri', 'faint', '1000 K'),
             ('diamond', 'niline', '1200 K')]                     # neutral greys; the shape tells T apart
    cases = ['H₂ only', 'H₂O', 'CO', 'all']
    panels = [('(a) forward TOF (log)', [[0.8, 0.82, 0.84, 0.85], [0.78, 0.8, 0.83, 0.84], [0.5, 0.6, 0.72, 0.8],
                                         [0.45, 0.55, 0.7, 0.79]]),
              ('(b) slowdown 1/Π', [[0.1, 0.1, 0.1, 0.1], [0.25, 0.22, 0.16, 0.12], [0.75, 0.6, 0.35, 0.18],
                                    [0.82, 0.68, 0.42, 0.22]]),
              ('(c) empty sites θ∗', [[0.85, 0.9, 0.93, 0.95], [0.7, 0.8, 0.88, 0.92], [0.3, 0.5, 0.75, 0.88],
                                      [0.25, 0.45, 0.7, 0.85]])]
    for p, (title, ys) in enumerate(panels):
        X, Y, Wd, Ht = 0.95 + p * 4.15, 2.15, 3.4, 2.3
        text(sl, X - 0.1, Y - 0.45, Wd + 0.2, 0.32, title, size=14, bold=True)
        axes(sl, X, Y, Wd, Ht, '', '')
        if p == 1:
            line(sl, X, Y + Ht * 0.9, X + Wd, Y + Ht * 0.9, color='muted', lw=1.0, dash=True)
        for i, case in enumerate(cases):
            text(sl, X + 0.15 + i * 0.82 - 0.3, Y + Ht + 0.05, 0.9, 0.28, case, size=12, align='c', color='muted')
            for j, (kind, c, _) in enumerate(marks):
                cx = X + 0.45 + i * 0.82 + (j - 1.5) * 0.12
                cy = Y + Ht * (1 - ys[i][j])
                shape(sl, kind, cx - 0.07, cy - 0.07, 0.14, 0.14, fill=c, line=None)
                if p == 1 and i > 0 and j == 0:                         # open marker (650 K only): fixed-site estimate
                    shape(sl, kind, cx - 0.08, cy - 0.28, 0.16, 0.16, fill=None, line=c, lw=1.0)
    for j, (kind, c, lab) in enumerate(marks):
        shape(sl, kind, 0.95 + j * 1.25, 5.0, 0.16, 0.16, fill=c, line=None)
        text(sl, 1.17 + j * 1.25, 4.91, 1.0, 0.3, lab, size=12)
    text(sl, 6.1, 4.88, 6.7, 0.95, 'Open circles in (b): the slowdown if the number of sites were fixed. _vs_P: one row '
         'per pressure-series temperature, a marker per p. _per_area copies show k_{d,eff} in (a).', size=13,
         color='muted')
    schematic(sl, 11.35, 1.62)
    actions(sl, ['Plot each case’s forward TOF on a log axis.',
                 'Show how many times slower than H₂ alone.',
                 'Show the empty sites that the rate needs.',
                 'Repeat per series, and per area.'])
    notes(sl, ['analyses/steady_state.py plot_rates() (lines 63-131): one row of panels per group of the series (the',
               'temperature series: one row, a marker per T; each pressure series: a row per pressure-series temperature,',
               'a marker per p), at most three rows per figure. slowdown_fixed_sites = slowdown_x (f/f_ref)^n (line 79).',
               'per_area=True: panel (a) is kd_eff = J_diss / p (Format 2). Files: rates_summary_vs_T.png, _vs_P.png and',
               'their _per_area copies. Docs: Stage 9 (line 807) and Table R2.',
               'The dots are a layout sketch, not results.'])


SLIDES = [s39_forward_net, s40_poisoning, s41_table, s42_composition, s43_rate_plots]
