"""Section E (slides 24-25): Stage 2, replaying RMG's own reactors in Cantera and reading the
flux diagrams. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PE = "E · Replaying RMG's reactors · Stage 2"
RUN = 'runs/try_2/Cu'


def s24_stage2(d):
    sl = d.slide("Stage 2 replays RMG's reactors in Cantera and draws where flux goes",
                 'A closed box like RMG’s, integrated from 10⁻⁸ to 1000 s, with one flux diagram per moment.',
                 part=PE, tag='package · analyses/flux.py')
    ins = [('input.py', "15 reactors, read by RMG's parser"), ('seed/seed/', "RMG's species and reactions"),
           ('chem_annotated.yaml', 'the model Cantera runs')]
    for i, (a, b) in enumerate(ins):
        step(sl, 0.55, 1.75 + i * 1.08, 2.7, 0.88, [('`' + a + '`', {'bold': True}), b], size=14)
        arrow(sl, 3.25, 2.19 + i * 1.08, 3.6, 3.27)
    chain = [['match species', 'by label, then by structure'], ['closed batch reactor', '10⁻⁸ to 1000 s', '100 time points'],
             ['one flux diagram', 'per time point'], ['one video per run', 'flux/<n>_<T>K/']]
    for i, t in enumerate(chain):
        x = 3.62 + i * 2.33
        step(sl, x, 2.8, 2.05, 0.95, t, size=14, fill='white' if i < 3 else 'panel2')
        if i < 3:
            arrow(sl, x + 2.05, 3.27, x + 2.33, 3.27)
    text(sl, 3.62, 1.72, 9.1, 0.95, 'Optional. Run with `bash run.sh runs/try_2/Cu --only flux`. Needs RMG-Py '
         "(flux['rmg_path'] in settings.py), RDKit, Graphviz and ffmpeg; one reactor at a time.", size=14,
         color='muted')
    text(sl, 3.62, 3.95, 9.1, 0.65, 'The reactor is sized like RMG’s: V = RT/P for 1 mol of gas, and catalyst area '
         'A = V × 10⁵ m⁻¹ from input.py.', size=14, color='muted')
    caveat(sl, 3.62, 4.78, 9.16, 1.0, 'Not run yet for any try_* mechanism. Every try_2/Cu reactor spans 650–1200 K, '
           'so it has no single T (RMG samples the range): each runs at both ends, 650 and 1200 K, 30 runs in all, '
           'or at flux[\'temperatures\'].', size=14)
    actions(sl, ['Read the 15 reactors with RMG’s own input parser.',
                 'Match RMG’s species and reactions to Cantera’s.',
                 'Run each reactor at 650 and 1200 K, 10⁻⁸ to 1000 s.',
                 'Draw a flux diagram at each time and stitch a video.'])
    notes(sl, ['docs/simulation_workflow.md Stage 2 (line 412) and Diagram 2. Code: analyses/flux.py: load_seed (79),',
               'mol_from_pt_smiles (88), build_matched_model (111), reaction_index (160), draw_species_images (172),',
               'reactor_temperatures (187), simulate (198-222: IdealGasReactor energy off, V = RT/P, A = V x',
               'surfaceVolumeRatio, logspace from 1e-8 s to terminationTime), condition_label (225),',
               'bake_condition_label (230), run (249-325), and RMG\'s rmgpy.tools.fluxdiagram.generate_flux_diagram.',
               'Settings: config.py flux section (lines 144-155: rmg_path, rmg_input, seed_dir, species_dir,',
               'reactor_indices, temperatures, n_points 100, rtol 1e-10, atol 1e-22, skip_existing). One reactor at a',
               'time in the main process: forking after RDKit/Cairo drawing crashes on macOS (flux.py:309-310; docs',
               'line 479). Species matching: base label first; otherwise the Cantera name is read as SMILES with [Pt]',
               'turned back into the site atom X, and matched by isomorphism.',
               f'Temperatures: {RUN}/input.py gives every surfaceReactor temperature=[(650,\'K\'),(1200,\'K\')].',
               'RMG-Py rmgpy/solver/surface.pyx lines 101-104 then set Trange and leave T unset; RMG itself samples',
               'the range (rmgpy/rmg/main.py RMG_Memory.generate_cond, line 2513: the first condition is the low end,',
               'later ones spread over the range). flux.reactor_temperatures runs such a reactor at both ends by',
               "default, or at flux['temperatures']; folders flux/<n>_<T>K/ (a one-temperature reactor keeps",
               'flux/<n>/). 15 reactors x 2 = 30 runs. No try_* run has a flux/ output yet. Docs Stage 2 step 6 (line 470).'])


def s25_flux_reading(d):
    sl = d.slide('Reading a flux diagram: arrow width is flux, on a log scale',
                 'Arrows join species that swap atoms. The thinnest arrow drawn is 10⁻⁶ of the largest.',
                 part=PE, tag='concept')
    card(sl, 0.55, 1.66, 6.3, 4.22, None)
    nodes = {'H₂': (1.55, 2.4, 'H'), 'H∗': (5.3, 2.4, 'H'), 'H₂∗': (1.55, 4.55, 'H'), 'CO': (3.6, 4.55, 'CO'),
             'CO∗': (5.6, 4.55, 'CO')}
    for t, (cx, cy, k) in nodes.items():
        shape(sl, 'rrect', cx - 0.5, cy - 0.27, 1.0, 0.54, fill='white', line=k, lw=2.0, radius=0.1, content=t,
              size=16, color=k, bold=True)
    arrow(sl, 2.08, 2.4, 4.77, 2.4, color='H', lw=7, head='sm')
    arrow(sl, 1.55, 2.69, 1.55, 4.26, color='H', lw=3, head='sm')
    arrow(sl, 4.12, 4.55, 5.08, 4.55, color='CO', lw=9, head='sm')
    text(sl, 2.6, 1.95, 2.0, 0.3, 'large flux', size=13, color='muted', align='c')
    text(sl, 1.65, 3.3, 1.6, 0.3, 'smaller', size=13, color='muted')
    schematic(sl, 0.7, 5.5, w=3.5, label='schematic, not a computed diagram')
    card(sl, 7.15, 1.66, 5.63, 2.6, 'Width rule')
    text(sl, 7.3, 2.1, 5.4, 0.45, 'w = w_{max} (1 − log₁₀ x / log₁₀ x_{tol})', size=18, font=MATH)
    text(sl, 7.3, 2.62, 5.4, 0.4, 'x = |R_{ab}| / max |R|,   x_{tol} = 10^{−6},   w_{max} = 10', size=16, font=MATH)
    for i, (xl, lw_, t) in enumerate([('x = 1', 8, 'widest'), ('x = 10^{−3}', 4, 'half width'),
                                      ('x < 10^{−6}', 0, 'not drawn')]):
        y = 3.2 + i * 0.33
        text(sl, 7.3, y - 0.14, 1.5, 0.3, xl, size=14, font=MATH)
        if lw_:
            line(sl, 8.9, y + 0.02, 10.6, y + 0.02, color='dark', lw=lw_)
        else:
            line(sl, 8.9, y + 0.02, 10.6, y + 0.02, color='faint', lw=1, dash=True)
        text(sl, 10.75, y - 0.14, 1.9, 0.3, t, size=14, color='muted')
    caveat(sl, 7.15, 4.45, 5.63, 1.43, 'Node borders put gas (kmol/m³) and surface (kmol/m²) concentrations on '
           'one scale. Compare borders only within one phase; arrow widths compare everywhere.')
    actions(sl, ['Each arrow joins two species that exchange atoms.',
                 'Its width grows with the log of the flux.',
                 'Below 10⁻⁶ of the largest flux, no arrow is drawn.',
                 'Compare node borders only within one phase.'])
    notes(sl, ['docs/simulation_workflow.md Stage 2, "How the diagram encodes it" (line 450, width rule line 454) and the',
               'Important box on node borders (line 460); Figure 3 (fig_flux_anatomy.png).',
               'RMG-Py rmgpy/tools/fluxdiagram.py (local checkout): max_node_count = max_edge_count = 50 (lines 55-56),',
               'concentration_tol = species_rate_tol = 1e-6 (57-58), max_node_pen_width = max_edge_pen_width = 10 (59-60);',
               'edge width = slope log10(x) + w_max with slope = -w_max / log10(tol) (line 250).',
               'The steady-state path diagrams of Stage 11 (Part K) use their own style and threshold.',
               'The diagram on this slide is schematic.'])


SLIDES = [s24_stage2, s25_flux_reading]
