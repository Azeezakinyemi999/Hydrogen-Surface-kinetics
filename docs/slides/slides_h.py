"""Section H (slides 36-38): Stages 5-6, two starting surfaces per condition, solve_condition(),
and the parallel, cached runner. Running example: runs/try_2/Cu."""
from deck_lib import *  # noqa: F401,F403

PH = 'H · Two starts, many conditions · Stages 5–6'
RUN = 'runs/try_2/Cu'


def s36_two_starts(d):
    sl = d.slide('Stage 5: two starting surfaces, because one gas can have two steady states',
                 'Start A is a clean surface. Start B already met the gas, or the impurities alone.',
                 part=PH, tag='concept')
    occA = {}
    occB = {(r, c): ('CO', 'CO', 1) for r in range(3) for c in range(4) if (r + c) % 3 != 0}
    occS = {(0, 0): ('H', 'H', 1), (0, 2): ('H', 'H', 1), (1, 1): ('H', 'H', 1), (1, 3): ('CO', 'CO', 1),
            (2, 0): ('H', 'H', 1), (2, 2): ('H', 'H', 1), (0, 3): ('H', 'H', 1)}
    for row, (lab, sub, occ, y) in enumerate([('Start A', 'clean surface', occA, 1.95),
                                              ('Start B', 'met CO first', occB, 3.85)]):
        text(sl, 0.55, y - 0.1, 1.6, 0.35, lab, size=16, bold=True)
        text(sl, 0.55, y + 0.25, 1.6, 0.6, sub, size=14, color='muted')
        grid(sl, 2.25, y - 0.15, 4, 3, pitch=0.4, d=0.3, occ=occ)
        arrow(sl, 4.25, y + 0.42, 5.15, y + 0.42, color='ink', lw=2)
        text(sl, 4.1, y + 0.02, 1.3, 0.3, 'solve', size=13, color='muted', align='c')
        grid(sl, 5.35, y - 0.15, 4, 3, pitch=0.4, d=0.3, occ=occS)
        text(sl, 5.2, y + 0.98, 2.0, 0.3, f'steady state {lab[-1]}', size=13, color='muted', align='c')
        arrow(sl, 7.35, y + 0.42, 7.9, 3.32 if row == 0 else 3.52, color='ink', lw=1.5)
    gate(sl, 7.9, 2.75, 2.6, 1.35, 'max |Δθ| ≤ 10^{−3} and TOF within 1\u00a0%?', size=12)
    good(sl, 10.75, 2.0, 2.03, 0.95, 'agree: one steady state found', size=13)
    caveat(sl, 10.75, 3.6, 2.03, 1.55, 'differ: two steady states; both are reported', size=13)
    arrow(sl, 10.5, 3.2, 10.75, 2.47, color='ok')
    arrow(sl, 10.5, 3.6, 10.75, 4.2, color='amber')
    text(sl, 7.7, 4.3, 3.0, 1.55, 'Start B, H₂ only: 99\u00a0% H∗, 1\u00a0% empty. With impurities: a clean surface held 1\u00a0s '
         'in the impurities alone, H₂ removed (e.g. 0.0526 bar H₂O).', size=13, color='muted')
    actions(sl, ['Start A: a clean surface, as for a fresh catalyst.',
                 'Start B: a surface already covered, by H∗ or by the impurities.',
                 'Solve both to steady state.',
                 'If they differ, this gas has more than one steady state.'])
    notes(sl, ['docs/simulation_workflow.md Stage 5, Figure 5 (fig_bistability.png) and Diagram 4.',
               'Agreement: max |θ_A − θ_B| ≤ 1e-3 and |TOF_A − TOF_B| / max ≤ 1e-2 (config.py:30-31, multi_ss_*).',
               'Start B: probe_covered_fraction 0.99 (config.py:28); pre_exposure_time 1 s (config.py:29); for try_2/Cu at',
               '1 bar the H2O case holds the surface in 0.0526 bar of H2O alone (impurity_ratio 10/190).',
               'Both are kept: start B\'s values go into the _B columns (tof_fwd_B, theta_free_B, poisoning_factor_B) and the',
               'poisoning curves draw its branch dashed (Part K).',
               'Agreement is evidence, not proof, of a unique steady state. One failed start: agreement unknown (None).',
               'The grids are schematic (the "steady state" grids are not computed results).'])


def s37_solve_condition(d):
    sl = d.slide('solve_condition(): one condition, two starts, one verdict',
                 'Start B is built to suit the case: H∗-covered for H₂ only, impurity-exposed otherwise.',
                 part=PH, tag='package · solver.py:233–254')
    step(sl, 0.55, 1.75, 2.25, 0.8, ['row → T, P, X', 'bar → Pa'], where='line 236')
    arrow(sl, 2.8, 2.15, 3.1, 2.15)
    step(sl, 3.1, 1.75, 2.35, 0.8, ['solve from start A', 'clean surface'], where='lines 237–238')
    arrow(sl, 5.45, 2.15, 5.75, 2.15)
    gate(sl, 5.75, 1.62, 1.95, 1.05, 'both starts? (240)', size=12)
    arrow(sl, 7.7, 2.15, 8.05, 2.15, color='muted')
    yes_no(sl, 7.72, 1.84, 'no', w=0.4)
    step(sl, 8.05, 1.82, 1.85, 0.65, 'return start A only', size=13)
    arrow(sl, 6.72, 2.67, 6.72, 2.97)
    yes_no(sl, 6.8, 2.66, 'yes', w=0.45)
    gate(sl, 5.75, 2.97, 1.95, 1.05, 'H₂ only in the gas? (242)', size=12)
    arrow(sl, 7.7, 3.5, 8.05, 3.5)
    yes_no(sl, 7.72, 3.18, 'yes', w=0.45)
    step(sl, 8.05, 3.15, 2.1, 0.7, ['start B:', '99 % H∗, 1 % empty'], size=13, where='lines 186–188')
    arrow(sl, 5.75, 3.5, 5.45, 3.5)
    yes_no(sl, 5.42, 3.18, 'no', w=0.4)
    step(sl, 0.55, 3.0, 4.9, 1.0, ['`impurity_covered_start()`',
                                    'impurities alone, at their own pressures · 1 s in the fixed-gas reactor · '
                                    'clip at 0, renormalise'], size=13, where='lines 194–218')
    poly(sl, [(3.0, 4.0), (3.0, 4.75), (8.05, 4.75)], end=True)
    arrow(sl, 8.55, 3.85, 8.55, 4.42)
    step(sl, 8.05, 4.42, 2.1, 0.65, 'solve from start B', size=13, where='lines 250–251')
    arrow(sl, 10.15, 4.75, 10.45, 4.75)
    step(sl, 10.45, 4.32, 2.33, 0.85, ['`compare_starts()`', 'agree · differ · None'], size=13,
         where='lines 221–230')
    text(sl, 0.55, 4.88, 2.4, 0.95, 'pre-exposure fails: start B = None, agreement unknown (246–248)', size=12,
         color='muted')
    text(sl, 3.15, 4.85, 4.8, 0.6, 'start B cached per mechanism, T and impurity pressures (202–204)',
         size=12, color='muted')
    actions(sl, ['Solve the condition from the clean start.',
                 'Build start B to suit the case.',
                 'Solve again from start B.',
                 'Compare: agree, differ, or unknown.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py (lines unchanged since the first commit): solve_condition() 233-254,',
               'probe_covered_start() 186-188 (probe_covered_fraction 0.99), impurity_covered_start() 194-218',
               '(pre_exposure_time 1 s, probe gas removed, impurity partial pressures unchanged; cache key at line 202),',
               'compare_starts() 221-230. Doc Stage 5 design choice: pre-exposed rather than "steady state in the',
               'impurity alone", because some impurity-only systems never settle.'])


def s38_runner(d):
    sl = d.slide('Stage 6: run every condition in parallel, each saved the moment it finishes',
                 'A fingerprint of everything that affects the answer names each saved result.',
                 part=PH, tag='package · runner.py:61–123')
    keys = ['fixed_probe | H2O | 800 K | 1 bar', 'fixed_probe | CO | 650 K | 1 bar',
            'fixed_probe | N2 | 1200 K | 1 bar', 'fixed_probe | all | 650 K | 1 bar']
    for i, k in enumerate(keys):
        shape(sl, 'rrect', 0.55, 1.75 + i * 0.5, 3.5, 0.4, fill='white', line='faint', radius=0.06, content=k,
              size=12, font=CODE, wrap=False)
    text(sl, 0.55, 3.75, 3.5, 0.3, 'jobs: one per condition (72 here)', size=13, color='muted')
    arrow(sl, 4.05, 2.65, 4.3, 2.65)
    gate(sl, 4.3, 2.05, 2.2, 1.2, 'saved with this fingerprint?', size=12, where='lines 75–79')
    step(sl, 4.57, 1.62, 1.65, 0.36, 'yes: load it', size=12)
    arrow(sl, 5.4, 2.05, 5.4, 1.98, color='ok')
    arrow(sl, 6.5, 2.65, 6.8, 2.65)
    yes_no(sl, 6.5, 2.32, 'no', w=0.4)
    for i in range(3):
        y = 1.85 + i * 0.62
        shape(sl, 'rrect', 6.8, y, 2.35, 0.48, fill='white', line='dark', lw=1.0, radius=0.06,
              content=f'worker process {i + 1}', size=12)
    text(sl, 6.8, 3.75, 2.4, 0.55, 'forked; half the logical CPUs: 8 here (line 58)', size=12, color='muted')
    arrow(sl, 9.15, 2.65, 9.4, 2.65)
    outs = [('ok', 'done: save .pkl at once', 'lines 113–115'), ('dark', 'worker died: error', 'lines 102–103'),
            ('amber', 'over 1800 s: stop, timeout', 'lines 104–107')]
    for i, (c, t, r) in enumerate(outs):
        y = 1.75 + i * 0.72
        shape(sl, 'rrect', 9.4, y, 3.38, 0.5, fill='white', line=c, lw=1.5, radius=0.08, content=t, size=13)
        ref(sl, 9.4, y + 0.49, 3.38, r, align='r')
    text(sl, 9.4, 3.95, 3.38, 0.3, 'main loop checks every 0.2 s (line 121)', size=12, color='muted', align='r')
    card(sl, 0.55, 4.35, 7.6, 1.55, 'The fingerprint (file name)')
    text(sl, 0.7, 4.75, 7.35, 0.6, '<key>_<first 12 hex of MD5 over the job + mechanism SHA-256, phase names, '
         'solver settings, Cantera version, both starts>.pkl', size=13)
    text(sl, 0.7, 5.4, 7.35, 0.35, 'e.g.  `fixed_probe_CO_650_K_1_bar_9d1585270cd4.pkl`', size=13, color='muted')
    caveat(sl, 8.35, 4.35, 4.43, 1.55, 'Time limits count wall-clock time, so a sleeping Mac gives false timeouts. '
           'run.sh runs under caffeinate \u2011dims: keep the lid open.', size=13)          # non-breaking hyphen
    actions(sl, ['Load any condition whose fingerprint is already saved.',
                 'Solve the rest in parallel worker processes.',
                 'Stop a worker that runs past the time limit.',
                 'Save each finished result at once.'])
    notes(sl, ['Hydrogen-Surface-kinetics runner.py run_jobs() lines 61-123: cache check 74-79 (a loaded result gets this',
               "run's job at line 78, so relabelled points need no new solve), fork context 71, workers started 86-91,",
               'error 102-103, timeout 104-107, saved 113-115, poll every 0.2 s at 121.',
               'cache_path() 42-45: MD5 of the job without its key and without fields starting with "_" (the labels',
               'P_ref_bar and series, line 43) plus the context; safe_name() 38-39 turns the key into the file name.',
               'Context: solver.py line 287 passes (model.fingerprint_context(), both_starts); model.py 220-223: mechanism',
               'SHA-256, phase names, solver settings, Cantera version. Defaults (config.py:34-38): workers = half the logical',
               'CPUs (16 on this Mac, so 8), job_timeout 1800 s, recompute False; run.sh options --workers, --timeout, --recompute.',
               f'Real file: {RUN}/output/cache/fixed_probe_CO_650_K_1_bar_9d1585270cd4.pkl.',
               'run.sh lines 54-55 (caffeinate -dims); doc Stage 6 and "Long runs and sleep" (line 1138): -i is not enough,',
               'because it does not stop maintenance sleep.'])


SLIDES = [s36_two_starts, s37_solve_condition, s38_runner]
