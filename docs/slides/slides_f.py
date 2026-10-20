"""Section F (slides 26-27): Stage 3, designing the gas conditions: the fixed-probe basis with
per-impurity ratios from RMG reactor 1, and the two series. Running example: runs/try_2/Cu."""
import math

from deck_lib import *  # noqa: F401,F403

PF = 'F · Designing the conditions · Stage 3'
RUN = 'runs/try_2/Cu'
# p(impurity) / p(H2) of RMG reactor 1 (H2 = 190), as in settings.py impurity_ratio
RATIOS = [('H₂O', 'H2O', 10 / 190), ('CO', 'CO', 1 / 190), ('N₂', 'N2', 0.85 / 190), ('CH₄', 'CH4', 0.8 / 190),
          ('O₂', 'O2', 0.15 / 190), ('NH₃', 'NH3', 0.1 / 190), ('CO₂', 'CO2', 0.05 / 190)]


def s26_basis(d):
    sl = d.slide('Stage 3: hold H₂ fixed and add each impurity on top, at RMG reactor 1’s ratios',
                 'Every case has the same p(H₂), so a change in rate comes from the impurity, not from less H₂.',
                 part=PF, tag='package · conditions.py')
    # left: the nine cases at p(H2) = 1 bar
    card(sl, 0.55, 1.66, 5.95, 4.22, 'Every case at p(H₂) = 1 bar')
    rows = [('H₂ only', None, 0.0)] + [(lab, key, r) for lab, key, r in RATIOS] + \
           [('all', 'dark', sum(r for _, _, r in RATIOS))]
    for i, (lab, key, r) in enumerate(rows):
        y = 2.08 + i * 0.4
        chip(sl, 0.75, y, lab, key=key, w=1.05, h=0.32, size=12)
        shape(sl, 'rect', 1.95, y + 0.02, 2.3, 0.28, fill='H', line=None, content='H₂  1 bar', size=12,
              color='white', bold=True)
        if r:
            txt = f'+ {r:.4f} bar' if lab != 'all' else f'+ {r:.4f} bar, all seven'
            text(sl, 4.35, y, 2.1, 0.32, txt, size=13, color=key if key != 'dark' else 'ink', bold=True, anchor='m')
    # right: the ratios on a log scale
    card(sl, 6.83, 1.66, 5.95, 4.22, 'p(impurity) / p(H₂), in ppm (log scale)')
    x0, x1, lo, hi = 8.05, 11.65, 2.0, 5.0                     # 100 to 100,000 ppm
    X = lambda ppm: x0 + (math.log10(ppm) - lo) / (hi - lo) * (x1 - x0)          # noqa: E731
    y_axis = 5.18
    for k in range(3, 6):
        xx = X(10 ** k / 10)
        line(sl, xx, 2.12, xx, y_axis, color='rule', lw=0.75)
    for k, lab in ((2, '100'), (3, '1,000'), (4, '10,000'), (5, '100,000')):
        text(sl, X(10 ** k) - 0.6, y_axis + 0.04, 1.2, 0.28, lab, size=12, color='muted', align='c')
    line(sl, x0, y_axis, x1, y_axis, color='muted', lw=1.0)
    for i, (lab, key, r) in enumerate(RATIOS):
        y = 2.2 + i * 0.42
        ppm = r * 1e6
        chip(sl, 7.0, y - 0.15, lab, key=key, w=0.85, h=0.3, size=12)
        line(sl, x0, y, X(ppm), y, color='rule', lw=0.75, dash=True)
        circle(sl, X(ppm), y, 0.2, fill=key, line=None)
        text(sl, X(ppm) + 0.14, y - 0.16, 1.1, 0.3, f'{ppm:,.0f}', size=12, color='ink')
    actions(sl, ['Hold p(H₂) at its series value in every case.',
                 'Add each impurity at its ratio to H₂, from RMG reactor 1.',
                 '“all” adds all seven: reactor 1’s own gas.',
                 'A change in rate then comes from the impurity, not from less H₂.'])
    notes(sl, [f'{RUN}/settings.py lines 49-54: impurity_ratio = {{H2O: 10/190, CO2: 0.05/190, CO: 1/190, CH4: 0.8/190,',
               "N2: 0.85/190, O2: 0.15/190, NH3: 0.1/190}}, mixture_mode 'same_each', include_input_basis False.",
               f'Ratios from RMG reactor 1 of {RUN}/input.py (lines 108-132: H2 190, H2O 10, CO 1, N2 0.85, CH4 0.8, O2 0.15,',
               'NH3 0.1, CO2 0.05). ppm = ratio x 1e6. The same ratios hold at every point of both series.',
               'conditions.py partial_pressures (lines 76-99): fixed_probe basis, p(H2) = P_ref and each impurity at',
               "ratio_for(r, label) x P_ref (ratio_for, lines 65-68); 'split' would divide one ratio over the mixture.",
               'check output (run.log): P_total 1.052632 bar for H2O ... 1.068158 bar for "all" at p(H2) = 1 bar.',
               'Docs Stage 3 (simulation_workflow.md line 488). The input_X basis (the RMG mole fractions, diluting H2)',
               'is a check, off here.'])


def s27_series(d):
    sl = d.slide('Two series, never mixed: temperature at 1 bar, pressure at 650 and 1200 K',
                 'Temperature and pressure effects stay apart. Both stay inside the range RMG built the mechanism for.',
                 part=PF, tag='package · conditions.py')
    # the map: T (K) across, log p(H2) up
    PX, PY, PW, PH = 1.35, 1.95, 5.6, 3.45
    t_lo, t_hi, p_lo, p_hi = 600.0, 1250.0, math.log10(0.004), math.log10(300.0)
    Xt = lambda T: PX + (T - t_lo) / (t_hi - t_lo) * PW                                # noqa: E731
    Yp = lambda p: PY + PH - (math.log10(p) - p_lo) / (p_hi - p_lo) * PH                # noqa: E731
    shape(sl, 'rect', Xt(650), Yp(100), Xt(1200) - Xt(650), Yp(0.01) - Yp(100), fill='panel', line='rule', lw=0.75)
    text(sl, Xt(650) + 0.3, Yp(100) + 0.04, 3.6, 0.28, 'RMG’s reactors: 650–1200 K, 0.01–100 bar', size=12,
         color='muted', italic=True)
    axes(sl, PX, PY, PW, PH, '', '')                                # labels placed below, clear of the ticks
    text(sl, PX, PY + PH + 0.34, PW, 0.3, 'T (K)', size=14, color='muted', align='c')
    text(sl, PX - 0.75, PY - 0.42, 2.6, 0.3, 'p(H₂), bar (log scale)', size=14, color='muted')
    for T in (650, 800, 1000, 1200):
        text(sl, Xt(T) - 0.4, PY + PH + 0.06, 0.8, 0.26, str(T), size=12, color='muted', align='c')
    for p, lab in ((0.01, '0.01'), (1, '1'), (100, '100')):
        text(sl, PX - 0.62, Yp(p) - 0.14, 0.55, 0.28, lab, size=12, color='muted', align='r')
    for T in (650, 800, 1000, 1200):                               # temperature series: filled circles at 1 bar
        circle(sl, Xt(T), Yp(1), 0.22, fill='dark', line=None)
    for T in (650, 1200):                                          # pressure series: open squares
        for p in (0.01, 1, 100):
            shape(sl, 'rect', Xt(T) - 0.17, Yp(p) - 0.17, 0.34, 0.34, fill=None, line='dark', lw=1.5)
    for T in (650, 1200):
        text(sl, Xt(T) - 0.45, Yp(1) - 0.5, 0.9, 0.28, 'both', size=12, color='ink', bold=True, align='c')
    # legend and the count
    circle(sl, 7.65, 2.1, 0.22, fill='dark', line=None)
    text(sl, 7.9, 1.95, 4.8, 0.3, 'temperature series: 650, 800, 1000, 1200 K at 1 bar', size=14)
    shape(sl, 'rect', 7.48, 2.43, 0.34, 0.34, fill=None, line='dark', lw=1.5)
    text(sl, 7.9, 2.45, 4.8, 0.3, 'pressure series: 0.01, 1, 100 bar at 650 and 1200 K', size=14)
    card(sl, 7.4, 3.0, 5.38, 2.88, None)
    text(sl, 7.6, 3.12, 5.0, 0.4, '4 + 2 × 3 − 2 shared = **8 points**', size=17)
    text(sl, 7.6, 3.58, 5.0, 0.4, '× 9 cases = **72 main conditions**', size=17)
    text(sl, 7.6, 4.08, 5.0, 0.3, 'each with a key, for example', size=13, color='muted')
    text(sl, 7.6, 4.36, 5.1, 0.3, '`fixed_probe | CO | 650 K | 1 bar`', size=13)
    text(sl, 7.6, 4.78, 5.0, 1.0, 'Every point-based analysis runs on both series; its own section can narrow that '
         '(series, temperatures, pressures, pressure_temperatures).', size=13, color='muted')
    actions(sl, ['Temperature series: vary T at p(H₂) = 1 bar.',
                 'Pressure series: vary p(H₂) at 650 K and at 1200 K.',
                 'A point in both is solved once and labelled “both”.',
                 '8 points × 9 cases = 72 main conditions.'])
    notes(sl, [f'{RUN}/settings.py lines 45-48: temperatures [650, 800, 1000, 1200] K at reference_pressure 1e5 Pa;',
               'pressures [1e3, 1e5, 1e7] Pa (0.01, 1, 100 bar) at pressure_temperatures [650, 1200] K.',
               'conditions.py series_points (lines 26-42): the temperature series first, then each pressure series; a',
               'point in both gets series "both". main_conditions (120-128) makes every case at every point;',
               'condition_key (137-139) writes "basis | case | T | p". The labels P_ref_bar and series are not part of',
               'the cache fingerprint (split_labels, 131-134), so relabelling never re-solves.',
               f'RMG range: {RUN}/input.py, temperature=[(650,\'K\'),(1200,\'K\')] in every reactor, initialPressure 1, 100',
               'and 0.01 bar. The 72 conditions are listed by check (run.log, 2026-10-15).',
               'Docs Stage 3 (line 488) and "Both series by default" (Stage 11 conventions).'])


SLIDES = [s26_basis, s27_series]
