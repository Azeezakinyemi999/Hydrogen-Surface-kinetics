"""Section K (slides 48-58): Stage 11, the further analyses: poisoning curves, site blocking, the
temperature and pressure sweeps, the Arrhenius form, equilibrium coverage, reaction orders, kinetic
and thermodynamic rate control, reaction paths and poisoning over time.
Running example: runs/try_2/Cu. The counts are from its check output; every figure is a layout sketch."""
import math

from deck_lib import *  # noqa: F401,F403

PK = 'K · Further analyses · Stage 11'
RUN = 'runs/try_2/Cu'
LEVELS = [1e-6, 3e-6, 1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.25]       # config.py curves ratios


def sig(t, mid, width, lo=0.0, hi=1.0):
    """A smooth step from hi (t well below mid) to lo (t well above mid): sketch curves only."""
    return lo + (hi - lo) / (1.0 + math.exp((t - mid) / width))


def at(f, t, t0, t1, X, Y, Wd, Ht):
    """Slide position of the point (t, f(t)) of a sketch curve drawn by curve()."""
    return X + Wd * (t - t0) / (t1 - t0), Y + Ht * (1 - f(t))


def panel(sl, X, Y, Wd, Ht, title, xlab=None, size=13):
    """A small plot: bold title above, L-shaped axes, an optional x label below."""
    text(sl, X - 0.05, Y - 0.42, Wd + 0.4, 0.32, title, size=size, bold=True)
    axes(sl, X, Y, Wd, Ht, '', '')
    if xlab:
        text(sl, X, Y + Ht + 0.04, Wd, 0.28, xlab, size=12, color='muted', align='c')


def half_mark(sl, cx, cy, color, d=0.2):
    """The open diamond the package draws on the 0.5 line (half-rate level, release T, ...)."""
    shape(sl, 'diamond', cx - d / 2, cy - d / 2, d, d, fill='white', line=color, lw=1.75)


def s48_curves(d):
    sl = d.slide('Poisoning curves: the rate against the impurity level',
                 'At each point, every impurity at 12 levels, each from both starts: 680 solves for try_2/Cu.',
                 part=PK, tag='package · curves.py:22–71')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'What is solved, and what is computed')
    x0, x1, ya = 0.95, 6.05, 2.58
    lo, hi = -6.0, math.log10(0.25)
    X = lambda r: x0 + (math.log10(r) - lo) / (hi - lo) * (x1 - x0)           # noqa: E731
    text(sl, 0.75, 2.06, 5.6, 0.3, 'the 12 levels r = p(impurity) / p(H₂), on a log scale', size=13, color='muted')
    line(sl, x0, ya, x1, ya, color='muted', lw=1.0)
    for r in LEVELS:
        circle(sl, X(r), ya, 0.15, fill='dark', line=None)
    for r, lab in ((1e-6, '10^{−6}'), (1e-4, '10^{−4}'), (1e-2, '10^{−2}'), (0.25, '0.25')):
        text(sl, X(r) - 0.4, ya + 0.1, 0.8, 0.3, lab, size=12, color='muted', align='c')
    text(sl, 0.75, 3.08, 5.6, 0.62, '8 points × (1 reference + 7 impurities × 12 levels) = **680 solves**, '
         'from both starts', size=14)
    for y, h, t in ((3.78, 0.36, 'Π(r) = TOF_{f}(r) / TOF_{f}(H₂ only), at the same point'),
                    (4.16, 0.62, 'θ_{imp} = 1 − θ_{∗} − θ_{H∗} − θ_{H₂∗}: the sites the impurity’s species hold'),
                    (4.62, 0.95, 'r_{½}: where Π = 0.5, by log–log interpolation; “at most this” if Π is already '
                                'below 0.5 at 10^{−6}, “never halved in the range tested” if it never gets there')):
        text(sl, 0.75, y, 5.6, h, t, size=14)
    # right: the figure, one column per point (rate on top, impurity-held sites below)
    card(sl, 6.83, 1.66, 5.95, 4.22, 'The figure: one column per point')
    schematic(sl, 11.55, 1.78, w=1.15, label='layout only')
    PX, PW, t0, t1 = 7.65, 4.6, -6.0, -0.6
    tx = lambda t: PX + (t - t0) / (t1 - t0) * PW                               # noqa: E731
    win = (-4.75, -3.55)                                                        # where the two starts disagree
    Y1, H1 = 2.3, 1.45
    shape(sl, 'rect', tx(win[0]), Y1, tx(win[1]) - tx(win[0]), H1, fill='panel2', line=None)
    axes(sl, PX, Y1, PW, H1, '', '')
    text(sl, 6.95, Y1 + H1 / 2 - 0.18, 0.6, 0.36, 'Π', size=15, color='muted', align='r')
    line(sl, PX, Y1 + H1 / 2, PX + PW, Y1 + H1 / 2, color='muted', lw=1.0, dash=True)     # Π = 0.5
    fA = lambda t: sig(t, -3.0, 0.42, 0.06, 0.9)                                # noqa: E731
    fB = lambda t: sig(t, -4.4, 0.3, 0.06, 0.9)                                 # noqa: E731
    curve(sl, fA, t0, t1, PX, Y1, PW, H1, color='CO', lw=2.0)
    wx, ww = tx(win[0]), tx(win[1]) - tx(win[0])                               # start B: only in the window
    curve(sl, fB, win[0], win[1], wx, Y1, ww, H1, color='CO', lw=1.5, dash=True)
    two = [t for t in map(math.log10, LEVELS) if win[0] < t < win[1]]
    for t in map(math.log10, LEVELS):
        circle(sl, *at(fA, t, t0, t1, PX, Y1, PW, H1), 0.12, fill='CO', line=None)
    for t in two:
        circle(sl, *at(fA, t, t0, t1, PX, Y1, PW, H1), 0.26, fill=None, line='ink', lw=1.25)
        bx, by = at(fB, t, t0, t1, PX, Y1, PW, H1)
        shape(sl, 'rect', bx - 0.06, by - 0.06, 0.12, 0.12, fill='white', line='CO', lw=1.25)
    t_half = -3.0 + 0.42 * math.log(0.84 / 0.44 - 1)                           # fA(t_half) = 0.5
    half_mark(sl, tx(t_half), Y1 + H1 / 2, 'CO')
    Y2, H2 = 4.05, 0.85
    shape(sl, 'rect', tx(win[0]), Y2, tx(win[1]) - tx(win[0]), H2, fill='panel2', line=None)
    axes(sl, PX, Y2, PW, H2, '', '')
    text(sl, 6.85, Y2 + H2 / 2 - 0.18, 0.7, 0.36, 'θ_{imp}', size=15, color='muted', align='r')
    gA = lambda t: 1 - sig(t, -3.3, 0.45, 0.08, 0.92)                           # noqa: E731
    gB = lambda t: 1 - sig(t, -4.6, 0.3, 0.08, 0.92)                            # noqa: E731
    curve(sl, gA, t0, t1, PX, Y2, PW, H2, color='CO', lw=2.0)
    curve(sl, gB, win[0], win[1], wx, Y2, ww, H2, color='CO', lw=1.5, dash=True)
    text(sl, PX, Y2 + H2 + 0.03, PW, 0.28, 'r = p(impurity) / p(H₂), log scale', size=12, color='muted', align='c')
    yl = 5.38                                                                   # legend
    circle(sl, 7.12, yl, 0.24, fill=None, line='ink', lw=1.25)
    text(sl, 7.3, yl - 0.15, 1.5, 0.3, 'two steady states', size=12, color='muted')
    line(sl, 8.9, yl, 9.3, yl, color='CO', lw=1.5, dash=True)
    text(sl, 9.36, yl - 0.15, 0.9, 0.3, 'start B', size=12, color='muted')
    half_mark(sl, 10.4, yl, 'CO', d=0.18)
    text(sl, 10.56, yl - 0.15, 1.6, 0.3, 'r_{½} (Π = 0.5)', size=12, color='muted')
    text(sl, 7.0, 5.57, 5.7, 0.28, 'poisoning_curves_vs_T.png, _vs_P_1.png, _vs_P_2.png', size=12, color='muted')
    actions(sl, ['Hold H₂ at the point’s pressure; add one impurity at level r.',
                 'Solve 12 levels, 10⁻⁶ to 0.25, from both starts.',
                 'Divide by the H₂-only rate at that point: Π(r).',
                 'Read off r½, the level that halves the rate.'])
    notes(sl, ['analyses/curves.py: curve_conditions (lines 22-30: per point a "curve r=0" reference and each impurity at',
               'each ratio, on the fixed-probe basis), half_poisoning_ratio (33-44: log-log interpolation and its two',
               'bounds), run (47-71: poisoning_curves.csv and poisoning_tolerance.csv with ppm_of_probe),',
               'plot_poisoning_curves / _poisoning_curves_page (74-141: per point a column, Π on top and θ_imp below,',
               'log-log; a ring where the starts disagree, the start-B branch dashed and the window shaded, an open',
               'diamond at each interpolated r½). θ_imp = 1 - θ* - θ_H - θ_H2* (common.role_coverages, line 164).',
               'Settings: config.py curves (line 83): ratios 1e-6 ... 0.25 (12), impurities None = every impurity (not',
               'the mixture). Figures hold at most 4 columns (plotting.py:24), so the 6 pressure points make _vs_P_1/_2.',
               f'{RUN} check (run.log line 8720): curves, 8 points, 680 solves; poisoning_tolerance.csv has 56 rows',
               '(7 impurities x 8 points). Docs: Stage 11.1 (line 933). The figure is a layout sketch, not results.'])


def s49_blocking(d):
    sl = d.slide('Site blocking: is the poisoning only lost sites?',
                 'At one point (650 K and 1 bar by default), for each impurity that lowers the rate or takes sites.',
                 part=PK, tag='package · curves.py:143–191')
    t0, t1, Y, Ht, Wd = -6.0, -0.6, 2.2, 1.65, 3.25
    xs = (1.0, 5.15, 9.3)
    titles = ('(a) sites', '(b) site-count factor f', '(c) poisoning factor Π')
    for X, title in zip(xs, titles):
        panel(sl, X, Y, Wd, Ht, title, xlab='r, log scale')
    schematic(sl, 11.55, 1.6, w=1.2, label='layout only')
    X = xs[0]                                                                   # (a) empty and impurity-held sites
    fe = lambda t: sig(t, -3.0, 0.5, 0.08, 0.86)                                # noqa: E731
    fi = lambda t: 1 - sig(t, -3.0, 0.5, 0.08, 0.92)                            # noqa: E731
    curve(sl, fe, t0, t1, X, Y, Wd, Ht, color='muted', lw=2.0)
    curve(sl, fi, t0, t1, X, Y, Wd, Ht, color='CO', lw=2.0)
    line(sl, 1.0, 4.3, 1.4, 4.3, color='muted', lw=2.0)
    text(sl, 1.47, 4.15, 1.1, 0.3, 'empty sites', size=12, color='muted')
    line(sl, 2.62, 4.3, 3.02, 4.3, color='CO', lw=2.0)
    text(sl, 3.09, 4.15, 1.3, 0.3, 'impurity-held', size=12, color='muted')
    X = xs[1]                                                                   # (b) f against its reference value
    yf = 0.3
    line(sl, X, Y + Ht * (1 - yf), X + Wd, Y + Ht * (1 - yf), color='muted', lw=1.0, dash=True)
    text(sl, X + 0.08, Y + Ht * (1 - yf) + 0.04, 1.3, 0.3, 'f of H₂ alone', size=12, color='muted')
    ff = lambda t: yf + 0.5 * (1 - sig(t, -3.0, 0.5))                           # noqa: E731
    curve(sl, ff, t0, t1, X, Y, Wd, Ht, color='CO', lw=2.0)
    X = xs[2]                                                                   # (c) the two predictions and the rates
    fs = lambda t: sig(t, -3.0, 0.45, 0.06, 0.9)                                # noqa: E731
    fc = lambda t: sig(t, -2.5, 0.45, 0.06, 0.9)                                # noqa: E731
    curve(sl, fc, t0, t1, X, Y, Wd, Ht, color='CO', lw=1.5, dash=True)
    curve(sl, fs, t0, t1, X, Y, Wd, Ht, color='CO', lw=2.0)
    for t in map(math.log10, LEVELS):
        circle(sl, *at(fs, t, t0, t1, X, Y, Wd, Ht), 0.13, fill='CO', line='white', lw=0.75)
    yl = 4.3                                                                    # legend of (c)
    line(sl, 8.0, yl, 8.45, yl, color='CO', lw=2.0)
    text(sl, 8.52, yl - 0.15, 1.3, 0.3, 'with f', size=12, color='muted')
    line(sl, 9.35, yl, 9.8, yl, color='CO', lw=1.5, dash=True)
    text(sl, 9.87, yl - 0.15, 1.5, 0.3, 'coverages only', size=12, color='muted')
    circle(sl, 11.4, yl, 0.13, fill='CO', line=None)
    text(sl, 11.52, yl - 0.15, 1.3, 0.3, 'from the rates', size=12, color='muted')
    card(sl, 0.55, 4.62, 6.0, 1.26, 'The prediction from the coverages (slide 40)')
    text(sl, 0.75, 5.02, 5.7, 0.8, ['coverages only:   (θ_{∗}/θ_{∗,ref})^{2}, as R0 needs two empty sites',
                                    'with the site count:   × (f/f_{ref})^{2}'], size=14, space_after=3)
    card(sl, 6.78, 4.62, 6.0, 1.26, 'Reading it')
    text(sl, 6.98, 5.02, 5.7, 0.8, ['Dots on the solid line: the impurity only takes sites.',
                                    'Dots off it: it also changes the H₂ step, e.g. by using up H∗.'], size=14,
         space_after=3)
    actions(sl, ['Take the poisoning curves at one point: 650 K, 1 bar.',
                 'Predict Π from the empty sites alone.',
                 'Correct for the change in the site count f.',
                 'Dots on the line: the poisoning is site blocking.'])
    notes(sl, ['analyses/curves.py: blocking_rows (lines 143-153: the curves at curves["blocking_point"], default the',
               'first series point, here 650 K and 1 bar; impurities with min Π < 0.99 or max θ_imp > 0.01 there),',
               'plot_site_blocking / _site_blocking_page (156-191: one row per impurity, at most three per figure:',
               '(a) empty and impurity-held sites, (b) site_count_factor, (c) blocking_full as a line, blocking_from_theta',
               'dashed, the poisoning factor as dots).',
               'The predictions: analyses/common.py add_poisoning (lines 217-233): blocking_from_theta = (θ*/θ*_ref)^n',
               '(p/p_ref), blocking_full = that x (f/f_ref)^n, blocking_check = Π / blocking_full; n = 2 empty sites',
               'for R0; p = p_ref on the fixed-probe basis. Docs: Stage 7 and Stage 11.1. Output: site_blocking.png.',
               'The panels are a layout sketch, not results.'])


def s50_tsweep(d):
    sl = d.slide('Temperature sweep: 12 temperatures at 1 bar, every case',
                 'The temperature series on a 50 K grid: 650 to 1200 K for 9 cases, 108 solves, always at 1 bar.',
                 part=PK, tag='package · sweep.py:20–86')
    Y, Ht, Wd = 2.2, 1.65, 3.25
    xs = (1.0, 5.15, 9.3)
    for X, title, xlab in zip(xs, ('(a) log TOF_{f} against 1000/T', '(b) E_{app} against T', '(c) Π against T'),
                              ('1000/T', 'T', 'T')):
        panel(sl, X, Y, Wd, Ht, title, xlab=xlab)
    schematic(sl, 11.55, 1.6, w=1.2, label='layout only')
    X = xs[0]
    curve(sl, lambda u: 0.88 - 0.38 * u, 0, 1, X, Y, Wd, Ht, color='dark', lw=2.0)
    curve(sl, lambda u: 0.88 - 0.38 * u - 0.42 * u * u, 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    X = xs[1]
    curve(sl, lambda u: 0.2, 0, 1, X, Y, Wd, Ht, color='dark', lw=2.0)
    curve(sl, lambda u: 0.22 + 0.66 * (1 - u) ** 2, 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    X = xs[2]
    line(sl, X, Y + Ht / 2, X + Wd, Y + Ht / 2, color='muted', lw=1.0, dash=True)
    curve(sl, lambda u: sig(u, 0.45, -0.12, 0.08, 0.92), 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    half_mark(sl, X + 0.45 * Wd, Y + Ht / 2, 'CO')
    text(sl, X + 0.45 * Wd + 0.16, Y + Ht / 2 + 0.04, 1.3, 0.3, 'release T', size=12, color='muted')
    yl = 4.3
    line(sl, 1.0, yl, 1.45, yl, color='dark', lw=2.0)
    text(sl, 1.52, yl - 0.15, 1.2, 0.3, 'H₂ only', size=12, color='muted')
    line(sl, 2.6, yl, 3.05, yl, color='CO', lw=2.0)
    text(sl, 3.12, yl - 0.15, 1.6, 0.3, 'an impurity case', size=12, color='muted')
    card(sl, 0.55, 4.55, 7.3, 1.33, 'Computed')
    text(sl, 0.75, 4.93, 7.0, 0.95, ['E_{app}(T) = −R ∂ln TOF/∂(1/T): central differences on the converged points (≥ 3)',
                                     'release T: the lowest T at which a poisoned rate recovers to Π = 0.5',
                                     'A temperature without a steady state is a gap, never bridged.'], size=14)
    card(sl, 8.05, 4.55, 4.73, 1.33, 'Outputs')
    text(sl, 8.25, 4.93, 4.45, 0.95, ['temperature_sweep.csv, release_temperature.csv',
                                      'temperature_sweep.png and _per_area.png',
                                      'sweep_coverages.png: who covers the surface'], size=13)
    actions(sl, ['Solve 9 cases at 12 temperatures, 650–1200\u00a0K, at 1\u00a0bar.',
                 'Plot log TOF against 1000/T.',
                 'Take E_{app} from the local slope.',
                 'Find where each poisoned rate recovers to half.'])
    notes(sl, ['analyses/sweep.py: sweep_temperatures (lines 20-25: None = from the lowest to the highest main temperature',
               'in steps of sweep["step"], 50 K), sweep_conditions (28-33: every case of sweep_cases, the reference, each',
               'impurity and the mixture, at impurity_ratio and reference_pressure), apparent_activation_energy (36-43:',
               'np.gradient in 1/T, cases with >= 3 converged points), release_temperature (46-60), run (63-86),',
               'plot_temperature_sweep (89-137: (a) Arrhenius plot, (b) E_app, (c) Π with release diamonds; per_area: (a)',
               'as kd_eff), plot_sweep_coverages (140-145). Settings: config.py sweep (line 89).',
               f'{RUN} check (run.log line 8721): sweep, 12 points, 108 solves, "temperature only, at reference_pressure".',
               'Docs: Stage 11.2 (line 947). The panels are a layout sketch, not results.'])


def s51_arrhenius(d):
    sl = d.slide('Arrhenius form: the clean metal’s constants, then each case’s',
                 'No solving: fresh surfaces for the H₂ step alone, and fits of the temperature sweep for each case.',
                 part=PK, tag='package · arrhenius.py:75–206')
    card(sl, 0.55, 1.66, 6.1, 3.24, 'Intrinsic: the H₂ step alone')
    grid(sl, 0.9, 2.25, 4, 2, pitch=0.38, d=0.3)
    grid(sl, 3.75, 2.25, 4, 2, pitch=0.38, d=0.3, occ={(r, c): ('H', 'H', 1) for r in range(2) for c in range(4)})
    text(sl, 0.75, 3.08, 2.85, 0.32, 'empty: k_{diss} = q_{f}/Γ', size=13)
    text(sl, 3.6, 3.08, 3.0, 0.32, 'full of H∗: k_{rec} = q_{r}/Γ', size=13)
    text(sl, 0.75, 3.48, 5.8, 0.36, 'K = (k_{diss}/k_{rec}) P°/p,    ln K = ΔS°/R − ΔH°/RT', size=15, font=MATH)
    text(sl, 0.75, 3.88, 5.8, 0.3, 'Per area: k_{d}, k_{r}, γ_{clean} and J_{0}', size=13, color='muted')
    good(sl, 0.72, 4.22, 5.76, 0.6, '`check_K_thermo` = 1 at all 12 temperatures; the γ fit gives back R0’s '
         'A\u00a0=\u00a00.03 and E_{a}\u00a0=\u00a05.0\u00a0kJ/mol.', size=13)
    card(sl, 6.85, 1.66, 5.93, 3.24, 'Observed: each case, from the sweep')
    for y, h, t, size in ((2.1, 0.62, 'Fits of TOF_{f}, TOF_{r} and TOF_{f}/TOF_{r}: A, E_{a}, R², the T range and the '
                                      'conditions', 14),
                          (2.74, 0.62, 'An impurity’s effect in two numbers: ΔE_{a} = E_{a} − E_{a,ref} and A/A_{ref}', 14),
                          (3.38, 0.62, 'R² below 0.99: flagged curved; read the local E_{a}(T) in arrhenius_local.csv', 14),
                          (4.02, 0.85, 'Per area: k_{d,eff}, J_{ex}, γ_{eff}. surface_coefficients.csv lists k_{d} and k_{r} '
                                       'as A\u00a0exp(−E_{a}/RT): the clean metal and the 9 cases.', 13)):
        text(sl, 7.05, y, 5.6, h, t, size=size)
    caveat(sl, 0.55, 5.0, W - 2 * LM, 0.88, 'For try_2/Cu the intrinsic k_{diss} fit is flagged curved (R² = 0.95), '
           'yet nothing changes on an empty surface: k_{diss} ∝ γ/√T, with γ = 0.03 e^{−5 kJ/mol / RT}. Its local '
           'E_{a} = 5.0 − RT/2 falls from 2.3 to 0 kJ/mol over 650–1200 K; the fitted 1.4 kJ/mol is an average.', size=14)
    actions(sl, ['Evaluate the H₂ step on an empty and a full surface.',
                 'Fit A and E_{a} of each rate; K gives ΔH° and ΔS°.',
                 'Fit each case’s rates from the temperature sweep.',
                 'Compare with H₂ alone: ΔE_{a} and A/A_{ref}.'])
    notes(sl, ['analyses/arrhenius.py: fit (lines 47-62: ln y = ln A - Ea/RT by least squares; R² left empty for a flat',
               'quantity), local_ea (65-72), intrinsic_table (75-108: on fresh surfaces at reference_pressure, k_diss on',
               'theta* = 1 and k_rec on theta_H = 1; K = k_diss/k_rec x P0/p; ΔS from the thermo converted from 1 atm to',
               '1 bar; check_K_thermo), run (145-206: fits, dEa_vs_reference and A_over_reference, curved = R² < min_r2 =',
               '0.99), intrinsic_quantities (216-228), surface_coefficients (237-257). Settings: config.py arrhenius',
               '(line 95): temperatures None = the sweep temperatures.',
               f'{RUN} (output/arrhenius_intrinsic.csv, arrhenius_fits.csv): check_K_thermo = 1.0 at 650 ... 1200 K;',
               'sticking probability fit A = 0.03, Ea = 5.000 kJ/mol, R² = 1 (R0 in the YAML: sticking-coefficient',
               'A 0.03, b 0, Ea 5000 J/mol); k_diss fit Ea = 1.37 kJ/mol, R² = 0.951, so "curved". k_diss = γ x the',
               'collision rate p/sqrt(2 pi M R T) per site, so ∝ γ/√T and the local Ea = 5.0 - RT/2 kJ/mol: 2.30 at 650 K,',
               '0.01 at 1200 K (RT/2 = 2.70 and 4.99 kJ/mol).',
               'surface_coefficients.csv: 10 rows. Docs: Stage 11.2a (line 968).'])


def s52_psweep(d):
    sl = d.slide('Pressure sweep: p(H₂) from 0.01 to 100 bar at fixed purity',
                 'Every impurity keeps its ratio to H₂, so only the total pressure changes: 2 × 9 × 9 = 162 solves.',
                 part=PK, tag='package · pressure_sweep.py:29–121')
    card(sl, 0.55, 1.66, 4.75, 4.22, 'The grid: 2 temperatures × 9 pressures')
    X0, X1 = 1.55, 4.9
    for j, (T, y) in enumerate(((650, 2.45), (1200, 2.95))):
        text(sl, 0.65, y - 0.16, 0.8, 0.32, f'{T} K', size=13, color='muted', align='r')
        line(sl, X0, y, X1, y, color='rule', lw=0.75)
        for k in range(9):
            circle(sl, X0 + k * (X1 - X0) / 8, y, 0.17, fill='dark', line=None)
    for k, lab in ((0, '0.01'), (4, '1'), (8, '100 bar')):
        text(sl, X0 + k * (X1 - X0) / 8 - 0.5, 3.12, 1.0, 0.28, lab, size=12, color='muted', align='c')
    text(sl, 0.75, 3.55, 4.4, 0.35, 'p(H₂) on a log grid, 2 per decade; × 9 cases = **162\u00a0solves**', size=14)
    text(sl, 0.75, 4.15, 4.4, 0.62, 'Each impurity keeps its ratio to H₂: the same purity at every pressure.', size=14)
    text(sl, 0.75, 4.85, 4.4, 0.95, 'pressure_sweep.csv, half_rate_pressure.csv; pressure_sweep.png and _per_area.png '
         '(J_{diss}); surface_exchange_vs_P.png; pressure_sweep_coverages_650K.png and _1200K.png', size=12,
         color='muted')
    Y, Ht, Wd = 2.2, 1.45, 1.95
    xs = (5.95, 8.3, 10.65)
    for X, title in zip(xs, ('(a) TOF, log–log', '(b) local slope', '(c) Π, and p_{½}')):
        panel(sl, X, Y, Wd, Ht, title, xlab='p(H₂), log')
    schematic(sl, 11.6, 3.98, w=1.15, label='layout only')
    X = xs[0]
    curve(sl, lambda u: 0.12 + 0.72 * u, 0, 1, X, Y, Wd, Ht, color='dark', lw=2.0)
    curve(sl, lambda u: 0.08 + 0.75 * u - 0.55 * u * u, 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    X = xs[1]
    for v, lab in ((1, '+1'), (0, '0'), (-1, '−1')):                              # slope axis from -1.2 to +1.2
        yv = Y + Ht * (1 - (v + 1.2) / 2.4)
        line(sl, X, yv, X + Wd, yv, color='rule' if v else 'muted', lw=0.75)
        text(sl, X - 0.5, yv - 0.15, 0.42, 0.3, lab, size=12, color='muted', align='r')
    curve(sl, lambda u: (1 - 0.35 * u + 1.2) / 2.4, 0, 1, X, Y, Wd, Ht, color='dark', lw=2.0)
    curve(sl, lambda u: (1 - 1.9 * u ** 1.4 + 1.2) / 2.4, 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    X = xs[2]
    line(sl, X, Y + Ht / 2, X + Wd, Y + Ht / 2, color='muted', lw=1.0, dash=True)
    curve(sl, lambda u: sig(u, 0.5, 0.13, 0.08, 0.92), 0, 1, X, Y, Wd, Ht, color='CO', lw=2.0)
    half_mark(sl, X + 0.5 * Wd, Y + Ht / 2, 'CO')
    yl = 4.12
    line(sl, 5.95, yl, 6.4, yl, color='dark', lw=2.0)
    text(sl, 6.47, yl - 0.15, 1.2, 0.3, 'H₂ only', size=12, color='muted')
    line(sl, 7.5, yl, 7.95, yl, color='CO', lw=2.0)
    text(sl, 8.02, yl - 0.15, 1.6, 0.3, 'an impurity case', size=12, color='muted')
    card(sl, 5.55, 4.45, 7.23, 1.43, 'For a two-site step, the slope tends to')
    text(sl, 5.75, 4.84, 6.9, 1.0, ['+1 on an almost empty surface',
                                    '0 when H∗ fills it, θ_{∗} ∝ p^{−1/2}',
                                    '−1 when a strongly bound impurity holds it at fixed purity, θ_{∗} ∝ 1/p'], size=14)
    actions(sl, ['At 650 and 1200 K, step p(H₂) over 9 pressures.',
                 'Scale every impurity with H₂: the same purity.',
                 'Take the local slope d\u00a0ln\u00a0TOF\u00a0/\u00a0d\u00a0ln\u00a0p.',
                 'Find where Π crosses 0.5, and in which direction.'])
    notes(sl, ['analyses/pressure_sweep.py: sweep_pressures (lines 29-41: a log grid with per_decade = 2 points per decade',
               'spanning the pressure series and reference_pressure: 1e3 to 1e7 Pa, 9 pressures), sweep_temperatures (44-46:',
               'pressure_temperatures, 650 and 1200 K), sweep_conditions (49-58), local_slope (61-69: np.gradient of ln TOF',
               'in ln p), half_rate_pressure (72-87: first crossing of Π = 0.5 going up, with "recovers above half" or',
               '"falls below half"), run (90-121), plot_pressure_sweep (124-179: one row per temperature).',
               'Settings: config.py pressure_sweep (line 99). The slope is the order in H2 for the H2-only case, and the',
               'order in the total pressure at fixed purity for an impurity case.',
               f'{RUN} check (run.log line 8723): pressure_sweep, 18 points, 162 solves; half_rate_pressure.csv 16 rows.',
               'Docs: Stage 11.2b (line 1004). The panels are a layout sketch, not results.'])


def s53_equilibrium(d):
    sl = d.slide('Equilibrium coverage and q: is the H₂ step at equilibrium?',
                 'From the thermochemistry alone, with no solving; then every steady state is checked against it.',
                 part=PK, tag='package · equilibrium.py:39–152')
    card(sl, 0.55, 1.66, 6.1, 4.22, 'Where H₂ + 2∗ ⇌ 2H∗ alone puts the surface')
    text(sl, 0.75, 2.1, 5.8, 0.38, 'K = e^{−ΔG°/RT} P°/P_{thermo}', size=16, font=MATH)
    text(sl, 0.75, 2.5, 5.8, 0.3, 'RMG’s thermochemistry is at 1 atm; K is at P° = 1 bar.', size=13, color='muted')
    for y, t in ((2.9, 'x = (K p/P°)^{1/2}'), (3.27, 'θ_{H,eq} = x/(1 + x),    θ_{v,eq} = 1/(1 + x)'),
                 (3.64, 'C_{H,eq} = θ_{H,eq} Γ,    J_{0} = k_{d} p θ_{v,eq}^{2}')):
        text(sl, 0.75, y, 5.8, 0.36, t, size=15, font=MATH)
    text(sl, 0.75, 4.02, 5.8, 0.32, 'T_{½}: where K p/P° = 1, searched from 200 to 3000 K', size=14)
    shape(sl, 'rrect', 0.72, 4.42, 5.76, 0.8, fill='white', line='rule', radius=0.08)
    text(sl, 0.85, 4.45, 5.55, 0.74, 'At 650 K, K = 5.5 × 10^{−4} (slide 16). At 1 bar, x = 0.023 and θ_{H,eq} = 0.023: '
         'clean Cu(111) holds little H.', size=14, anchor='m')
    text(sl, 0.75, 5.35, 5.8, 0.45, 'Grid: 4 temperatures × 17 pressures, 0.01 to 100 bar, 4 per decade', size=13,
         color='muted')
    card(sl, 6.85, 1.66, 5.93, 4.22, 'The check q, at every converged steady state')
    text(sl, 7.05, 2.1, 5.6, 0.4, 'q = (θ_{H}/θ_{v}) / (K p/P°)^{1/2}', size=16, font=MATH)
    text(sl, 7.05, 2.55, 5.6, 0.62, 'An impurity lowers θ_{H} and θ_{v} together; at equilibrium their ratio stays fixed.',
         size=13)
    for i, (q, key, t) in enumerate((('q = 1', 'ok', 'the H₂ step is at equilibrium: Π = (θ_{∗}/θ_{∗,ref})^{2}'),
                                     ('q < 1', 'dark', 'H∗ is used faster than it is made: the H₂ step limits the rate'),
                                     ('q > 1', 'dark', 'other reactions make H∗, and it leaves as H₂'))):
        y = 3.05 + i * 0.6
        chip(sl, 7.05, y + 0.06, q, key=key, w=0.85, h=0.34, size=13)
        text(sl, 8.05, y, 4.6, 0.55, t, size=13, anchor='m')
    text(sl, 7.05, 4.95, 5.6, 0.32, 'try_2/Cu: 342 steady states, 72 main and 108 + 162 from the sweeps', size=13,
         color='muted')
    text(sl, 7.05, 5.35, 5.6, 0.32, 'equilibrium_coverage.png: (a) isotherms, (b) against T, (c) q', size=12,
         color='muted')
    actions(sl, ['Take K of the H₂ step from the thermochemistry.',
                 'Solve θ_{H}/θ_{v} = √(Kp/P°) on a T and p grid.',
                 'Compute q for every steady state of the run.',
                 'q = 1: at equilibrium; otherwise the H₂ step is driven.'])
    notes(sl, ['analyses/equilibrium.py: K_thermo (lines 39-45: exp(-ΔG/RT) of H2 + 2* <=> 2 H*, converted from the thermo',
               'reference pressure to P0 = 1 bar), coverage (48-51), half_coverage_temperature (54-73: scan 200-3000 K in 10 K',
               'steps, then bisection), grid (76-87: temperatures = the temperature series and pressure_temperatures;',
               'pressures = a log grid spanning the pressure series with per_decade = 4), check_table (102-117: q for every',
               'converged steady state of steady_state, sweep and pressure_sweep), run (120-152), plot_equilibrium (155-218).',
               'J0 uses kd from per_area.clean_surface (per_area.py line 88). Settings: config.py equilibrium (line 106).',
               f'{RUN}: equilibrium_coverage.csv 68 rows (4 x 17); half_coverage_temperature.csv at 0.01, 1 and 100 bar;',
               'equilibrium_check.csv 342 rows. Worked example: K = 5.5e-4 at 650 K (slide 16), x = sqrt(5.5e-4) = 0.0235,',
               'θ_H,eq = 0.0229 (arrhenius_intrinsic.csv theta_H_eq at 650 K: 0.02295). Docs: Stage 11.2c (line 1015).'])


def s54_orders(d):
    sl = d.slide('Reaction orders: nudge one partial pressure by ±10 %',
                 'n_{j} = ∂ln TOF/∂ln p_{j} by central differences, every other pressure fixed: 304 solves for try_2/Cu.',
                 part=PK, tag='package · orders.py:17–68')
    card(sl, 0.55, 1.66, 6.1, 4.22, 'Central differences, per case and point')
    text(sl, 0.75, 2.1, 5.8, 0.4, 'n_{j} = Δ ln TOF_{f} / Δ ln p_{j}', size=16, font=MATH)
    text(sl, 0.75, 2.52, 5.8, 0.3, 'p_{j} × 1.1 and × 0.9 (δ = 0.10), every other pressure fixed', size=13, color='muted')
    rows = (('H₂ only', [('base', 0.7), ('H₂ +10 %', 1.05), ('H₂ −10 %', 1.05)], '3 solves'),
            ('each impurity', [('base', 0.7), ('H₂ ±10 %', 1.05), ('impurity ±10 %', 1.45)], '5 solves'))
    for i, (lab, chips, n) in enumerate(rows):
        y = 3.0 + i * 0.55
        text(sl, 0.75, y, 1.4, 0.34, lab, size=13, bold=True, anchor='m')
        x = 2.15
        for t, w in chips:
            chip(sl, x, y, t, w=w, h=0.34, size=12)
            x += w + 0.1
        text(sl, x + 0.05, y, 1.0, 0.34, n, size=13, color='muted', anchor='m')
    text(sl, 0.75, 4.2, 5.8, 0.62, '8 points × (3 + 7 × 5) = **304 solves**, from both starts; none for “all”', size=14)
    text(sl, 0.75, 4.9, 5.8, 0.85, 'The log span is 0.2, so a TOF tolerance of 10^{−4} gives n to about 10^{−3}.', size=14)
    card(sl, 6.85, 1.66, 5.93, 4.22, 'Reading them: TOF_{f} ∝ p_{H₂} θ_{∗}^{2}')
    text(sl, 7.05, 2.15, 3.4, 0.3, 'the surface', size=13, color='muted', bold=True)
    text(sl, 10.45, 2.15, 1.0, 0.3, 'n_{H₂}', size=13, color='muted', bold=True, align='c')
    text(sl, 11.55, 2.15, 1.0, 0.3, 'n_{I}', size=13, color='muted', bold=True, align='c')
    table = (('almost empty', '+1', '0'), ('H∗ fills it: θ_{∗} ∝ p_{H₂}^{−1/2}', '0', '0'),
             ('the impurity holds it: θ_{∗} ∝ 1/p_{I}', '+1', '−2'), ('they share it', 'between', 'between'))
    for i, (what, nh, ni) in enumerate(table):
        y = 2.55 + i * 0.6
        line(sl, 7.05, y - 0.04, 12.6, y - 0.04, color='rule', lw=0.75)
        text(sl, 7.05, y, 3.4, 0.55, what, size=14, anchor='m')
        text(sl, 10.45, y, 1.0, 0.55, nh, size=14, bold=True, align='c', anchor='m')
        text(sl, 11.55, y, 1.0, 0.55, ni, size=14, bold=True, align='c', anchor='m')
    text(sl, 7.05, 5.1, 5.6, 0.7, 'reaction_orders_vs_T.png and _vs_P.png: a panel per point; ● the order in H₂, '
         '■ in the impurity', size=12, color='muted')
    actions(sl, ['Raise one partial pressure by 10\u00a0%, then lower it by 10\u00a0%.',
                 'Hold every other pressure fixed.',
                 'n = the change in ln TOF over the change in ln p.',
                 'Read the regime from n in H₂ and in the impurity.'])
    notes(sl, ['analyses/orders.py: order_conditions (lines 17-35: per point and case, base and probe +/- delta, and for an',
               'impurity case also impurity +/- delta; cases = the reference and each impurity, not the mixture), run',
               '(38-68: n = (ln TOF+ - ln TOF-) / (ln(1+δ) - ln(1-δ)); all_starts_agree), plot_orders (71-116: one panel per',
               'point, the same y scale everywhere, values labelled where |n| >= 0.1). Settings: config.py orders (line',
               '111): delta 0.10. Error estimate: docs Stage 11.3 (line 1040), with the reactor-method TOF tolerance 1e-4.',
               f'{RUN} check (run.log line 8726): orders, 8 points, 304 solves = 8 x (3 + 7 x 5); reaction_orders.csv 64',
               'rows (8 cases x 8 points). The two-site readings follow from TOF_f ∝ p_H2 θ*^2 (docs 11.3).'])


def s55_drc(d):
    sl = d.slide('Degree of rate control: which step’s speed sets the rate',
                 'Speed one reaction up and down by 5 %, forward and reverse together, so its K does not change.',
                 part=PK, tag='package · rate_control.py:40–141')
    card(sl, 0.55, 1.66, 6.1, 4.22, 'Kinetic DRC of reaction i')
    text(sl, 0.75, 2.1, 5.8, 0.42, 'X_{DRC,i} = ∂ln r / ∂ln k_{i}', size=17, font=MATH)
    text(sl, 0.75, 2.55, 5.8, 0.32, 'k_{f,i} and k_{r,i} × 1.05 and × 0.95 together: K_{i} unchanged', size=13,
         color='muted')
    text(sl, 0.75, 3.0, 5.8, 0.62, '**Screened:** R0, the H₂ adsorption #1 (H₂ + ∗ ⇌ H₂∗), and the 20 largest '
         '|net rates| at the base', size=14)
    text(sl, 0.75, 3.7, 5.8, 0.32, '**Two rates r:**', size=14)
    text(sl, 0.75, 4.04, 5.8, 0.62, ['the forward TOF, a gross rate: its X need not sum to 1',
                                     'the net H₂ consumption: X sums to 1; blank near equilibrium'], size=14)
    text(sl, 0.75, 5.0, 5.8, 0.75, 'A _DRC.png figure only where the DRC is not trivial (trivial: X = 1 on R0, '
         'about 0 elsewhere)', size=13, color='muted')
    card(sl, 6.85, 1.66, 5.93, 4.22, 'Every perturbed solve starts from a clean surface')
    tick(sl, 7.2, 2.36, ok=False)
    text(sl, 7.5, 2.12, 5.15, 0.62, 'A warm start from the base passed the plateau test after about '
         '3\u00a0×\u00a010^{−8}\u00a0s, before slow species moved: X ≈ 0, wrongly.', size=13)
    tick(sl, 7.2, 3.0, ok=True)
    text(sl, 7.5, 2.85, 5.15, 0.32, 'A clean start passes only once the slow species have settled.', size=13)
    line(sl, 7.05, 3.42, 12.6, 3.42, color='rule', lw=0.75)
    text(sl, 7.05, 3.55, 5.6, 0.62, 'H₂ only + 7 impurities (not “all”) × 8 points = **64 base solves**, then 2 per '
         'screened reaction and species, from one start each', size=14)
    text(sl, 7.05, 4.3, 5.6, 0.62, 'check estimated 3,392; the run did 64 + 2,992 = 3,056.', size=14)
    text(sl, 7.05, 5.2, 5.6, 0.32, 'Output: sensitivity.csv, 1,496 rows (DRC and TRC)', size=12, color='muted')
    actions(sl, ['Solve the base steady state from a clean surface.',
                 'Scale one reaction’s k_{f} and k_{r} by 1.05, then by 0.95.',
                 'Solve each again from a clean surface.',
                 'X = the change in ln r over the change in ln k.'])
    notes(sl, ['analyses/rate_control.py: _solve (lines 29-37: every solve from S.clean_surface, with the job\'s multipliers',
               'and h_shifts), drc_reactions (40-46: the probe and molecular-adsorption reactions plus the drc_top_n = 20',
               'largest |net rates of progress| at the base), run (63-141: cases = reference + each impurity; 2 jobs per',
               'screened reaction (1 +/- drc_delta = 0.05) and per screened species; X_probe_consumed left NaN where the net',
               'consumption is below 1e-6 of the forward TOF), drc_nontrivial (144-152), plot_rate_control (168-241).',
               'Multipliers: surf.set_multiplier scales forward and reverse together (model.py new_surface, 185-211).',
               f'{RUN} run.log: "64 to solve" for the base states (line 1849), then "2992 perturbed solves" (line 1915);',
               'check (line 8727): rate_control, 8 points, 3392, "about 53 solves per case and point; estimate".',
               'sensitivity.csv: 1,496 rows (DRC and TRC). Warm-start finding: docs Stage 11.4 (line 1057), "Why a warm',
               'start failed". R0 = #0 and the molecular adsorption #1 H2(1) + vacantX(8) <=> [H][H].[Pt](34) (run.log l. 6-7).'])


def s56_trc(d):
    sl = d.slide('Thermodynamic rate control: which species’ stability sets it',
                 'Poisoning is set by how strongly species bind, not by how fast steps are, and DRC misses it.',
                 part=PK, tag='package · rate_control.py:49–60')
    card(sl, 0.55, 1.66, 6.1, 4.22, 'Shift one species, keep the barriers')
    level(sl, 1.0, 3.35, 1.1, 'species n', color='ink', label_side='above')
    level(sl, 2.85, 2.6, 1.1, 'transition state', color='ink', label_side='above')
    level(sl, 4.75, 3.05, 1.1, 'product', color='ink', label_side='above')
    poly(sl, [(2.1, 3.35), (2.85, 2.6)], color='faint', lw=1.25)
    poly(sl, [(3.95, 2.6), (4.75, 3.05)], color='faint', lw=1.25)
    line(sl, 1.0, 3.67, 2.1, 3.67, color='H', lw=2.0, dash=True)
    arrow(sl, 1.55, 3.41, 1.55, 3.62, color='H', lw=1.25, head='sm')
    text(sl, 2.2, 3.52, 2.4, 0.3, '±δ = 1 kJ/mol', size=12, color='H', bold=True)
    text(sl, 2.6, 2.66, 1.6, 0.3, 'held fixed', size=12, color='muted', align='c')
    text(sl, 0.75, 4.0, 5.8, 0.62, '1. Shift H_{n} by ±δ: a_{6} of both NASA ranges (new_surface with h_shifts)', size=14)
    text(sl, 0.75, 4.65, 5.8, 0.62, '2. Hold the barriers: each reaction with n as a reactant × e^{ν′δ/RT}', size=14)
    text(sl, 0.75, 5.33, 5.8, 0.42, 'X_{TRC,n} = ∂ln r / ∂(−G_{n}/RT)', size=16, font=MATH)
    card(sl, 6.85, 1.66, 5.93, 4.22, 'Why it is needed, and how to read it')
    text(sl, 7.05, 2.1, 5.6, 0.95, 'DRC keeps every K fixed, so it gives about 0 for an adsorption step at '
         'equilibrium: the kind of step that decides poisoning.', size=14)
    text(sl, 7.05, 2.95, 5.6, 0.62, 'Screened: ∗, H∗, H₂∗, and every species with θ ≥ 10^{−3} at the base', size=14)
    shape(sl, 'rrect', 7.02, 3.6, 5.6, 0.95, fill='white', line='rule', radius=0.08)
    text(sl, 7.17, 3.63, 5.35, 0.9, 'If a species only takes sites, X_{TRC,k} ≈ −2θ_{k}: the values add up to about '
         '−2(1 − θ_{∗}).', size=14, anchor='m')
    text(sl, 7.05, 4.75, 5.6, 0.85, 'rate_control_<T>K_<P>bar.png: a panel per case, species with |X| ≥ 0.01, the '
         'sum in the axis label', size=12, color='muted')
    actions(sl, ['Pick a role species or one with θ ≥ 10⁻³.',
                 'Shift its enthalpy by ±1\u00a0kJ/mol.',
                 'Hold the transition states with rate multipliers.',
                 'X_{TRC} ≈ −2θ when it only takes sites.'])
    notes(sl, ['analyses/rate_control.py: trc_species (lines 49-53: the empty site, H*, H2* and every species with coverage',
               '>= trc_min_coverage = 1e-3 at the base), trc_multipliers (56-60: each reaction with the species as a reactant',
               'x exp(nu\' dG / RT), so the transition state stays where it was), run (63-141: X = (ln r+ - ln r-) / span,',
               'span = -2 trc_delta / RT, trc_delta = 1e6 J/kmol = 1 kJ/mol), _trc_page (185-210).',
               'The enthalpy shift: model.py new_surface (lines 185-211): a6 of both NASA7 ranges += dh/R (h0 for constant cp),',
               'then modify_species. Settings: config.py rate_control (line 112). Docs: Stage 11.4 (line 1057), the',
               'reading X_TRC,i ≈ -n θ_i with n = 2 for R0; why a 1 kJ/mol shift is safe, unlike editing thermo in place',
               '(Stage 10.5). The energy diagram is schematic.'])


def s57_paths(d):
    sl = d.slide('Reaction paths: where the H atoms go at each steady state',
                 'Cantera’s element-flux diagram of each solved main steady state. No new solving.',
                 part=PK, tag='package · paths.py:137–203')
    card(sl, 0.55, 1.66, 5.95, 4.22, 'An H-atom flux diagram, one case')
    nodes = {'H₂': (1.75, 2.78, 'H', 4.5), 'H∗': (5.0, 2.78, 'H', 2.25), 'H₂∗': (1.75, 4.55, 'H', 0.75),
             'H₂O∗': (5.0, 4.55, 'H2O', 1.75)}
    for t, (cx, cy, k, lw_) in nodes.items():
        shape(sl, 'rrect', cx - 0.5, cy - 0.26, 1.0, 0.52, fill='white', line='ink', lw=lw_, radius=0.08, content=t,
              size=15, color=k, bold=True)
    arrow(sl, 1.6, 3.1, 1.6, 4.23, color='ink', lw=6, head='sm')             # H2 <-> H2*: the large exchange
    arrow(sl, 1.9, 4.23, 1.9, 3.1, color='ink', lw=6, head='sm')
    arrow(sl, 2.3, 2.67, 4.45, 2.67, color='ink', lw=0.75, head='sm')        # R0, both ways
    arrow(sl, 4.45, 2.89, 2.3, 2.89, color='ink', lw=0.75, head='sm')
    arrow(sl, 4.85, 3.1, 4.85, 4.23, color='ink', lw=2.5, head='sm')         # H* <-> H2O*
    arrow(sl, 5.15, 4.23, 5.15, 3.1, color='ink', lw=2.5, head='sm')
    text(sl, 0.7, 3.52, 0.8, 0.3, '1', size=12, color='muted', align='r')
    text(sl, 2.85, 2.3, 1.1, 0.3, '10^{−5}', size=12, color='muted', align='c')
    text(sl, 5.27, 3.52, 0.9, 0.3, '10^{−2}', size=12, color='muted')
    text(sl, 0.75, 5.0, 5.6, 0.32, 'labels: flux / largest flux; outline width ~ log amount', size=12, color='muted')
    schematic(sl, 0.75, 5.45, w=5.5, label='schematic, not a computed diagram')
    card(sl, 6.83, 1.66, 5.95, 4.22, 'The defaults, and why')
    rows = (('element', 'H, the element of H₂'),
            ('one-way flow', 'both directions: near equilibrium the net flows cancel'),
            ('threshold', '10^{−8} of the largest flux: the H₂ step can carry only 10^{−5}–10^{−6}'),
            ('rmg style', 'black arrows, width ~ log flux; outline ~ log coverage or mole fraction'),
            ('images', 'structures from species/ next to settings.py'),
            ('files', '72 diagrams, 9 cases × 8 points: .dot and .png'))
    for i, (k, t) in enumerate(rows):
        y = 2.1 + i * 0.62
        text(sl, 7.0, y, 1.55, 0.58, k, size=13, bold=True)
        text(sl, 8.6, y, 4.05, 0.58, t, size=13)
    actions(sl, ['Rebuild a solved main steady state.',
                 'Ask Cantera for the H-atom fluxes between species.',
                 'Keep both directions above 10⁻⁸ of the largest.',
                 'Draw in RMG’s style, with structure images.'])
    notes(sl, ['analyses/paths.py: default_conditions (lines 22-24: every case at every point of the series), raw_dot',
               '(137-158: common.load_state, ct.ReactionPathDiagram(surf, element), threshold, flow_type; amounts = coverage',
               'or mole fraction), style_dot (75-134), render (161-171: Graphviz dot), diagram (174-183), run (186-203:',
               'elements None = those of the probe gas). Settings: config.py paths (line 121): threshold 1e-8, flow_type',
               '"OneWayFlow", style "rmg", species_images None = species/ next to settings.py, edge_labels True.',
               f'{RUN} run.log line 10171: structure images from runs/try_2/Cu/species; output/reaction_paths holds 144',
               'files = 72 diagrams x (.dot, .png), path_H_<case>_<T>K_<P>bar_OneWayFlow. The "Scale" under each diagram',
               'is the largest flux in kmol m-2 s-1. Docs: Stage 11.5 (line 1096). The diagram on the slide is schematic.'])


def s58_transients(d):
    sl = d.slide('Poisoning over time: switch the gas at t = 0, follow the rate',
                 'Start from a catalyst working in clean H₂, add the impurity, and record how fast the rate falls.',
                 part=PK, tag='package · transients.py:24–150')
    card(sl, 0.55, 1.66, 6.45, 4.22, 'One run: the rate after the gas switch')
    schematic(sl, 5.75, 1.78, w=1.15, label='layout only')
    X, Y, Wd, Ht = 1.45, 2.55, 5.15, 2.65
    axes(sl, X, Y, Wd, Ht, '', '')
    text(sl, X + 0.08, Y - 0.32, 3.4, 0.28, 'TOF / TOF(t = 0), log scale', size=12, color='muted')
    text(sl, X, Y + Ht + 0.04, Wd, 0.28, 'time after the switch, log: 10^{−9} to 10^{8} s', size=12, color='muted',
         align='c')
    logistic = lambda v: 1 / (1 + math.exp(-v))                                 # noqa: E731
    L = lambda u: -2.4 * logistic((u - 0.45) / 0.07) + 0.92 * logistic((u - 0.72) / 0.04)    # noqa: E731
    fr = lambda u: 0.92 + 0.3 * L(u)                                # noqa: E731  (log10 of the ratio, to 0..1)
    for v, lab in ((0.5, '½'), (0.1, '1/10')):
        yv = Y + Ht * (1 - (0.92 + 0.3 * math.log10(v)))
        line(sl, X, yv, X + Wd, yv, color='muted', lw=1.0, dash=True)
        text(sl, X + Wd - 0.6, yv - 0.29, 0.55, 0.28, lab, size=12, color='muted', align='r')
    curve(sl, fr, 0, 1, X, Y, Wd, Ht, n=120, color='CO', lw=2.25)
    marks = []
    for v, lab in ((0.5, 't_{½}'), (0.1, 't_{1/10}')):                         # first crossings, from the curve
        u = 0.45 + 0.07 * math.log((-math.log10(v) / 2.4) / (1 + math.log10(v) / 2.4))
        marks.append((u, lab))
    marks.append((0.84, 't_{settle}'))
    for u, lab in marks:
        cx, cy = at(fr, u, 0, 1, X, Y, Wd, Ht)
        line(sl, cx, cy, cx, Y + Ht - 0.34, color='faint', lw=0.75, dash=True)
        text(sl, cx - 0.4, Y + Ht - 0.32, 0.8, 0.28, lab, size=12, color='ink', align='c')
    u_min = min((i / 200 for i in range(201)), key=L)
    cx, cy = at(fr, u_min, 0, 1, X, Y, Wd, Ht)
    circle(sl, cx, cy, 0.16, fill='ink', line=None)
    text(sl, cx - 0.6, cy + 0.1, 1.2, 0.28, 'deepest drop', size=12, color='ink', align='c')
    text(sl, X + Wd - 0.75, Y + Ht * (1 - fr(1.0)) + 0.06, 0.7, 0.28, 'end', size=12, color='ink', align='r')
    card(sl, 7.2, 1.66, 5.58, 4.22, 'The method, for try_2/Cu')
    for y, h, t, size in ((2.1, 0.36, 'Start: the H₂-only steady state at the point, from the main run.', 14),
                          (2.55, 0.62, 't = 0: the gas becomes H₂ + impurity; the fixed-gas reactor integrates the '
                                       'surface.', 14),
                          (3.25, 0.62, '69 times from 10^{−9} to 10^{8} s, 4 per decade; TOF divided by its value at '
                                       't = 0.', 14),
                          (3.95, 0.36, '8 points × 7 impurities × 4 levels (10^{−6}, 10^{−4}, 10^{−2}, 0.25) = '
                                       '**224 runs**', 14),
                          (4.45, 0.62, 'Quote times to 2 significant figures: a 10^{−10} change in the start moved one '
                                       'half-time by 1.7 %.', 13),
                          (5.2, 0.36, 'poisoning_over_time.csv, _trajectories.csv', 12)):
        text(sl, 7.4, y, 5.25, h, t, size=size, color='muted' if size < 14 else 'ink')
    actions(sl, ['Start from the H₂-only steady state at the point.',
                 'At t = 0, add the impurity to the gas.',
                 'Integrate the surface; record 69 times.',
                 'Read t½, t₁/₁₀, the deepest drop and the end.'])
    notes(sl, ['analyses/transients.py: times (lines 24-27: logspace from t_start 1e-9 to t_stop 1e8 s, per_decade 4: 69',
               'times), run_transient (30-56: set_state, solver.fixed_gas_reactor with the gas multiplier 0, tof0 at t = 0',
               'before integration, then net.advance through the times), crossing_time (59-69), settle_time (72-79: within',
               '10 % of the final value), reference_start (82-94: the reference steady state of the main run at the point),',
               'run (97-150: poisoning_over_time.csv with t_half_s, t_tenth_s, min_over_t0, t_min_s, end_over_t0,',
               't_settle_s), plot_poisoning_over_time (153-206: one figure per point, at most 4 impurities each: _1, _2).',
               'Settings: config.py transients (line 136): ratios [1e-6, 1e-4, 1e-2, 0.25], impurities None = each one.',
               f'{RUN} check (run.log line 8728): transients, 8 points, 224 = 8 x 7 x 4; poisoning_over_time.csv 224 rows.',
               'Precision: docs Stage 11.6 (line 1116), "Precision of the times". The curve is a layout sketch.'])


SLIDES = [s48_curves, s49_blocking, s50_tsweep, s51_arrhenius, s52_psweep, s53_equilibrium, s54_orders, s55_drc,
          s56_trc, s57_paths, s58_transients]
