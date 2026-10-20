"""Section G (slides 28-35): Stage 4, solving one steady state: the problem, the method choice,
set_state, the direct and reactor methods, the plateau test, the convergence test, and the
convergence traces. Running example: runs/try_2/Cu."""
import math

from deck_lib import *  # noqa: F401,F403

PG = 'G · Solving one steady state · Stage 4'
RUN = 'runs/try_2/Cu'


def s28_problem(d):
    sl = d.slide('Stage 4: a steady state is 71 coverages that make every net rate zero',
                 'Two ways to find it: Newton’s method jumps there; time integration walks there.',
                 part=PG, tag='concept')
    card(sl, 0.55, 1.66, 5.5, 4.22, None)
    text(sl, 0.75, 1.75, 2.4, 0.35, '71 unknowns', size=16, bold=True)
    text(sl, 3.3, 1.75, 2.6, 0.35, '71 equations', size=16, bold=True)
    names = [('vacantX(8)', 'vacantX'), ('[PtH](33)', '[PtH]'), ('[H][H].[Pt](34)', '[H][H]'), ('OCX(15)', 'OCX'),
             ('⋮', None)]
    for i, (n, short) in enumerate(names):
        y = 2.25 + i * 0.5
        text(sl, 0.75, y, 2.5, 0.4, ('θ  ' + n) if short else '   ⋮', size=14, font=CODE)
        text(sl, 3.3, y, 2.6, 0.4, (f'ṡ({short}) = 0') if short else '   ⋮', size=14, font=CODE)
    text(sl, 0.75, 4.8, 5.1, 0.9, 'Σ θ = 1 replaces one of the equations: the site balance makes one of them '
         'redundant.', size=14, color='muted')
    X, Y, Wd, Ht = 7.4, 1.95, 4.9, 2.9
    axes(sl, X, Y, Wd, Ht, 'one coverage θ (of 71)', 'net rate ṡ(θ)', ylabel_w=1.25)
    zero = Y + Ht * 0.6
    line(sl, X, zero, X + Wd, zero, color='rule', lw=1.0, dash=True)
    f = lambda t: 0.9 * math.exp(-3 * t) + 0.05          # noqa: E731  (0..1 for t in 0..1)
    fp = lambda t: -2.7 * math.exp(-3 * t)               # noqa: E731
    curve(sl, f, 0.0, 1.0, X, Y, Wd, Ht, color='H', lw=2.5)
    root = math.log(0.9 / 0.35) / 3
    rx = X + Wd * root
    t = 0.0
    for k in range(2):
        ft = f(t)
        t1 = t + (0.4 - ft) / fp(t)
        px, py = X + Wd * t, Y + Ht * (1 - ft)
        if k:
            poly(sl, [(px, zero), (px, py)], color='dark', lw=1.0, dash=True)
        poly(sl, [(px, py), (X + Wd * t1, zero)], color='dark', lw=1.75, end=True)
        t = t1
    circle(sl, rx, zero, 0.16, fill='ok', line=None)
    text(sl, rx + 0.12, zero - 0.42, 1.4, 0.3, 'steady state', size=13, color='ok', bold=True)
    text(sl, X + 1.4, Y + 0.1, 3.0, 0.3, 'Newton: follow the tangent', size=13, color='dark')
    for k in range(7):
        tt = 0.02 + (root - 0.02) * (1 - 0.62 ** k)
        circle(sl, X + Wd * tt, zero + 0.5, 0.09, fill='faint', line=None)
    text(sl, X + 0.15, zero + 0.62, 3.2, 0.3, 'time steps: smaller and smaller', size=13, color='muted')
    schematic(sl, X + Wd - 1.0, Y - 0.32)
    actions(sl, ['Unknowns: the coverage of each of the 71 surface species.',
                 'Equations: each species’ net production rate is zero.',
                 'Newton follows tangents: fast, but it can fail.',
                 'Time integration follows the real dynamics: slower, robust.'])
    notes(sl, ['docs/simulation_workflow.md section 2.4 (steady state: ṡ_k(θ) = 0 for every surface species, Σθ = 1,',
               f'one equation redundant) and 2.6 (Newton solve, stiffness). 71 = surface species in {RUN}/chem_annotated.yaml',
               '(including the empty site vacantX(8)).',
               'The one-dimensional picture is schematic: the real problem has 71 coupled unknowns.'])


def s29_solve(d):
    sl = d.slide('solve_steady_state(): choose a method by temperature, accept only checked results',
                 'Below 400 K the direct solve is skipped: there it tends to fail, and slowly.',
                 part=PG, tag='package · solver.py:146–165')
    step(sl, 0.55, 1.72, 2.2, 0.75, ['inputs', '`T, P, X, θ⁰`'], where='line 146')
    arrow(sl, 2.75, 2.1, 3.1, 2.1)
    gate(sl, 3.1, 1.62, 2.3, 0.95, 'T < 400 K?', size=14)
    ref(sl, 3.1, 2.57, 2.3, 'line 143', align='c')
    step(sl, 6.1, 1.62, 3.0, 0.42, 'yes: methods = (reactor,)', size=13)
    step(sl, 6.1, 2.15, 3.0, 0.42, 'no: methods = (direct, reactor)', size=13)
    arrow(sl, 5.4, 2.0, 6.1, 1.83)
    arrow(sl, 5.4, 2.2, 6.1, 2.36)
    shape(sl, 'rrect', 0.45, 3.0, 9.05, 1.95, fill=None, line='faint', lw=1.0, dash=True, radius=0.1)
    text(sl, 2.0, 3.04, 4.0, 0.3, 'for each method, in order', size=13, color='muted', italic=True)
    poly(sl, [(7.6, 2.57), (7.6, 2.8), (1.75, 2.8), (1.75, 3.45)], end=True)
    step(sl, 0.65, 3.45, 2.2, 0.85, ['`set_state()`', 'fresh surface · line 154'], size=13)
    arrow(sl, 2.85, 3.87, 3.2, 3.87)
    step(sl, 3.2, 3.45, 2.0, 0.85, ['run the method', 'line 155'], size=13)
    arrow(sl, 5.2, 3.87, 5.5, 3.87)
    gate(sl, 5.5, 3.3, 3.3, 1.15, 'finished, and check_convergence() passes? (156–159)', size=12)
    arrow(sl, 8.8, 3.87, 9.85, 3.87)
    yes_no(sl, 8.95, 3.55, 'yes')
    step(sl, 9.85, 3.45, 2.93, 0.85, ['return θ, info', 'info names the method'], fill='okbg', line='ok',
         where='lines 161–163')
    poly(sl, [(7.15, 4.45), (7.15, 4.78), (1.75, 4.78), (1.75, 4.3)], end=True, color='muted')
    text(sl, 2.0, 4.48, 4.5, 0.3, 'no: try the next method', size=13, color='muted')
    step(sl, 9.85, 4.65, 2.93, 0.6, 'none left: return None, failed', size=13, where='lines 164–165')
    arrow(sl, 9.5, 4.95, 9.85, 4.95, color='muted')
    cases = ['A method can finish but fail a check: that counts as a failure, and the next method runs (157–159).',
             'Every attempt, with its time and error, is kept in info["attempts"] (line 160).',
             'try_2/Cu is at 650 K or above, so direct runs first: it solved start A at all 72 conditions, start B at 71.']
    for i, t in enumerate(cases):
        x = 0.55 + i * 4.13
        shape(sl, 'rrect', x, 5.3, 3.97, 0.62, fill='white', line='rule', radius=0.08)
        text(sl, x + 0.1, 5.31, 3.77, 0.6, t, size=12, anchor='m')
    actions(sl, ['Below 400 K, use the reactor method only.',
                 'Start every attempt from a freshly set state.',
                 'Accept a result only if check_convergence() passes.',
                 'Otherwise try the next method, or report failure.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py: default_methods() lines 142-143, solve_steady_state() lines',
               '146-165 (unchanged since the first commit). direct_min_T = 400 K (config.py:21, SOLVER_DEFAULTS).',
               f'{RUN}/run.log, main steady states: 71 conditions "A direct, B direct" and 1 "A direct, B reactor":',
               'the reactor method solved the one start B the direct solve did not.',
               'Doc Stage 4, "Choosing the method" and Diagram 3.'])


def s30_set_state(d):
    sl = d.slide('set_state(): a fresh surface, and the gas set but never integrated',
                 'Cantera’s steady-state solver keeps state on the surface object, so every attempt gets a new one.',
                 part=PG, tag='package · solver.py:31–52')
    heads = ['1 · fresh surface', '2 · set the gas', '3 · check the start', '4 · set and confirm']
    refs_ = ['line 36', 'lines 37–38', 'lines 39–44', 'lines 45–52']
    xs0 = [0.55 + i * (2.94 + 0.15) for i in range(4)]
    for i, x in enumerate(xs0):
        card(sl, x, 1.66, 2.94, 3.25, heads[i], size=15)
        ref(sl, x + 0.15, 4.58, 2.6, refs_[i])
        if i < 3:
            arrow(sl, x + 2.94, 3.25, x + 3.09, 3.25, color='ink', lw=1.5)
    x = xs0[0]
    shape(sl, 'rrect', x + 0.3, 2.25, 2.34, 1.25, fill='white', line='faint', dash=True, radius=0.08)
    text(sl, x + 0.3, 2.3, 2.34, 0.3, 'new Interface object', size=13, color='muted', align='c')
    for k in range(4):
        site(sl, x + 0.62 + k * 0.5, 3.1, d=0.36)
    text(sl, x + 0.15, 3.65, 2.65, 0.85, '`model.new_surface()`', size=13)
    x = xs0[1]
    shape(sl, 'rrect', x + 0.3, 2.2, 2.34, 1.35, fill='panel2', line=None, radius=0.08)
    h2(sl, x + 0.85, 2.6, d=0.3)
    h2(sl, x + 1.95, 2.75, d=0.3)
    gas(sl, x + 1.35, 3.15, 'CO', w=0.55, h=0.3, size=12)
    text(sl, x + 0.15, 3.65, 2.65, 0.85, '`gas.TPX = T, P, X`', size=13)
    x = xs0[2]
    for k, (t, y) in enumerate([('names known?', 2.18), ('sums to 1?', 3.0)]):
        shape(sl, 'diamond', x + 0.2, y, 1.75, 0.72, fill='white', line='dark', lw=1.25, content=t, size=12,
              margin=0.0)
        text(sl, x + 2.0, y + 0.2, 0.9, 0.35, 'no: error', size=12, color='ink', bold=True)
    text(sl, x + 0.15, 3.78, 2.7, 0.7, 'np.isclose: within about 10^{−5}', size=12, color='muted')
    x = xs0[3]
    slab(sl, x + 0.37, 2.55, 4, pitch=0.55, d=0.38, metal_label='')
    text(sl, x + 0.15, 3.1, 2.7, 0.3, 'clean start: θ_{vacantX} = 1', size=13)
    text(sl, x + 0.15, 3.45, 2.7, 1.1, ['assert T, P, X applied', 'return a copy of gas X'], size=13)
    cases = ['θ⁰ = {vacantX(8): 0.99, [PtH](33): 0.02} sums to 1.01: ValueError. Cantera would otherwise rescale it '
             'and run a different start.',
             'Why a fresh object: on one reused Interface, a second steady-state solve failed with a singular-matrix '
             'error (Cantera 3.0.1).']
    for i, t in enumerate(cases):
        x = 0.55 + i * 6.17
        shape(sl, 'rrect', x, 5.05, 6.03, 0.85, fill='white', line='rule', radius=0.08)
        text(sl, x + 0.1, 5.06, 5.83, 0.83, t, size=14, anchor='m')
    actions(sl, ['Make a fresh surface object for this attempt.',
                 'Set the gas temperature, pressure and composition.',
                 'Refuse starting coverages that are unknown or don’t sum to 1.',
                 'Return the gas composition, to check later that it never moved.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py set_state() lines 31-52; Model.new_surface() model.py lines 185-210',
               '(also applies optional rate multipliers and enthalpy shifts for Stages 10-11).',
               'The sum check (line 43) is np.isclose(total, 1.0, atol=1e-8) with NumPy\'s default rtol = 1e-5, so sums',
               'within about 1e-5 of 1 pass; 1.01 fails. (The docs\' "1 ± 1e-8" is tighter than the code.)',
               'Doc Stage 4 "Setting the state" and the Important box "A fresh surface object for every solve".',
               'The θ⁰ example in the edge case is illustrative input, not from a run.'])


def s31_direct(d):
    sl = d.slide('Method 1, direct: Cantera’s Newton solve, fast when it works',
                 'advance_coverages_to_steady_state() solves the rate equations directly, with damping and time steps.',
                 part=PG, tag='package · solver.py:88–95')
    X, Y, Wd, Ht = 1.9, 1.95, 4.0, 1.75
    axes(sl, X, Y, Wd, Ht, 'θ', 'ṡ(θ)', ylabel_w=1.1)
    zero = Y + Ht * 0.5
    line(sl, X, zero, X + Wd, zero, color='rule', dash=True)
    f = lambda t: 0.95 - 0.9 * t ** 0.6                 # noqa: E731
    curve(sl, f, 0.0, 1.0, X, Y, Wd, Ht, color='H', lw=2.5)
    t = 0.05
    for k in range(2):
        ft = f(t)
        slope = -0.9 * 0.6 * t ** (-0.4)
        t1 = min(0.98, t + (0.5 - ft) / slope)
        poly(sl, [(X + Wd * t, zero), (X + Wd * t, Y + Ht * (1 - ft))], color='dark', lw=1.0, dash=True)
        poly(sl, [(X + Wd * t, Y + Ht * (1 - ft)), (X + Wd * t1, zero)], color='dark', lw=1.75, end=True)
        t = t1
    schematic(sl, X + Wd - 0.9, Y - 0.32)
    cases = [('Singular matrix', 'a reused surface object, or species switched off by multipliers: their Jacobian '
              'rows are all zero.'),
             ('Slow failure at low T', 'expected below about 400 K, so solve_steady_state() skips it there.')]
    for i, (h, t) in enumerate(cases):
        y = 1.72 + i * 1.12
        card(sl, 6.6, y, 6.18, 1.0, None, fill='white', line='rule')
        text(sl, 6.75, y + 0.06, 5.9, 0.3, h, size=14, bold=True)
        text(sl, 6.75, y + 0.38, 5.9, 0.6, t, size=14, color='muted')
    code(sl, 0.55, 4.05, 7.35, 1.82, ['t0 = time.time()', 'try:',
                                       '    model.surf.advance_coverages_to_steady_state()',
                                       '    return True, None, time.time() - t0, None',
                                       'except ct.CanteraError as e:',
                                       '    return False, _cantera_error(e), time.time() - t0, None'],
         size=13, first=90)
    text(sl, 8.15, 4.15, 4.6, 1.6, ['The gas is never touched.', 'A Cantera error comes back as a short message '
                                     '(line 95), so one bad condition cannot stop a whole batch.'], size=14,
         color='muted', space_after=6)
    actions(sl, ['Linearise the rate equations at the current coverages.',
                 'Step to where the linear model is zero.',
                 'Repeat until the change is tiny.',
                 'On a Cantera error, return the message instead of crashing.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py method_direct() lines 88-95. Cantera calls',
               'solvePseudoSteadyStateProblem: damped Newton with time stepping (doc section 2.6).',
               'Edge cases from doc Stage 4 checks table and section 10.2 (switched-off species).',
               'The Newton picture is schematic and one-dimensional.'])


def s32_reactor(d):
    sl = d.slide('Method 3, reactor: age the surface in a gas so large it cannot change',
                 '1 m³ of gas over 10⁻¹⁰ m² of catalyst, gas reactions off, integrated in half-decade steps.',
                 part=PG, tag='package · solver.py:98–136')
    shape(sl, 'rrect', 0.75, 1.75, 5.0, 3.55, fill='panel', line='faint', lw=1.25, radius=0.08)
    text(sl, 0.9, 1.82, 4.0, 0.3, '1 m³ of gas, energy off: T fixed', size=14, bold=True)
    for (px, py) in [(1.4, 2.5), (2.6, 2.75), (3.9, 2.45), (1.9, 3.4), (3.3, 3.55), (4.6, 3.2), (1.3, 4.25),
                     (2.8, 4.3), (4.3, 4.2)]:
        h2(sl, px, py, d=0.28)
    shape(sl, 'rect', 3.1, 5.12, 0.3, 0.1, fill='dark', line=None)
    text(sl, 3.55, 4.95, 2.1, 0.32, '10^{−10} m² catalyst', size=13)
    shape(sl, 'rrect', 4.1, 1.78, 1.55, 0.38, fill='dark', line=None, radius='pill',
          content='gas reactions × 0', size=12, color='white', bold=True)
    text(sl, 0.75, 5.35, 5.0, 0.5, 'A/V = 10^{−10} m^{−1}: the surface cannot change the gas.', size=14,
         color='muted')
    X, Y, Wd = 6.4, 2.6, 6.2
    arrow(sl, X, Y, X + Wd, Y, color='ink', lw=1.5)
    for k in range(33):
        x = X + Wd * k / 33
        tall = k % 2 == 0
        line(sl, x, Y - (0.16 if tall else 0.09), x, Y + (0.16 if tall else 0.09), color='dark', lw=1.0)
    for k, lab in [(0, '10^{−8}'), (8, '10^{−4}'), (16, '10^{0}'), (24, '10^{4}'), (32, '10^{8} s')]:
        text(sl, X + Wd * k / 33 - 0.4, Y + 0.22, 0.8, 0.3, lab, size=13, align='c')
    text(sl, X, Y - 0.75, 6.2, 0.35, '33 checkpoints, half a decade apart (line 121)', size=14, bold=True)
    rows = ['IdealGasReactor, volume 1 m³, energy off (lines 102–103)',
            'ReactorSurface with area A = V × 10^{−10} m^{−1} (line 104)',
            'tolerances: rtol 10^{−10}, atol 10^{−22} (line 106)',
            'gas reactions back on afterwards, even after an error (lines 135–136)']
    for i, t in enumerate(rows):
        text(sl, X, 3.35 + i * 0.5, 6.3, 0.45, t, size=14)
    actions(sl, ['Put the surface in 1 m³ of gas with 10⁻¹⁰ m² of catalyst.',
                 'Switch gas reactions off, so only the surface reacts.',
                 'Integrate through 33 checkpoints, 10⁻⁸ to 10⁸ s.',
                 'Switch gas reactions back on afterwards, even after an error.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py fixed_gas_reactor() lines 98-107 and method_reactor() lines 110-136.',
               'n_checkpoints = round(2 (log10 1e8 + 8)) + 1 = 33. Settings in config.py SOLVER_DEFAULTS (lines 22-25):',
               'reactor_area_per_volume 1e-10, reactor_rtol 1e-10, reactor_atol 1e-22, reactor_t_max 1e8 s.',
               '"Method 2" (repeated advance_coverages) was dropped as too slow; numbering kept (doc Stage 4).'])


def s33_plateau(d):
    sl = d.slide('The plateau test: stop only when nothing changed since the last checkpoint',
                 'A small residual at one moment is not enough: the surface may still be drifting slowly.',
                 part=PG, tag='package · solver.py:122–132')
    X, Y, Wd, Ht = 1.75, 1.95, 4.9, 3.1
    axes(sl, X, Y, Wd, Ht, 'time (log), one dot per checkpoint', 'forward TOF', ylabel_w=1.1)
    vals = [0.05, 0.12, 0.28, 0.5, 0.66, 0.74, 0.78, 0.8, 0.81]
    for i, v in enumerate(vals):
        cx, cy = X + 0.25 + i * 0.5, Y + Ht * (1 - v)
        if i < 4:
            circle(sl, cx, cy, 0.2, fill='white', line='dark', lw=1.25)
        elif i < 8:
            circle(sl, cx, cy, 0.2, fill='amber', line=None)
        else:
            circle(sl, cx, cy, 0.26, fill='ok', line=None)
    for i in range(len(vals) - 1):
        a, b = vals[i], vals[i + 1]
        line(sl, X + 0.25 + i * 0.5, Y + Ht * (1 - a), X + 0.25 + (i + 1) * 0.5, Y + Ht * (1 - b), color='faint')
    leg = [('white', 'dark', 'residual too large'), ('amber', None, 'residual ok, still changing'),
           ('ok', None, 'accepted: ok and unchanged')]
    for i, (fill, ln, t) in enumerate(leg):
        cx = 0.75 + i * 2.05
        circle(sl, cx, 5.6, 0.2, fill=fill, line=ln, lw=1.25)
        text(sl, cx + 0.16, 5.44, 1.9, 0.35, t, size=12, color='muted')
    schematic(sl, X + Wd - 0.9, Y - 0.32)
    step(sl, 7.3, 1.72, 2.6, 0.55, ['advance to the next checkpoint', '(line 123)'], size=12)
    arrow(sl, 8.6, 2.27, 8.6, 2.52)
    step(sl, 7.3, 2.52, 2.6, 0.55, ['`check_convergence()`', '(line 124)'], size=12)
    arrow(sl, 8.6, 3.07, 8.6, 3.3)
    gate(sl, 7.35, 3.3, 2.5, 0.95, 'residual ok? (126)', size=13)
    arrow(sl, 8.6, 4.25, 8.6, 4.42)
    gate(sl, 6.95, 4.42, 3.3, 1.15, 'TOF and θ changed < 10^{−4} since the last? (127–129)', size=12)
    arrow(sl, 10.25, 5.0, 11.0, 5.0, color='ok')
    yes_no(sl, 10.35, 4.66, 'yes')
    step(sl, 11.0, 4.7, 1.78, 0.6, 'accept, return t', size=13, fill='okbg', line='ok', where='line 130')
    poly(sl, [(7.35, 3.77), (6.95, 3.77), (6.95, 2.0), (7.3, 2.0)], end=True, color='muted')
    poly(sl, [(6.95, 5.0), (6.95, 3.77)], color='muted')
    yes_no(sl, 6.98, 3.45, 'no', w=0.4)
    yes_no(sl, 6.98, 4.65, 'no', w=0.4)
    text(sl, 10.3, 1.75, 2.5, 0.9, 'after 10^{8} s: “no plateau by t\u00a0=\u00a01e+08\u00a0s” (line 132)', size=13,
         color='ink')
    actions(sl, ['At each checkpoint, test the residuals.',
                 'Compare TOF and coverages with the previous checkpoint.',
                 'Accept when both changed by less than 10⁻⁴.',
                 'If 10⁸ s comes first, report “no plateau”.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py method_reactor() lines 120-136; plateau_tof_rtol = 1e-4 and',
               'plateau_theta_atol = 1e-4 (config.py:26-27). Doc Stage 4: plateau criterion and Figure 4 (fig_convergence.png).',
               'Design choice in the doc: a small residual can occur while the surface still drifts slowly.',
               'The checkpoint plot is schematic. Planned in the doc: a separate "plateau" label for high-T cases',
               'where fast species keep residuals above the limit because of integrator error.'])


def s34_check(d):
    sl = d.slide('check_convergence(): four tests every steady state must pass',
                 'Each species is judged against its own turnover, so fast and slow species are both held to account.',
                 part=PG, tag='package · solver.py:55–85')
    tests = [('residual', '|dθ_{k}/dt| ≤ 10^{−8} + 10^{−12} G_{k}', 'for every species k', 'lines 59–62, 70'),
             ('site balance', '|Σ θ − 1| ≤ 10^{−8}', 'coverages add to 1', 'lines 66, 71'),
             ('no negatives', 'min θ ≥ −10^{−12}', 'no coverage below zero', 'lines 67, 72'),
             ('gas fixed', 'max |x − x⁰| ≤ 10^{−9}', 'the gas did not move', 'lines 68, 73')]
    for i, (h, f, sub, r) in enumerate(tests):
        x = 0.55 + i * 2.6
        card(sl, x, 1.72, 2.4, 2.1, None, fill='white', line='dark', lw=1.25)
        tick(sl, x + 0.3, 2.0, d=0.32)
        text(sl, x + 0.55, 1.83, 1.8, 0.35, h, size=15, bold=True)
        text(sl, x + 0.12, 2.35, 2.2, 0.75, f, size=13, font=MATH, align='c', anchor='m')   # each on one line
        text(sl, x + 0.12, 3.12, 2.2, 0.32, sub, size=13, color='muted', align='c')
        ref(sl, x + 0.12, 3.47, 2.2, r, align='c')
        if i < 3:
            text(sl, x + 2.4, 2.55, 0.2, 0.4, '+', size=16, bold=True, align='c')
    arrow(sl, 10.95, 2.77, 11.25, 2.77, color='ok', lw=2)
    step(sl, 11.25, 2.35, 1.53, 0.85, ['converged', '= all four'], fill='okbg', line='ok', where='line 76')
    card(sl, 0.55, 4.05, 7.7, 1.85, 'Why both an absolute and a relative term')
    text(sl, 0.7, 4.5, 7.4, 1.35, 'A species turning over 10⁹ times per second carries round-off near 10⁻⁷ s⁻¹ '
         '(2.2 × 10⁻¹⁶ × G): an absolute limit of 10⁻⁸ alone could never pass. A relative limit alone would '
         'accept a large residual on a slow species.', size=14)
    card(sl, 8.5, 4.05, 4.28, 1.85, 'It also records')
    text(sl, 8.65, 4.5, 4.0, 1.35, 'the species closest to its limit, its residual and allowed value, and the forward '
         'TOF (lines 63–65, 78–81).', size=14)
    actions(sl, ['Every species’ net rate must be tiny against its own turnover.',
                 'The coverages must still add to 1.',
                 'No coverage may be negative.',
                 'The gas must not have moved.'])
    notes(sl, ['Hydrogen-Surface-kinetics solver.py check_convergence() lines 55-85. Tolerances in config.py',
               'SOLVER_DEFAULTS (lines 16-20): abs_tol 1e-8 1/s, rel_gross_tol 1e-12, sum_tol 1e-8, neg_tol -1e-12,',
               'gas_tol 1e-9. G_k = σ_k (creation + destruction)/Γ (lines 60-61). Doc Stage 4 "Convergence test" and',
               'section 2.6 (residual scaling). A test relative to the probe TOF was tried and rejected (doc design choices).'])


def _panel(sl, x, y, w, h, title, ylab, curves, accept=0.62, limit=None):
    """A small schematic trace: curves are (f, colour, width) with f: 0..1 -> 0..1 (0 = bottom)."""
    text(sl, x - 0.1, y - 0.58, w + 0.2, 0.32, title, size=14, bold=True)
    axes(sl, x, y, w, h, '', '')
    text(sl, x, y + h + 0.04, w, 0.28, 'time (log)', size=12, color='muted', align='c')
    text(sl, x + 0.08, y - 0.27, w, 0.26, ylab, size=12, color='muted')       # above the plot, clear of the curves
    if limit is not None:
        line(sl, x, y + h * (1 - limit), x + w, y + h * (1 - limit), color='muted', lw=1.0, dash=True)
    for f, col, lw_ in curves:
        curve(sl, f, 0.0, 1.0, x, y, w, h, color=col, lw=lw_)
    line(sl, x + w * accept, y, x + w * accept, y + h, color='ok', lw=1.25, dash=True)


def s35_traces(d):
    sl = d.slide('Convergence traces: follow chosen conditions to 10⁸ s to check the state is final',
                 'The same reactor, tolerances and checkpoints as the reactor method, but every checkpoint is kept.',
                 part=PG, tag='package · analyses/convergence.py')
    sig = lambda t, c, k: 1 / (1 + math.exp(-k * (t - c)))                      # noqa: E731
    _panel(sl, 0.95, 2.25, 3.3, 1.45, '(a) coverages', 'θ',
           [(lambda t: 0.08 + 0.82 * sig(t, 0.35, 14), 'H', 2.25), (lambda t: 0.92 - 0.8 * sig(t, 0.35, 14), 'faint', 2.0)])
    _panel(sl, 5.05, 2.25, 3.3, 1.45, '(b) forward TOF', 'TOF_{f}',
           [(lambda t: 0.1 + 0.75 * sig(t, 0.3, 12), 'dark', 2.25)])
    _panel(sl, 9.15, 2.25, 3.3, 1.45, '(c) normalized residual (log)', 'residual / allowed',
           [(lambda t: max(0.05, 0.95 - 1.1 * t), 'dark', 2.25)], limit=0.35)
    text(sl, 11.45, 2.25 + 1.45 * 0.65 - 0.3, 1.0, 0.28, 'limit = 1', size=12, color='muted')
    text(sl, 0.95, 4.05, 11.5, 0.3, 'dashed green: the checkpoint where the reactor method accepted the state',
         size=13, color='ok')
    schematic(sl, 11.0, 4.05)
    card(sl, 0.55, 4.45, 5.95, 1.45, 'Which conditions')
    text(sl, 0.7, 4.85, 5.65, 1.0, 'The main conditions whose clean start failed (at most 6); if all converged, the '
         'reference at the first point. For try_2/Cu: H₂ only, 650 K, 1 bar.', size=13)
    card(sl, 6.83, 4.45, 5.95, 1.45, 'What it writes; a good trace')
    text(sl, 6.98, 4.85, 5.65, 1.0, 'convergence.csv (one row per checkpoint) and convergence.png. Good: everything '
         'flat after the dashed line, the residual below 1 to the end.', size=13)
    actions(sl, ['Pick the conditions: failures first, else the reference.',
                 'Integrate from a clean surface exactly as the reactor method does.',
                 'Record every checkpoint, on to 10⁸ s.',
                 'Check that nothing moves after the accepted checkpoint.'])
    notes(sl, ['analyses/convergence.py: docstring (lines 1-13); checkpoints() (26-30): per_decade x (log10 t_max + 8) + 1;',
               'trace() (32), default_conditions() (65), run() (75), plot_convergence() (112).',
               "Settings: config.py convergence (lines 77-82): conditions None (the failed main conditions, or the reference",
               'at the first point), max_conditions 6, per_decade 2 (the reactor method\'s own checkpoints).',
               'convergence.csv columns: t_s, tof_fwd, norm_residual, worst_species, residual_ok, plateau_ok, accepted_here,',
               'case, T_K, P_ref_bar, accepted_s, error.',
               f'{RUN}: every main condition converged, so the default traces "H2 only" at 650 K, 1 bar (output/convergence.csv).',
               'Docs: Stage 4 "On real data" and Table R2 (convergence). The three panels here are schematic.'])


SLIDES = [s28_problem, s29_solve, s30_set_state, s31_direct, s32_reactor, s33_plateau, s34_check, s35_traces]
