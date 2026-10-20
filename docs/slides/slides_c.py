"""Section C (slides 11-18): the theory Cantera relies on, with the numbers of runs/try_2/Cu:
two phases, coverages, site counting, mass action, sticking, reverse rates, steady state, and
rates per site and per area."""
import math

from deck_lib import *  # noqa: F401,F403

PC = 'C · Into Cantera: the theory'
RUN = 'runs/try_2/Cu'


def s11_two_phases(d):
    sl = d.slide('Cantera loads two phases: the gas above and the sites below',
                 'Surface reactions connect them: adsorption takes gas, desorption gives it back.',
                 part=PC, tag='package · model.py')
    shape(sl, 'rrect', 0.55, 1.68, 7.6, 2.0, fill='panel', line=None, radius=0.1)
    text(sl, 0.7, 1.74, 6.0, 0.32, '**gas phase**: T, P and mole fractions x', size=14, color='ink')
    for cx, cy in [(1.4, 2.55), (3.2, 2.85), (5.1, 2.45), (6.95, 2.95)]:
        h2(sl, cx, cy)
    gas(sl, 2.3, 3.3, 'CO')
    gas(sl, 6.0, 3.35, 'H2O', 'H₂O')
    xs = slab(sl, 0.7, 4.75, 12, occ={2: ('H', 'H', 1), 3: ('H', 'H', 1), 7: ('CO', 'CO', 1), 9: ('H', 'H', 1)})
    text(sl, 0.7, 5.36, 7.4, 0.32, '**surface phase**: sites, each empty or held, with coverages θ', size=14)
    cx = (xs[2] + xs[3]) / 2
    arrow(sl, 3.2, 3.05, cx, 4.38, color='H', lw=2)
    text(sl, 1.45, 3.75, 1.4, 0.3, 'adsorption', size=14, color='H')
    arrow(sl, xs[9], 4.38, 6.95, 3.15, color='H', lw=2)
    text(sl, 7.0, 3.85, 1.2, 0.3, 'desorption', size=14, color='H')
    code(sl, 8.55, 1.72, 4.23, 1.5, ['self.gas = ct.Solution(path,', '    self.cfg.gas_phase,', '    transport_model=None)',
                                      'self.surf = ct.Interface(path,', '    self.cfg.surface_phase,', '    [self.gas])'],
         size=13)
    ref(sl, 8.55, 3.25, 4.2, 'model.py:166–167')
    text(sl, 8.55, 3.68, 4.23, 0.62, '**Interface**: the surface phase plus its reactions with the gas.', size=14)
    text(sl, 8.55, 4.38, 4.23, 0.62, '**transport_model=None**: a uniform gas needs no transport data.', size=14)
    text(sl, 8.55, 5.08, 4.23, 0.62, 'Phase names: "gas" and "surface", the settings defaults.', size=14, color='muted')
    actions(sl, ['Load the gas phase: species, T, P and composition.',
                 'Load the surface phase, with the gas as its neighbour.',
                 'Surface reactions move species between the two phases.',
                 'Skip transport: nothing flows through space here.'])
    notes(sl, ['Hydrogen-Surface-kinetics/hydrogen_surface_kinetics/model.py: Model._load_phases, lines 165-167, called from',
               '_load_mechanism (lines 137-163) inside the warnings catcher of __init__ (lines 91-94).',
               'Phase names from config.py defaults (gas_phase "gas", surface_phase "surface"), as in the try_2/Cu YAML.',
               'If Cantera refuses the file only because of unmatched duplicate flags, a repaired copy is loaded (Stage 0).',
               'Doc Stage 1 procedure step 1 and design choice "Transport off". Picture is schematic.'])


_OCC10 = {(0, 1): ('H', 'H', 1), (0, 3): ('N2', 'N₂', 2), (1, 0): ('H', 'H', 1), (1, 1): ('N2', 'N₂', 2),
          (1, 3): ('H', 'H', 1)}                     # 3 empty, 3 H*, 2 N2** (4 sites): 10 sites


def s12_coverage(d):
    sl = d.slide('A coverage is the share of sites a species holds; the shares add to 1',
                 'Ten sites, counted one by one. A species that sits on two sites counts twice.',
                 part=PC, tag='concept')
    grid(sl, 1.0, 2.25, 5, 2, pitch=0.9, d=0.72, occ=_OCC10)
    text(sl, 0.9, 4.05, 5.2, 0.3, 'example surface with 10 sites (not a result)', size=14, color='muted', italic=True)
    rows = [('site', 'empty sites', '3 sites', 'θ_{∗} = 0.3'),
            ('H', 'H∗', '3 sites', 'θ_{H∗} = 0.3'),
            ('N2', 'N₂∗∗', '2 molecules × 2 sites = 4', 'θ_{N₂∗∗} = 0.4')]
    y = 1.95
    for icon, name, count, th in rows:
        if icon == 'site':
            site(sl, 6.95, y + 0.25, d=0.42)
        elif icon == 'H':
            adsorbate(sl, 6.95, y + 0.25, 'H', 'H', d=0.42)
        else:
            adsorbate(sl, 7.1, y + 0.25, 'N2', 'N₂', d=0.42, n_sites=2, pitch=0.38, size=12)
        text(sl, 7.55, y + 0.05, 1.3, 0.4, name, size=16, bold=True)
        text(sl, 8.8, y + 0.05, 2.5, 0.4, count, size=16)
        text(sl, 11.2, y + 0.02, 1.6, 0.45, th, size=18, font=MATH)
        y += 0.66
    line(sl, 6.7, y + 0.05, 12.75, y + 0.05, color='faint', lw=1.0)
    text(sl, 7.55, y + 0.12, 1.3, 0.4, 'total', size=16, bold=True)
    text(sl, 8.8, y + 0.12, 2.5, 0.4, '10 sites', size=16)
    text(sl, 11.2, y + 0.1, 1.6, 0.45, 'Σ θ = 1', size=18, font=MATH)
    text(sl, 0.55, 4.75, W - 2 * LM, 0.5,
         'θ_{k} = (sites held by species k) / (all sites)          Σ_{k} θ_{k} = 1, empty sites included',
         size=20, font=MATH, align='c')
    actions(sl, ['Count the sites each species holds.',
                 'Divide by the total number of sites.',
                 'Empty sites count too, so the shares add to 1.',
                 'A two-site species such as N₂∗∗ counts twice.'])
    notes(sl, ['docs/simulation_workflow.md section 2.1 and Figure 2 (ten sites hold eight entities).',
               'Each species k occupies σ_k sites ("sites:" in the YAML, default 1). N₂∗∗ here stands for',
               f'[Pt]N=N[Pt](57) ("sites: 2"), one of the 18 two-site species of {RUN}/chem_annotated.yaml.',
               'The surface is an example, not a result.'])


def s13_site_factor(d):
    sl = d.slide('Cantera counts entities, so two-site species make the surface look bigger',
                 'The same ten sites hold 8 entities. Cantera sets those 8 to Γ, which implies 1.25 Γ sites.',
                 part=PC, tag='concept')
    pos = grid(sl, 0.85, 2.2, 5, 2, pitch=0.78, d=0.62, occ=_OCC10)
    ents = [(0, 0), (0, 1), (0, 2), (0, 3), (1, 0), (1, 1), (1, 3), (1, 4)]
    for n, key in enumerate(ents, 1):
        cx, cy = pos[key]
        if key in ((0, 3), (1, 1)):
            cx += 0.39
        badge(sl, cx + 0.3, cy - 0.32, n, d=0.28)
    text(sl, 0.75, 3.9, 4.6, 0.6, '8 entities: 3 empty sites, 3 H∗, 2 N₂∗∗', size=14, color='muted')
    eqs = ['S = Σ_{k} θ_{k} / σ_{k}', '   = 0.3 + 0.3 + 0.4 / 2 = 0.8', 'f = 1 / S = 1.25', 'C_{∗} = f θ_{∗} Γ']
    for i, e in enumerate(eqs):
        text(sl, 5.45, 1.95 + i * 0.55, 3.6, 0.5, e, size=20, font=MATH)
    card(sl, 9.25, 1.8, 3.53, 2.7, 'What it changes')
    text(sl, 9.4, 2.25, 3.25, 1.3, 'A step that needs n empty sites runs f^{n} times faster than θ_{∗} alone '
         'suggests.', size=14)
    text(sl, 9.4, 3.45, 3.25, 0.9, 'R0 needs n = 2:  1.25² ≈ 1.56 ×', size=16, bold=True)
    caveat(sl, 0.55, 4.85, W - 2 * LM, 0.82, '18 of the 71 surface species here take two sites; where they dominate, '
           'f approaches 2. The workflow reports f for every condition (site_count_factor, Stage 7), and Stage 9 shows '
           'a fixed-site-count estimate.')
    actions(sl, ['Divide each coverage by its species’ site count σ.',
                 'Cantera holds these entities at Γ in total.',
                 'The implied number of sites is f Γ, here 1.25 Γ.',
                 'Rates that need n empty sites pick up f^{n}.'])
    notes(sl, ['docs/simulation_workflow.md section 2.1: C_k = Γ θ_k / (σ_k S), S = Σ θ_j / σ_j, f = 1/S,',
               'C_∗ = f θ_∗ Γ; verified against Cantera in Appendix E (docs/figures/check_equations.py).',
               'Consequences box of section 2.1; site_count_factor in Stage 7 (analyses/common.py); open markers in Stage 9.',
               f'{RUN}: 53 one-site and 18 two-site surface species (counted with a YAML parser).',
               'The convention could change between Cantera versions: re-run check_equations.py after upgrading.'])


def s14_mass_action(d):
    sl = d.slide('Rates follow mass action: a rate constant times concentrations',
                 'For R0 the forward rate needs one H₂ and two empty sites; the reverse needs two H∗.',
                 part=PC, tag='concept')
    xl = slab(sl, 0.85, 3.45, 2, metal_label='')
    h2(sl, (xl[0] + xl[1]) / 2, 2.5)
    slab(sl, 4.95, 3.45, 2, occ={0: ('H', 'H', 1), 1: ('H', 'H', 1)}, metal_label='')
    arrow(sl, 2.35, 2.95, 4.75, 2.95, color='ink', lw=2.5)
    arrow(sl, 4.75, 3.55, 2.35, 3.55, color='ink', lw=2.5)
    text(sl, 2.1, 2.4, 2.9, 0.4, 'q_{f} = k_{f} C_{H₂} C_{∗}²', size=18, font=MATH, align='c')
    text(sl, 2.1, 3.65, 2.9, 0.4, 'q_{r} = k_{r} C_{H∗}²', size=18, font=MATH, align='c')
    text(sl, 0.75, 4.45, 5.6, 0.45, 'net rate of progress:  q = q_{f} − q_{r}', size=20, font=MATH)
    text(sl, 0.75, 5.0, 5.6, 0.4, 'units: kmol per m² of surface per second', size=14, color='muted')
    card(sl, 6.83, 1.66, 5.95, 4.15, 'Every reaction i')
    eqs = [('q_{i} = k_{f,i} Π_{k} C_{k}^{ν′} − k_{r,i} Π_{k} C_{k}^{ν″}', 20),
           ('k_{f} = A T^{b} exp(−E_{a} / RT)', 20),
           ('surface:  C_{k} = Γ θ_{k} / (σ_{k} S)   kmol/m²', 18),
           ('gas:  C_{k} = x_{k} P / RT   kmol/m³', 18)]
    for i, (e, s) in enumerate(eqs):
        text(sl, 7.0, 2.2 + i * 0.75, 5.7, 0.6, e, size=s, font=MATH)
    text(sl, 7.0, 5.15, 5.6, 0.5, 'ν′ and ν″: how many of species k react and form.', size=14, color='muted')
    actions(sl, ['Multiply the forward rate constant by each reactant’s concentration.',
                 'Surface concentrations come from coverages, gas ones from pressure.',
                 'Do the same for the products, with the reverse constant.',
                 'The net rate is forward minus reverse.'])
    notes(sl, ['docs/simulation_workflow.md section 2.2 (rate of progress, Arrhenius form) and 2.1 (concentrations).',
               f'R0 as Cantera loads it from {RUN}: H2(1) + 2 vacantX(8) <=> 2 [PtH](33). Pictures are schematic.'])


def s15_sticking(d):
    sl = d.slide('Adsorption uses a sticking probability: about 1 % of H₂ hits split on bare Cu at 650 K',
                 'At 650 K and 1 bar, each site is hit by H₂ about 4.1 × 10⁸ times per second.',
                 part=PC, tag='chem_annotated.yaml')
    step(sl, 0.75, 1.8, 5.5, 0.95, ['H₂ collisions per site per second (kinetic theory)',
                                     ('**4.1 × 10^{8} s^{−1}**', {'size': 18})])
    arrow(sl, 3.5, 2.75, 3.5, 3.12)
    shape(sl, 'rrect', 2.35, 3.12, 2.3, 0.5, fill='dark', line=None, radius='pill', content='× γ = 0.0119',
          size=16, color='white', bold=True)
    arrow(sl, 3.5, 3.62, 3.5, 3.98)
    step(sl, 0.75, 3.98, 5.5, 0.95, ['splitting events per site per second (R0, bare surface)',
                                      ('**4.9 × 10^{6} s^{−1}**', {'size': 18})], line='H', lw=1.5)
    arrow(sl, 4.65, 3.37, 5.7, 3.37, color='faint')
    text(sl, 5.15, 3.08, 1.6, 0.3, '98.8 % bounce', size=13, color='muted')
    card(sl, 6.83, 1.66, 5.95, 3.3, 'From γ to a rate constant')
    eqs = [('γ = a T^{b} exp(−c / RT)', 20), ('k_{f} = (γ / Γ^{m}) √(RT / 2πW)', 20),
           ('m = number of surface reactants: 2 for R0', 14), ('W = molar mass of the gas: 2.016 kg/kmol for H₂', 14),
           ('R0: a = 0.03, b = 0, c = 5 kJ/mol, so γ grows with T:', 14),
           ('0.0119 at 650 K, 0.0182 at 1200 K', 14),
           ('R1, H₂ + site ⇌ H₂∗: a = 0.008 (vdW template)', 14)]
    yy = 2.12
    for e, s in eqs:
        text(sl, 7.0, yy, 5.7, 0.5, e, size=s, font=MATH if s > 14 else BODY)
        yy += 0.5 if s > 14 else 0.36
    caveat(sl, 0.55, 5.1, W - 2 * LM, 0.78, 'Cantera warns when γ exceeds 1 at some temperature: 6 reactions in this '
           'mechanism do, e.g. N₂ + 2 ∗ ⇌ N₂∗∗ with a = 1.66. The package records these warnings at load (model.py:91–94).')
    actions(sl, ['Kinetic theory gives how often H₂ hits each site.',
                 'A fraction γ of hits on a bare surface reacts.',
                 'γ becomes a rate constant through Γ and the molar mass W.',
                 'Cantera flags any γ above 1 when it loads the file.'])
    notes(sl, ['Impingement Z = p / sqrt(2 pi M R T) = 1.209e4 mol m-2 s-1 at 650 K, 1 bar (M = 2.016e-3 kg/mol);',
               'per site Z / Γ with Γ = 2.943e-5 mol/m2: 4.11e8 per second. γ(650 K) = 0.03 exp(-5000/(R 650)) = 0.0119,',
               'so bare-surface TOF_f = γ Z / Γ = 4.88e6 per second (θ_∗ = 1, f = 1). At 1200 K: γ = 0.0182, 5.49e6.',
               'These follow from the mechanism parameters alone; they are not analysis results. They agree with the',
               f'package\'s own evaluation of the mechanism ({RUN}/output/arrhenius_intrinsic.csv: gamma_clean 0.011894,',
               'k_diss_clean 4.885e6 1/s at 650 K).',
               'Formula: doc section 2.2 (sticking coefficients; Motz-Wise off by default).',
               f'R0: {RUN}/chem_annotated.yaml line 3494 (A 0.03, b 0, Ea 5000 J/mol); R1 = surface reaction #1 (A 0.008).',
               'Six load warnings "Sticking coefficient is greater than 1" (run.log line 9); N2 + 2 vacantX <=> [Pt]N=N[Pt](57),',
               'A = 1.66, YAML line 3533.'])


def _hump(x0, y0, x1, y1, h, n=40):
    pts = []
    for i in range(n + 1):
        t = i / n
        s = 3 * t * t - 2 * t ** 3
        pts.append((x0 + (x1 - x0) * t, y0 + (y1 - y0) * s - h * math.sin(math.pi * t) ** 1.4))
    return pts


def s16_thermo(d):
    sl = d.slide('Reverse rates come from thermodynamics, so equilibrium is always respected',
                 'Cantera never reads a reverse rate from the file: it divides the forward rate constant by K.',
                 part=PC, tag='concept')
    # free energy at 650 K: 2 H* lies ABOVE H2 + 2*, because H2 loses its gas entropy
    level(sl, 0.85, 4.35, 1.5, 'H₂(g) + 2 ∗', size=14)
    level(sl, 4.35, 2.75, 1.5, '2 H∗', size=14, label_side='above')
    poly(sl, _hump(2.35, 4.35, 4.35, 2.75, 0.18), color='ink', lw=1.75)
    line(sl, 2.35, 4.35, 6.3, 4.35, color='faint', dash=True)
    line(sl, 5.85, 2.75, 6.3, 2.75, color='faint', dash=True)
    arrow(sl, 6.15, 2.8, 6.15, 4.3, color='ink', lw=1.5, start=True)
    text(sl, 3.95, 3.3, 2.1, 0.4, 'ΔG° = +40.5 kJ/mol', size=16, font=MATH, align='r')
    text(sl, 0.75, 1.72, 5.6, 0.62, 'At 650 K for R0 on Cu: ΔH° = −38 kJ/mol, but H₂ loses its gas entropy '
         '(ΔS° = −121 J mol^{−1} K^{−1}), so ΔG° > 0 and K = 5.5 × 10^{−4}.', size=14)
    text(sl, 0.75, 4.62, 4.0, 0.3, 'schematic free-energy diagram, 650 K, 1 bar', size=12, color='faint', italic=True)
    card(sl, 6.83, 1.66, 5.95, 3.3, 'What Cantera computes')
    eqs = ['H°/RT = a_{1} + a_{2}T/2 + … + a_{6}/T', 'S°/R = a_{1} ln T + a_{2}T + … + a_{7}', 'G° = H° − T S°',
           'K_{c} = exp(−ΔG°/RT) Π_{k} (C_{k}°)^{ν(k)}', 'k_{r} = k_{f} / K_{c}']
    for i, e in enumerate(eqs):
        text(sl, 7.0, 2.1 + i * 0.52, 5.7, 0.5, e, size=18, font=MATH)
    caveat(sl, 0.55, 5.1, W - 2 * LM, 0.78, 'Shift one species’ enthalpy by δH and every reverse rate it appears in '
           'changes by e^{δH/RT}: 50 kJ/mol at 650 K is a factor of 1 × 10^{4}. Forward rates stay as they were.')
    actions(sl, ['Each species’ H and S come from its NASA polynomials.',
                 'ΔG° of the reaction gives its equilibrium constant K.',
                 'The reverse rate constant is k_{f} divided by K.',
                 'So editing an energy changes reverse rates, not forward ones.'])
    notes(sl, ['docs/simulation_workflow.md section 2.3. NASA 7-coefficient polynomials, two temperature ranges.',
               'Standard concentrations: P°/RT for gas, Γ/σ_k for surface species.',
               f'R0 numbers from the NASA polynomials of H2(1), vacantX(8) and [PtH](33) in {RUN}/chem_annotated.yaml,',
               'gas entropy moved from the file\'s 1 atm to P° = 1 bar: at 650 K dH = -38.04 kJ/mol, dS = -120.9 J/(mol K),',
               'dG = +40.5 kJ/mol, K = exp(-dG/RT) = 5.5e-4 (at 1200 K: dG = +104 kJ/mol, K = 2.8e-5). The package gives',
               'the same (output/arrhenius_intrinsic.csv; check_K_thermo = 1: rates and thermo agree).',
               'Pitfall from the doc: standard_gibbs_RT for an ideal gas is at the current pressure, not P°.',
               'exp(50e3 / (8.314 x 650)) = 1.04e4. This is why Stage 10 never edits thermochemistry in place.',
               'The energy diagram is schematic; the small hump stands for R0\'s 5 kJ/mol sticking barrier.'])


def s17_balance(d):
    sl = d.slide('Coverages change at their net production rates; steady state is when all stop',
                 'Every reaction adds to or takes from each species. At steady state the sums are all zero.',
                 part=PC, tag='concept')
    cx, cy = 3.1, 3.55                                   # labels stay inside the left margin
    circle(sl, cx, cy, 0.95, fill='H', line=None, content='H∗', size=20, color='white', bold=True)
    spokes = [(-2.05, -1.25, '#0  H₂ + 2 ∗ ⇌ 2 H∗'), (1.0, -1.35, '#19  CH₄ + 2 ∗ ⇌ H∗ + CH₃∗'),
              (-2.2, 1.25, '#26  H + ∗ ⇌ H∗'), (1.05, 1.2, '#37  NH₃ + 2 ∗ ⇌ H∗ + NH₂∗')]
    for dx, dy, t in spokes:
        ex, ey = cx + dx * 0.55, cy + dy * 0.6
        arrow(sl, ex, ey, cx + dx * 0.2, cy + dy * 0.28, color='H', lw=1.5, start=True)
        text(sl, ex - 1.4 if dx < 0 else ex - 0.2, ey - 0.42 if dy < 0 else ey + 0.02, 3.0, 0.4, t, size=13,
             color='ink', align='l')
    text(sl, 0.65, 5.2, 4.6, 0.4, 'H∗ appears in 126 surface reactions', size=14, color='muted')
    text(sl, 5.15, 2.05, 3.6, 0.5, 'ṡ_{k} = Σ_{i} ν_{k,i} q_{i}', size=22, font=MATH, align='c')
    text(sl, 5.15, 2.75, 3.6, 0.5, 'dθ_{k}/dt = σ_{k} ṡ_{k} / Γ', size=22, font=MATH, align='c')
    text(sl, 5.15, 3.5, 3.6, 1.2, 'for every surface species k. Σ θ stays 1, because every reaction conserves sites.',
         size=14, color='muted', align='c')
    X, Y, Wd, Ht = 9.6, 1.95, 3.0, 2.6
    axes(sl, X, Y, Wd, Ht, 'time (log scale)', 'coverage', ylabel_w=1.0)
    curve(sl, lambda t: 0.08 + 0.9 / (1 + math.exp(3.2 * (t - 1.2))), 0, 3, X, Y, Wd, Ht, color='faint', lw=2)
    curve(sl, lambda t: 0.84 / (1 + math.exp(-3.2 * (t - 1.2))), 0, 3, X, Y, Wd, Ht, color='H', lw=2.5)
    text(sl, X + 2.0, Y + 0.05, 1.0, 0.3, 'θ_{H∗}', size=14, color='H', font=MATH)
    text(sl, X + 2.55, Y + 2.0, 0.6, 0.3, 'θ_{∗}', size=14, color='muted', font=MATH)
    line(sl, X + 2.05, Y, X + 2.05, Y + Ht, color='ok', lw=1.5, dash=True)
    text(sl, X + 0.6, Y + Ht + 0.32, 2.6, 0.5, 'steady state: every ṡ_{k} = 0', size=14, color='ok', bold=True)
    schematic(sl, X, Y - 0.32)
    actions(sl, ['Add up what every reaction makes and uses of one species.',
                 'That sum is its net production rate ṡ.',
                 'Its coverage changes at σ\u00a0ṡ\u00a0/\u00a0Γ.',          # the formula never splits
                 'Steady state: every ṡ is zero, and the coverages still add to 1.'])
    notes(sl, ['docs/simulation_workflow.md section 2.4 (conservation equations and the steady-state problem).',
               f'Reactions shown are real: surface reactions #0, #19, #26 and #37 of {RUN}/chem_annotated.yaml, four of',
               'the 126 that contain [PtH](33) (H = H(16), CH3* = [CH3][Pt](52), NH2* = N[Pt](61)).',
               '#26 is the H-atom route: H atoms adsorb from, and desorb to, the gas.',
               'The plot is schematic: shape only, not a result.'])


def s18_units(d):
    sl = d.slide('Every rate in two units: per site (TOF) and per area (the literature’s)',
                 'The site density Γ converts one into the other. Nothing is solved again.',
                 part=PC, tag='package · per_area.py')
    # left: the bridge
    step(sl, 0.75, 1.85, 2.45, 1.0, ['**per site**', 'TOF, s^{−1}'], size=15)
    step(sl, 4.05, 1.85, 2.45, 1.0, ['**per area**', 'J, mol H₂ m^{−2} s^{−1}'], size=15)
    arrow(sl, 3.2, 2.35, 4.05, 2.35, color='ink', lw=2)
    text(sl, 3.08, 1.92, 1.1, 0.36, '× Γ', size=16, font=MATH, align='c', bold=True)
    text(sl, 0.75, 3.05, 5.75, 0.62, 'Γ = 2.943 × 10^{−5} mol sites per m² of Cu(111), from the YAML', size=14)
    text(sl, 0.75, 3.62, 5.75, 0.9, 'Clean surface at 650 K: k_{d} = γ Z / p = 1.44 × 10^{−3} mol H₂ m^{−2} s^{−1} '
         'Pa^{−1}, a property of the mechanism.', size=14, color='muted')
    # right: what is added
    card(sl, 6.83, 1.66, 5.95, 3.26, 'Format 2: added to every results table')
    rows = [('J_{diss}, J_{rec}, J_{net}', 'TOF × Γ', 'mol H₂ m^{−2} s^{−1}'),
            ('J_{ex}', 'min(J_{diss}, J_{rec})', 'mol H₂ m^{−2} s^{−1}'),
            ('k_{d,eff}', 'J_{diss} / p', 'mol H₂ m^{−2} s^{−1} Pa^{−1}'),
            ('γ_{eff}', 'J_{diss} / Z', '–'),
            ('C_{H}', 'θ_{H} Γ', 'mol H m^{−2}')]
    for i, (q, f, u) in enumerate(rows):
        yy = 2.14 + i * 0.5
        if i % 2 == 0:
            shape(sl, 'rect', 6.93, yy - 0.04, 5.75, 0.46, fill='white', line=None)
        text(sl, 7.0, yy, 1.75, 0.4, q, size=14, font=MATH)
        text(sl, 8.75, yy, 1.85, 0.4, f, size=14, font=MATH)
        text(sl, 10.55, yy, 2.2, 0.4, u, size=13, color='muted')
    text(sl, 7.0, 4.62, 5.6, 0.3, 'Z = p / √(2π M R T): H₂ hitting each m² per second', size=13, color='muted')
    caveat(sl, 0.55, 5.1, W - 2 * LM, 0.78, 'Per atomic area of an ideal Cu(111) facet: literature values per geometric '
           'area can be higher by a roughness factor. Converting the units does not change any activation energy.')
    actions(sl, ['Count events per site: the TOF.',
                 'Multiply by Γ: the same events per m².',
                 'Divide by p, or by Z: the literature’s k_{d} and the sticking probability.',
                 'Keep both: per-site tables and figures, and their _per_area copies.'])
    notes(sl, ['hydrogen_surface_kinetics/per_area.py: definitions (lines 4-13), add_columns (lines 61-73: J = TOF x Γ,',
               'J_exchange = min, kd_eff = J_diss/p, gamma_eff = J_diss/Z, C_H = θ_H Γ), impingement_flux (56-58),',
               'clean_surface (88-104). Doc Stage 7 ("Per area", Format 2) and docs/rates_units_equilibrium.md.',
               f'Γ: {RUN}/chem_annotated.yaml site-density 2.943e-05 mol/m2 (line 42).',
               'kd = γ Z / p at 650 K, 1 bar: 0.0119 x 1.209e4 / 1e5 = 1.44e-3 mol m-2 s-1 Pa-1; the package gives',
               f'kd_clean_mol_m2_s_Pa = 0.001438 ({RUN}/output/arrhenius_intrinsic.csv). A mechanism property, not a result.',
               'Every figure with a rate panel has a _per_area copy; the per-site versions stay as they are.'])


SLIDES = [s11_two_phases, s12_coverage, s13_site_factor, s14_mass_action, s15_sticking, s16_thermo, s17_balance,
          s18_units]
