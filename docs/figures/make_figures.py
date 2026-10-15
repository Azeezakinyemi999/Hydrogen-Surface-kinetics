"""Draw the figures for docs/simulation_workflow.md.

Everything here comes from the fictitious toy mechanism (toy_mechanism.yaml) or is a
schematic, so the figures illustrate the workflow without showing results from any
real mechanism. Re-run after changing the toy mechanism or the plotting style:

    python make_figures.py

Needs Cantera, NumPy and matplotlib, plus Graphviz's `dot` for the flux-diagram figure.
"""
import os
import shutil
import subprocess

import cantera as ct
import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
YAML = os.path.join(HERE, 'toy_mechanism.yaml')

# Validated categorical palette (fixed order) and neutrals
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    '#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948')
INK, SECONDARY, MUTED, GRID, PALE = '#33322f', '#55544f', '#8a8983', '#e4e3dd', '#f4f3ef'
SPECIES_COLORS = {'S': PALE, 'A(s)': BLUE, 'I(s)': ORANGE, 'I(ss)': AQUA}
SPECIES_LABELS = {'S': 'empty site', 'A(s)': 'A*', 'I(s)': 'I* (1 site)', 'I(ss)': 'I** (2 sites)'}

T_REF = 500.0      # K
P_A2 = 1.0e5       # Pa, probe-gas partial pressure held fixed in every toy calculation
R_SURF = 3         # index of the surface reaction A(s) + I(s) <=> AI + 2 S

# Same convergence settings as the steady-state notebook
ABS_TOL, REL_GROSS_TOL = 1e-8, 1e-12
PLATEAU_TOF_RTOL, PLATEAU_THETA_ATOL = 1e-4, 1e-4

plt.rcParams.update({
    'font.size': 9, 'axes.titlesize': 10, 'axes.labelsize': 9, 'axes.edgecolor': MUTED,
    'axes.labelcolor': INK, 'xtick.color': SECONDARY, 'ytick.color': SECONDARY, 'text.color': INK,
    'axes.titlelocation': 'left', 'figure.facecolor': 'white', 'savefig.facecolor': 'white',
    'legend.frameon': False, 'legend.fontsize': 8,
})


def style(ax, grid_axis='y'):
    ax.spines[['top', 'right']].set_visible(False)
    if grid_axis:
        ax.grid(axis=grid_axis, color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    path = os.path.join(HERE, name)
    fig.savefig(path, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', path)


# ----------------------------------------------------------------------------------
# Toy-model helpers (same approach as the notebook's fixed-gas reactor method)
# ----------------------------------------------------------------------------------
class FixedGasSurface:
    """Toy surface held in a constant gas: tiny surface in a huge volume, gas reactions off."""

    def __init__(self, T, p_I, theta0, multipliers=None, p_A2=P_A2):
        self.gas = ct.Solution(YAML, 'gas')
        self.surf = ct.Interface(YAML, 'surface', [self.gas])
        for i, m in (multipliers or {}).items():
            self.surf.set_multiplier(m, i)
        P = p_A2 + p_I
        X = {'A2': p_A2 / P}
        if p_I > 0:
            X['I'] = p_I / P
        self.gas.TPX = T, P, X
        self.surf.TP = T, P
        self.surf.coverages = theta0
        self.gamma = self.surf.site_density
        self.size = np.array([s.size for s in self.surf.species()])
        self.idx = np.array([self.surf.kinetics_species_index(n) for n in self.surf.species_names])
        reactor = ct.IdealGasReactor(self.gas, energy='off')
        reactor.volume = 1.0
        ct.ReactorSurface(self.surf, reactor, A=1e-10)
        self.net = ct.ReactorNet([reactor])
        self.net.rtol, self.net.atol = 1e-10, 1e-22

    def advance(self, t):
        self.net.advance(t)

    @property
    def theta(self):
        return dict(zip(self.surf.species_names, self.surf.coverages))

    @property
    def tof_probe(self):  # forward rate of A2 + 2 S -> 2 A(s), per site
        return self.surf.forward_rates_of_progress[0] / self.gamma

    @property
    def rate_AI(self):    # net formation of gas AI, per site
        return self.surf.net_production_rates[self.surf.kinetics_species_index('AI')] / self.gamma

    @property
    def site_count_factor(self):
        return float(np.sum(self.surf.concentrations * self.size) / self.gamma)

    def residual_ratio(self):
        """max_k |r_k| / (ABS_TOL + REL_GROSS_TOL * G_k); converged when <= 1."""
        to_theta = self.size / self.gamma
        r = self.surf.net_production_rates[self.idx] * to_theta
        gross = (self.surf.creation_rates[self.idx] + self.surf.destruction_rates[self.idx]) * to_theta
        return float(np.max(np.abs(r) / (ABS_TOL + REL_GROSS_TOL * gross)))


def steady(T, p_I, theta0, multipliers=None, t_max=1e8, p_A2=P_A2):
    """Integrate until TOF and coverages stop changing between half-decade checkpoints."""
    s = FixedGasSurface(T, p_I, theta0, multipliers, p_A2=p_A2)
    prev = None
    for t in np.logspace(-8, np.log10(t_max), int(2 * (np.log10(t_max) + 8)) + 1):
        s.advance(t)
        cur = (s.tof_probe, s.surf.coverages.copy())
        if prev is not None:
            d_tof = abs(cur[0] - prev[0]) / max(abs(prev[0]), 1e-300)
            if d_tof < 1e-7 and np.abs(cur[1] - prev[1]).max() < 1e-8:
                break
        prev = cur
    return s


# ----------------------------------------------------------------------------------
# Figure: closed batch reactor vs fixed-gas surface (schematic)
# ----------------------------------------------------------------------------------
def fig_reactor_models():
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.5))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 7)
        ax.axis('off')

    # Left: closed batch reactor
    ax = axes[0]
    ax.set_title('a) Closed batch reactor  (Stage 2)', color=INK)
    ax.add_patch(FancyBboxPatch((1.0, 1.6), 8.0, 4.6, boxstyle='round,pad=0.05', fc='#eef4fc', ec=BLUE, lw=1.4))
    ax.add_patch(Rectangle((1.0, 1.15), 8.0, 0.45, fc=ORANGE, ec='none'))
    ax.text(5.0, 0.7, 'catalyst surface, area $A$', ha='center', color=INK)
    ax.text(5.0, 5.55, 'gas volume $V$  (fixed)', ha='center', color=INK)
    for x, y, lab in [(2.2, 4.2, 'A$_2$'), (4.0, 3.4, 'I'), (6.3, 4.4, 'AI'), (7.6, 3.2, 'A$_2$'), (3.0, 2.5, 'AI')]:
        ax.text(x, y, lab, ha='center', va='center', color=INK)
    ax.add_patch(FancyArrowPatch((5.2, 2.8), (5.2, 1.75), arrowstyle='-|>', mutation_scale=10, color=MUTED))
    ax.add_patch(FancyArrowPatch((5.8, 1.75), (5.8, 2.8), arrowstyle='-|>', mutation_scale=10, color=MUTED))
    ax.text(6.0, 2.25, 'adsorb / desorb', color=SECONDARY, va='center', fontsize=8)
    ax.text(5.0, 6.55, 'gas composition changes as the surface reacts', ha='center', color=SECONDARY, fontsize=8)

    # Right: fixed-gas surface
    ax = axes[1]
    ax.set_title('b) Fixed-gas surface  (Stages 4-5)', color=INK)
    ax.add_patch(FancyBboxPatch((0.6, 1.6), 8.8, 4.6, boxstyle='round,pad=0.05', fc='#eef4fc', ec=BLUE, lw=1.4,
                                linestyle=(0, (4, 3))))
    ax.add_patch(Rectangle((4.4, 1.15), 1.2, 0.45, fc=ORANGE, ec='none'))
    ax.text(5.0, 0.7, 'tiny surface:  $A/V = 10^{-10}$ m$^{-1}$', ha='center', color=INK)
    ax.text(5.0, 5.55, 'huge gas reservoir, gas reactions off', ha='center', color=INK)
    ax.text(5.0, 3.7, '$T$, $p_k$ held constant', ha='center', va='center', color=INK, fontsize=10)
    ax.add_patch(FancyArrowPatch((4.8, 2.7), (4.8, 1.75), arrowstyle='-|>', mutation_scale=10, color=MUTED))
    ax.add_patch(FancyArrowPatch((5.2, 1.75), (5.2, 2.7), arrowstyle='-|>', mutation_scale=10, color=MUTED))
    ax.text(5.0, 6.55, 'only the coverages evolve, towards steady state', ha='center', color=SECONDARY, fontsize=8)
    fig.tight_layout()
    save(fig, 'fig_reactor_models.png')


# ----------------------------------------------------------------------------------
# Figure: counting sites with one- and two-site species (schematic + numbers)
# ----------------------------------------------------------------------------------
def fig_site_counting():
    # 10 sites; occupants listed left to right with the number of sites they take
    layout = [('A(s)', 1), ('S', 1), ('I(s)', 1), ('I(ss)', 2), ('S', 1), ('A(s)', 1), ('A(s)', 1), ('I(ss)', 2)]
    n_sites = sum(n for _, n in layout)
    fig = plt.figure(figsize=(9.2, 3.4))
    ax = fig.add_axes([0.02, 0.42, 0.96, 0.5])
    ax.set_xlim(-0.2, n_sites + 0.2)
    ax.set_ylim(-0.9, 1.6)
    ax.axis('off')
    x = 0
    for sp, n in layout:
        ax.add_patch(FancyBboxPatch((x + 0.06, 0.06), n - 0.12, 0.88, boxstyle='round,pad=0.02,rounding_size=0.12',
                                    fc=SPECIES_COLORS[sp], ec=MUTED if sp == 'S' else 'white', lw=1.0))
        ax.text(x + n / 2, 0.5, SPECIES_LABELS[sp].split(' (')[0], ha='center', va='center',
                color=INK if sp == 'S' else 'white', fontsize=9, fontweight='bold' if sp != 'S' else None)
        x += n
    for i in range(n_sites + 1):
        ax.plot([i, i], [-0.25, -0.05], color=MUTED, lw=0.8)
    for i in range(n_sites):
        ax.text(i + 0.5, -0.5, str(i + 1), ha='center', va='center', color=MUTED, fontsize=7)
    ax.text(n_sites / 2, 1.35, f'{n_sites} surface sites, {len(layout)} surface "entities" (adsorbates + empty sites)',
            ha='center', color=INK)

    counts = {}
    for sp, n in layout:
        counts.setdefault(sp, [0, n])[0] += 1
    theta = {sp: c * n / n_sites for sp, (c, n) in counts.items()}
    entities = sum(c for c, _ in counts.values())
    lines = ['Coverage (fraction of sites):  ' + ',  '.join(
        f'$\\theta$({SPECIES_LABELS[sp].split(" (")[0]}) = {theta[sp]:.1f}' for sp in ['A(s)', 'S', 'I(s)', 'I(ss)'])
        + f'    (sum = {sum(theta.values()):.1f})',
        f'Cantera keeps the number of entities at $\\Gamma$:  $\\sum_k C_k = \\Gamma$  '
        f'$\\Rightarrow$  $C_k = \\Gamma\\,\\cdot$(entities of $k$)/{entities}',
        f'Sites implied:  $\\sum_k \\sigma_k C_k = \\Gamma \\cdot {n_sites}/{entities}$  '
        f'$\\Rightarrow$  site-count factor $f = {n_sites / entities:.2f}$   '
        f'(empty-site concentration $C_* = f\\,\\theta_*\\,\\Gamma$)']
    for i, line in enumerate(lines):
        fig.text(0.03, 0.31 - i * 0.11, line, color=INK if i else SECONDARY, fontsize=9)
    save(fig, 'fig_site_counting.png')


# ----------------------------------------------------------------------------------
# Figure: approach to steady state, checkpoints and the acceptance test
# ----------------------------------------------------------------------------------
def fig_convergence():
    p_I = 10.0
    s = FixedGasSurface(T_REF, p_I, 'S:1')
    t_dense = np.logspace(-9, 1, 161)
    names = ['S', 'A(s)', 'I(s)', 'I(ss)']
    cov, tof, ratio = [], [], []
    for t in t_dense:
        s.advance(t)
        th = s.theta
        cov.append([th[n] for n in names])
        tof.append(s.tof_probe)
        ratio.append(s.residual_ratio())
    cov = np.array(cov)

    # Notebook acceptance logic on half-decade checkpoints
    s = FixedGasSurface(T_REF, p_I, 'S:1')
    t_chk = np.logspace(-8, 1, 19)
    prev, accept, chk = None, None, []
    for t in t_chk:
        s.advance(t)
        cur = (s.tof_probe, s.surf.coverages.copy())
        rr = s.residual_ratio()
        ok_plateau = False
        if prev is not None:
            ok_plateau = (abs(cur[0] - prev[0]) / max(abs(prev[0]), 1e-300) < PLATEAU_TOF_RTOL
                          and np.abs(cur[1] - prev[1]).max() < PLATEAU_THETA_ATOL)
        chk.append((t, rr, cur[0], ok_plateau))
        if accept is None and prev is not None and rr <= 1 and ok_plateau:
            accept = (t, rr, cur[0])
        prev = cur

    fig, axes = plt.subplots(3, 1, figsize=(7.6, 7.4), sharex=True)
    ax = axes[0]
    for j, n in enumerate(names):
        ax.plot(t_dense, cov[:, j], color=SPECIES_COLORS[n] if n != 'S' else MUTED, lw=2,
                label=SPECIES_LABELS[n])
    ax.set_ylabel('coverage $\\theta_k$')
    ax.set_title('a) Coverages, starting from a clean surface')
    ax.legend(loc='upper right', ncol=2, fontsize=8)
    ax.set_ylim(0, 1.05)
    style(ax)

    ax = axes[1]
    ax.plot(t_dense, tof, color=BLUE, lw=2)
    ax.set_yscale('log')
    ax.set_ylabel('forward TOF of probe\nreaction (s$^{-1}$)')
    ax.set_title('b) Probe-reaction rate')
    style(ax)

    ax = axes[2]
    ax.plot(t_dense, ratio, color=INK, lw=1.5, label='$\\max_k |r_k| \\,/\\, (\\varepsilon_{abs} + \\varepsilon_{rel} G_k)$')
    ax.axhline(1.0, color=RED, lw=1, ls='--')
    ax.text(t_dense[0], 1.6, 'residual test passes below this line', color=RED, fontsize=8)
    for t, rr, _, ok_plateau in chk:
        ax.plot(t, rr, 'o', ms=6, mfc=AQUA if ok_plateau else 'white', mec=AQUA if ok_plateau else MUTED, mew=1.2)
    if accept:
        ax.annotate(f'accepted at t = {accept[0]:.0e} s\n(residual test and plateau both pass)',
                    xy=(accept[0], accept[1]), xytext=(2e-7, 3e-2),
                    arrowprops=dict(arrowstyle='-|>', color=INK, lw=0.8), fontsize=8, va='center')
    ax.plot([], [], 'o', mfc='white', mec=MUTED, label='checkpoint (half decade)')
    ax.plot([], [], 'o', mfc=AQUA, mec=AQUA, label='checkpoint where TOF and $\\theta$ have plateaued')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.set_xlabel('surface time  t (s)')
    ax.set_ylabel('normalised residual')
    ax.set_title('c) Convergence test')
    ax.legend(loc='center left', fontsize=8)
    style(ax)
    fig.suptitle(f'Toy model, {T_REF:.0f} K, $p_{{A_2}}$ = 1 bar, $p_I$ = {p_I:g} Pa', x=0.01, ha='left',
                 fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_convergence.png')


# ----------------------------------------------------------------------------------
# Figure: poisoning curve and the site-blocking decomposition (surface reaction off)
# ----------------------------------------------------------------------------------
def fig_poisoning_curve():
    mult = {R_SURF: 0.0}           # pure site blocking: I cannot react away
    p_list = np.logspace(-6, 3, 28)
    ref = steady(T_REF, 0.0, 'S:1', mult)
    th_ref, tof_ref, f_ref = ref.theta['S'], ref.tof_probe, ref.site_count_factor
    rows = []
    for p in p_list:
        s = steady(T_REF, p, 'S:1', mult)
        th = s.theta
        rows.append((p, th['S'], th['A(s)'], th['I(s)'], th['I(ss)'], s.site_count_factor, s.tof_probe))
    rows = np.array(rows)
    p, free, a, i1, i2, f, tof = rows.T
    pf = tof / tof_ref
    pred_theta = (free / th_ref) ** 2
    pred_full = pred_theta * (f / f_ref) ** 2

    fig, axes = plt.subplots(1, 3, figsize=(11.5, 3.6))
    ax = axes[0]
    for y, n in [(free, 'S'), (a, 'A(s)'), (i1, 'I(s)'), (i2, 'I(ss)')]:
        ax.plot(p, y, color=SPECIES_COLORS[n] if n != 'S' else MUTED, lw=2, label=SPECIES_LABELS[n])
    ax.set_xscale('log')
    ax.set_ylim(0, 1.02)
    ax.set_xlabel('impurity pressure $p_I$ (Pa)')
    ax.set_ylabel('coverage')
    ax.set_title('a) Who takes the sites')
    ax.legend(fontsize=7.5, loc='center right')
    style(ax)

    ax = axes[1]
    ax.plot(p, f, color=AQUA, lw=2)
    ax.axhline(1.0, color=MUTED, lw=1, ls='--')
    ax.set_xscale('log')
    ax.set_xlabel('impurity pressure $p_I$ (Pa)')
    ax.set_ylabel('$f = \\sum_k \\sigma_k C_k / \\Gamma$')
    ax.set_title('b) Site-count factor')
    style(ax)

    ax = axes[2]
    ax.plot(p, pred_theta, color=MUTED, lw=1.5, ls='--', label='$(\\theta_*/\\theta_{*,\\mathrm{ref}})^2$  (coverages only)')
    ax.plot(p, pred_full, color=BLUE, lw=2, label='$\\times\\,(f/f_\\mathrm{ref})^2$  (full prediction)')
    ax.plot(p, pf, 'o', ms=4.5, mfc='white', mec=INK, mew=1, label='poisoning factor $\\Pi$ from rates')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('impurity pressure $p_I$ (Pa)')
    ax.set_ylabel('$\\Pi = \\mathrm{TOF}_f / \\mathrm{TOF}_{f,\\mathrm{ref}}$')
    ax.set_title('c) Poisoning factor')
    ax.legend(fontsize=7.5, loc='lower left')
    style(ax)
    fig.suptitle(f'Toy model, pure site blocking (surface reaction switched off), {T_REF:.0f} K, $p_{{A_2}}$ = 1 bar',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_poisoning_curve.png')


# ----------------------------------------------------------------------------------
# Figure: two starting surfaces, two steady states
# ----------------------------------------------------------------------------------
def fig_bistability():
    p_list = np.logspace(0.5, 2.0, 25)
    rows = []
    for p in p_list:
        a = steady(T_REF, p, 'S:1', t_max=1e7)
        b = steady(T_REF, p, 'I(s):1', t_max=1e7)
        ia = a.theta['I(s)'] + a.theta['I(ss)']
        ib = b.theta['I(s)'] + b.theta['I(ss)']
        rows.append((p, ia, ib, a.tof_probe, b.tof_probe))
    rows = np.array(rows)
    p, ia, ib, ta, tb = rows.T
    differ = (np.abs(ia - ib) > 1e-3) | (np.abs(ta - tb) / np.maximum(ta, tb) > 1e-2)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.7))
    for ax, ya, yb, ylabel, title, log in [
            (axes[0], ia, ib, 'impurity coverage  $\\theta_{I*} + \\theta_{I**}$', 'a) Surface state', False),
            (axes[1], ta, tb, 'forward TOF of probe reaction (s$^{-1}$)', 'b) Probe-reaction rate', True)]:
        if differ.any():
            lo, hi = p[differ].min(), p[differ].max()
            ax.axvspan(lo / 1.03, hi * 1.03, color=YELLOW, alpha=0.18, lw=0)
            ax.text(np.sqrt(lo * hi), 0.42 if not log else 0.55, 'starts\ndisagree', transform=ax.get_xaxis_transform(),
                    ha='center', va='center', fontsize=8, color=SECONDARY)
        ax.plot(p, ya, 'o-', color=BLUE, ms=4.5, lw=1.5, label='start A: clean surface')
        ax.plot(p, yb, 's--', color=ORANGE, mfc='white', ms=4.5, lw=1.2, label='start B: impurity-covered surface')
        ax.set_xscale('log')
        if log:
            ax.set_yscale('log')
        ax.set_xlabel('impurity pressure $p_I$ (Pa)')
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        style(ax)
    axes[0].set_ylim(0, 1.1)
    axes[1].legend(loc='lower left', fontsize=8)
    fig.suptitle(f'Toy model with the surface reaction A* + I* $\\to$ AI, {T_REF:.0f} K, $p_{{A_2}}$ = 1 bar',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_bistability.png')


# ----------------------------------------------------------------------------------
# Figure: how to read a flux diagram (Graphviz, same width rule as RMG)
# ----------------------------------------------------------------------------------
def fig_flux_anatomy():
    if shutil.which('dot') is None:
        print('skipping fig_flux_anatomy.png: Graphviz `dot` not found')
        return
    s = steady(T_REF, 10.0, 'S:1')
    q = s.surf.net_rates_of_progress / s.gamma          # per site
    # Flux pairs (reactant -> product sharing atoms), as RMG stores them; empty sites omitted
    pairs = {0: [('A2', 'A(s)')], 1: [('I', 'I(s)')], 2: [('I(s)', 'I(ss)')], 3: [('A(s)', 'AI'), ('I(s)', 'AI')]}
    flux = {}
    for i, plist in pairs.items():
        for a, b in plist:
            flux[(a, b)] = flux.get((a, b), 0.0) + q[i]
    fmax = max(abs(v) for v in flux.values())
    tol, wmax = 1e-6, 10.0
    labels = {'A2': 'A2\n(gas)', 'I': 'I\n(gas)', 'AI': 'AI\n(gas)', 'A(s)': 'A*', 'I(s)': 'I*', 'I(ss)': 'I**'}
    lines = ['digraph flux {', '  rankdir=LR; bgcolor="white"; pad=0.2; nodesep=0.7; ranksep=1.6;',
             '  node [shape=box, style="rounded,filled", fillcolor="#f4f3ef", color="#8a8983", '
             'fontname="Helvetica", fontsize=12, penwidth=1.2];',
             '  edge [fontname="Helvetica", fontsize=10, fontcolor="#55544f", color="#2a78d6"];']
    for n, lab in labels.items():
        lines.append(f'  "{n}" [label="{lab}"];')
    for (a, b), v in flux.items():
        x = abs(v) / fmax
        if x < tol:
            lines.append(f'  "{a}" -> "{b}" [dir=none, penwidth=1.0, style=dotted, color="#8a8983", '
                         f'label=<<font point-size="9">net flux below tolerance:<br/>not drawn in RMG frames</font>>];')
            continue
        w = wmax * (1 - np.log10(x) / np.log10(tol))
        direction = 'forward' if v > 0 else 'back'
        lines.append(f'  "{a}" -> "{b}" [dir={direction}, penwidth={w:.2f}, arrowsize={0.6 + 0.04 * w:.2f}, '
                     f'label=<<table border="0" cellborder="0" cellpadding="2" bgcolor="white"><tr><td>'
                     f'x = {x:.1g}, width {w:.1f}</td></tr></table>>];')
    lines.append('  labelloc="t"; fontname="Helvetica"; fontsize=12; '
                 'label="Toy flux diagram at steady state (x = |flux| / max flux)";')
    lines.append('}')
    dot_path = os.path.join(HERE, 'fig_flux_anatomy.dot')
    with open(dot_path, 'w') as f:
        f.write('\n'.join(lines))
    subprocess.run(['dot', '-Tpng', '-Gdpi=170', dot_path, '-o', os.path.join(HERE, 'fig_flux_anatomy.png')], check=True)
    os.remove(dot_path)
    print('wrote', os.path.join(HERE, 'fig_flux_anatomy.png'))


# ----------------------------------------------------------------------------------
# Figure: extended analyses on the toy (temperature sweep and degree of rate control)
# ----------------------------------------------------------------------------------
def fig_extended():
    temps = np.arange(450.0, 701.0, 25.0)
    series = {}
    for label, p_I in [('A$_2$ only', 0.0), ('with $p_I$ = 1 Pa', 1.0)]:
        series[label] = np.array([steady(T, p_I, 'S:1', {R_SURF: 0.0}).tof_probe for T in temps])

    # Degree of rate control for the net formation of AI (an overall rate), clean start
    T_drc, p_drc, delta = T_REF, 10.0, 0.01
    base = steady(T_drc, p_drc, 'S:1').rate_AI
    xrc = []
    for i in range(4):
        up = steady(T_drc, p_drc, 'S:1', {i: 1 + delta}).rate_AI
        dn = steady(T_drc, p_drc, 'S:1', {i: 1 - delta}).rate_AI
        xrc.append((np.log(up) - np.log(dn)) / (np.log(1 + delta) - np.log(1 - delta)))
    xrc = np.array(xrc)

    fig, axes = plt.subplots(1, 2, figsize=(10.0, 3.7))
    ax = axes[0]
    inv = 1000.0 / temps
    j = int(np.argmin(np.abs(temps - 550.0)))       # local E_app at 550 K by central difference
    for (label, tof), color, marker, dy in zip(series.items(), [BLUE, ORANGE], ['o', 's'], [3.0, 1 / 3.0]):
        ax.plot(inv, tof, marker=marker, color=color, lw=1.5, ms=4.5, label=label)
        dlnr = np.log(tof[j + 1]) - np.log(tof[j - 1])
        e_app = -ct.gas_constant * dlnr / (1.0 / temps[j + 1] - 1.0 / temps[j - 1]) / 1e6
        ax.annotate(f'$E_{{app}}$(550 K) = {e_app:.0f} kJ/mol', xy=(inv[j], tof[j]), xytext=(inv[j] + 0.08, tof[j] * dy),
                    color=color, fontsize=8, va='center', arrowprops=dict(arrowstyle='-', color=color, lw=0.6))
    ax.set_yscale('log')
    ax.set_xlabel('1000 / T  (K$^{-1}$)')
    ax.set_ylabel('forward TOF of probe reaction (s$^{-1}$)')
    ax.set_title('a) Temperature sweep (Arrhenius plot)')
    ax.legend(loc='lower left', fontsize=8)
    style(ax)

    ax = axes[1]
    names = ['R1  A$_2$ + 2* $\\rightleftharpoons$ 2A*', 'R2  I + * $\\rightleftharpoons$ I*',
             'R3  I* + * $\\rightleftharpoons$ I**', 'R4  A* + I* $\\rightleftharpoons$ AI + 2*']
    y = np.arange(len(names))[::-1]
    ax.barh(y, xrc, color=[BLUE if v >= 0 else ORANGE for v in xrc], height=0.55)
    ax.axvline(0, color=MUTED, lw=0.8)
    for yi, v in zip(y, xrc):
        ax.text(v + (0.03 if v >= 0 else -0.03), yi, f'{v:+.2f}', va='center', ha='left' if v >= 0 else 'right', fontsize=8)
    ax.set_yticks(y, names, fontsize=8)
    ax.set_xlabel('degree of rate control $X_{RC,i}$')
    ax.set_title(f'b) Degree of rate control of AI formation  (sum = {xrc.sum():.2f})')
    lim = max(1.0, np.abs(xrc).max() * 1.35)
    ax.set_xlim(-lim, lim)
    style(ax, grid_axis='x')
    fig.suptitle(f'Toy model: (a) surface reaction switched off; (b) {T_drc:.0f} K, $p_I$ = {p_drc:g} Pa, $\\delta$ = {delta}',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_extended_analyses.png')


# ----------------------------------------------------------------------------------
# Format 2 (per area) and equilibrium coverage: docs/rates_units_equilibrium.md
# ----------------------------------------------------------------------------------
P0 = 1e5           # Pa, standard pressure of K
R_MOL = ct.gas_constant / 1e3


def per_area(s, p_A2=P_A2):
    """Format-2 numbers of a solved toy surface: fluxes per m2, kd_eff, sticking probability,
    exchange flux and hydrogen on the surface."""
    G = s.gamma * 1e3                                              # mol sites / m2
    M = s.gas.molecular_weights[s.gas.species_index('A2')] / 1e3   # kg / mol
    fwd = s.surf.forward_rates_of_progress[0] * 1e3                # mol A2 / m2 / s
    rev = s.surf.reverse_rates_of_progress[0] * 1e3
    T = s.surf.T
    Z = p_A2 / np.sqrt(2 * np.pi * M * R_MOL * T)
    th = s.theta
    return {'tof': fwd / G, 'J_diss': fwd, 'kd_eff': fwd / p_A2, 'gamma': fwd / Z, 'J_ex': min(fwd, rev),
            'C_H': th['A(s)'] * G, 'theta_H': th['A(s)'], 'theta_v': th['S'], 'G': G}


def K_eq(T):
    """Dimensionless K = theta_A^2 / (theta_v^2 p/P0) of A2 + 2 S <=> 2 A(s), from the thermochemistry."""
    gas = ct.Solution(YAML, 'gas')
    surf = ct.Interface(YAML, 'surface', [gas])
    g = lambda sp: sp.thermo.h(T) - T * sp.thermo.s(T)                          # noqa: E731
    dG = 2 * g(surf.species('A(s)')) - g(gas.species('A2')) - 2 * g(surf.species('S'))
    return np.exp(-dG / (ct.gas_constant * T)) * P0 / gas.species('A2').thermo.reference_pressure


def theta_eq(K, p):
    x = np.sqrt(K * np.asarray(p, dtype=float) / P0)
    return x / (1 + x)


def fig_per_area():
    """Format 1 against Format 2: the same temperature sweep per site and per area, and the
    effective adsorption coefficient against pressure."""
    temps = np.arange(450.0, 651.0, 25.0)
    cases = [('A$_2$ only', 0.0, INK, 'o'), ('with $p_I$ = 1 Pa', 1.0, BLUE, 's')]
    sweep = {label: [per_area(steady(T, p_I, 'S:1')) for T in temps] for label, p_I, _, _ in cases}
    pressures = np.logspace(np.log10(3e3), np.log10(3e6), 7)
    psweep = [per_area(steady(T_REF, 0.0, 'S:1', p_A2=p), p) for p in pressures]
    kd_clean = psweep[0]['J_diss'] / pressures[0] / psweep[0]['theta_v'] ** 2      # kd = J / (p theta_v^2)
    G = psweep[0]['G']
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    inv = 1000.0 / temps
    for (label, _, color, marker) in cases:
        for ax, key in zip(axes[:2], ('tof', 'kd_eff')):
            y = np.array([d[key] for d in sweep[label]])
            slope = np.polyfit(1 / temps, np.log(y), 1)[0]
            ax.plot(inv, y, marker=marker, color=color, lw=1.5, ms=4.5, label=f'{label}: $E_a$ = {-slope * R_MOL / 1e3:.1f} kJ/mol')
    axes[0].set_ylabel('TOF (s$^{-1}$, per site)')
    axes[0].set_title('a) Format 1: per site')
    axes[1].set_ylabel('$k_{d,\\mathrm{eff}}$ (mol A$_2$ m$^{-2}$ s$^{-1}$ Pa$^{-1}$)')
    axes[1].set_title(f'b) Format 2: per area ($\\Gamma$ = {G:.2g} mol m$^{{-2}}$)')
    for ax in axes[:2]:
        ax.set_yscale('log')
        ax.set_xlabel('1000 / T  (K$^{-1}$)')
        ax.legend(loc='lower left', fontsize=7.5)
        style(ax)
    ax = axes[2]
    ax.plot(pressures / 1e5, [d['kd_eff'] for d in psweep], 'o-', color=INK, lw=1.5, ms=4.5, label='$k_{d,\\mathrm{eff}} = k_d\\,\\theta_v^2$, A$_2$ only')
    ax.axhline(kd_clean, color=MUTED, ls='--', lw=1.2, label='$k_d$, clean surface (first order)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('p(A$_2$) (bar)')
    ax.set_ylabel('mol A$_2$ m$^{-2}$ s$^{-1}$ Pa$^{-1}$')
    ax.set_title(f'c) Against pressure, {T_REF:.0f} K')
    ax.legend(loc='lower left', fontsize=7.5)
    style(ax)
    fig.suptitle('Toy model: the same rates per site and per area; activation energies are unchanged, only the prefactor and units differ',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_per_area.png')


def fig_sticking_exchange():
    """Sticking probability, exchange flux and hydrogen on the surface against temperature."""
    temps = np.arange(450.0, 651.0, 25.0)
    cases = [('A$_2$ only', 0.0, INK, 'o'), ('with $p_I$ = 1 Pa', 1.0, BLUE, 's')]
    data = {label: [per_area(steady(T, p_I, 'S:1')) for T in temps] for label, p_I, _, _ in cases}
    clean = [per_area(FixedGasSurface(T, 0.0, 'S:1')) for T in temps]          # empty surface: gamma of the mechanism
    G = clean[0]['G']
    th = theta_eq(np.array([K_eq(T) for T in temps]), P_A2)
    J0 = np.array([c['J_diss'] for c in clean]) * (1 - th) ** 2                # kd p theta_v,eq^2
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    for label, _, color, marker in cases:
        d = data[label]
        axes[0].plot(temps, [x['gamma'] for x in d], marker=marker, color=color, lw=1.5, ms=4.5, label=label)
        axes[1].plot(1000 / temps, [x['J_ex'] for x in d], marker=marker, color=color, lw=1.5, ms=4.5, label=label)
        axes[2].plot(temps, [x['C_H'] for x in d], marker=marker, color=color, lw=1.5, ms=4.5, label=label)
    dash = dict(color=MUTED, ls='--', lw=1.2)
    axes[0].plot(temps, [c['gamma'] for c in clean], label='clean surface (mechanism)', **dash)
    axes[1].plot(1000 / temps, J0, label='clean surface at equilibrium, $J_0$', **dash)
    axes[2].plot(temps, th * G, label='clean surface at equilibrium', **dash)
    axes[0].set(yscale='log', xlabel='T (K)', ylabel='$\\gamma_\\mathrm{eff} = J_\\mathrm{diss}/Z$', title='a) Sticking probability')
    axes[1].set(yscale='log', xlabel='1000 / T  (K$^{-1}$)', ylabel='$J_\\mathrm{ex}$ (mol A$_2$ m$^{-2}$ s$^{-1}$)', title='b) Exchange flux')
    axes[2].set(xlabel='T (K)', ylabel='$C_\\mathrm{H}$ (mol H m$^{-2}$)', title='c) Hydrogen on the surface')
    for ax in axes:
        ax.legend(fontsize=7.5)
        style(ax)
    fig.suptitle(f'Toy model, p(A$_2$) = 1 bar: Format 2 quantities; $\\Gamma$ = {G:.2g} mol m$^{{-2}}$',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_sticking_exchange.png')


def fig_equilibrium():
    """Equilibrium coverage of a surface holding only A, and the equilibrium check of steady states."""
    p_line = np.logspace(np.log10(3e3), np.log10(3e6), 80)
    T_line = np.linspace(380.0, 720.0, 80)
    K_line = np.array([K_eq(T) for T in T_line])
    fig, axes = plt.subplots(1, 3, figsize=(13.0, 3.8))
    for T, color in ((450.0, ORANGE), (500.0, INK), (550.0, BLUE)):
        axes[0].plot(p_line / 1e5, theta_eq(K_eq(T), p_line), color=color, lw=1.5, label=f'{T:.0f} K')
    for p in (3e4, 3e5):
        d = per_area(steady(T_REF, 0.0, 'S:1', p_A2=p), p)
        axes[0].plot(p / 1e5, d['theta_H'], 'o', mfc='white', mec=INK, mew=1.4, ms=6)
    axes[0].plot([], [], 'o', mfc='white', mec=INK, mew=1.4, label='steady states, 500 K')
    for p, color in ((3e4, YELLOW), (1e5, MAGENTA), (3e5, GREEN)):
        axes[1].plot(T_line, theta_eq(K_line, p), color=color, lw=1.5, label=f'{p / 1e5:g} bar')
        i = np.nonzero(np.diff(np.sign(np.log(K_line * p / P0))))[0][0]           # K p/P0 = 1
        a, b = T_line[i], T_line[i + 1]
        fa, fb = np.log(K_line[i] * p / P0), np.log(K_line[i + 1] * p / P0)
        axes[1].plot(a - fa * (b - a) / (fb - fa), 0.5, 'D', mfc='white', mec=color, mew=1.5, ms=7)
    axes[0].set(xscale='log', xlabel='p(A$_2$) (bar)', ylabel='$\\theta_\\mathrm{A,eq}$', ylim=(0, 1.02),
                title='a) Isotherms, surface holding only A')
    axes[1].set(xlabel='T (K)', ylabel='$\\theta_\\mathrm{A,eq}$', ylim=(0, 1.02), title='b) Against T; diamonds: $T_{1/2}$')
    pts = (3e4, 1e5, 3e5)
    for j, (label, p_I, color) in enumerate((('A$_2$ only', 0.0, INK), ('with I, r = 10$^{-5}$', None, BLUE))):
        q = []
        for p in pts:
            d = per_area(steady(T_REF, 1e-5 * p if p_I is None else p_I, 'S:1', p_A2=p), p)
            q.append(d['theta_H'] / d['theta_v'] / np.sqrt(K_eq(T_REF) * p / P0))
        axes[2].plot(np.arange(3) + (j - 0.5) * 0.12, q, 'o', color=color, ms=6, label=label)
    axes[2].axhline(1.0, color=MUTED, ls='--', lw=1)
    axes[2].set_xticks(range(3), [f'{T_REF:.0f} K\n{p / 1e5:g} bar' for p in pts])
    axes[2].set(ylim=(0.998, 1.002), ylabel='$q = (\\theta_\\mathrm{A}/\\theta_v)\\,/\\,(K p/P^\\circ)^{1/2}$',
                title='c) Equilibrium check of steady states')
    for ax in axes:
        ax.legend(fontsize=7.5)
        style(ax)
    fig.suptitle('Toy model: equilibrium of A$_2$ + 2* $\\rightleftharpoons$ 2A* from the thermochemistry',
                 x=0.01, ha='left', fontsize=10, color=SECONDARY)
    fig.tight_layout()
    save(fig, 'fig_equilibrium.png')


if __name__ == '__main__':
    fig_reactor_models()
    fig_site_counting()
    fig_convergence()
    fig_poisoning_curve()
    fig_bistability()
    fig_flux_anatomy()
    fig_extended()
    fig_per_area()
    fig_sticking_exchange()
    fig_equilibrium()
