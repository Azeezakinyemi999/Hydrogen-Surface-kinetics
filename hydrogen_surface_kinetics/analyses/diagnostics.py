"""Is the mechanism believable? (documentation: Stage 10)

1. adsorption_checks: reaction enthalpies from the mechanism's own thermo vs literature.
2. switch-off what-ifs: for each entry of `switch_off` (a species: every surface reaction
   that makes or uses it; a reaction number; or a list of those together), the reactions
   get rate multiplier 0 and the main cases are solved again. [] (default) = none.
3. top_reactions: largest sources and sinks of chosen species at chosen conditions.

Outputs: mechanism_binding_check.csv + binding_check.png, mechanism_switch_off_<label>.csv
+ switch_off_<label>.png per what-if and mechanism_switch_off.csv with all of them,
top_reactions.csv.
"""
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from .. import solver as S
from ..labels import probe_label, species_label
from . import common


def reaction_thermo(model, reactants, products, T):
    """dH and dG (kJ/mol) of reactants -> products at T (reference-state thermo)."""
    def obj(name):
        return model.gas.species(name) if name in model.gas.species_names else model.surf.species(name)

    def total(side, fn):
        return sum(n * fn(obj(sp)) for sp, n in side.items())
    dH = total(products, lambda sp: sp.thermo.h(T)) - total(reactants, lambda sp: sp.thermo.h(T))
    dG = (total(products, lambda sp: sp.thermo.h(T) - T * sp.thermo.s(T))
          - total(reactants, lambda sp: sp.thermo.h(T) - T * sp.thermo.s(T)))
    return dH / 1e6, dG / 1e6


def binding_check(ctx):
    rows = []
    temps = sorted(ctx.cfg.temperatures)
    for label, reac, prod, (lo, hi), note in ctx.cfg.diagnostics['adsorption_checks']:
        dH298, _ = reaction_thermo(ctx.model, reac, prod, 298.15)
        verdict = 'OK' if lo <= dH298 <= hi else ('TOO STRONG' if dH298 < lo else 'TOO WEAK')
        out = {'adsorption': label, 'dH_298_mech': dH298, 'dH_literature': f'{lo} to {hi}',
               'lit_low': float(lo), 'lit_high': float(hi), 'verdict': verdict}
        for T in temps:
            out[f'dG_{T:.0f}'] = reaction_thermo(ctx.model, reac, prod, T)[1]
        out['literature_note'] = note
        rows.append(out)
    return pd.DataFrame(rows)


VERDICT_STYLE = {                       # color + marker shape, so the verdict is never color alone
    'OK':         dict(color=P.INK, marker='o', label='within range'),
    'TOO STRONG': dict(color=P.DIVERGING_HIGH, marker='<', label='too strong (more negative)'),
    'TOO WEAK':   dict(color=P.DIVERGING_LOW, marker='>', label='too weak'),
}


def plot_binding_check(df):
    """Literature range (gray bar) and the mechanism's value (marker) per adsorption."""
    n = len(df)
    y = np.arange(n)[::-1]                                   # first check at the top
    fig, ax = plt.subplots(figsize=(8.0, 0.30 * n + 1.1))
    lo_all = min(df['lit_low'].min(), df['dH_298_mech'].min())
    hi_all = max(df['lit_high'].max(), df['dH_298_mech'].max())
    pad = 0.06 * (hi_all - lo_all)
    for yi, r in zip(y, df.itertuples()):
        ax.plot([r.lit_low, r.lit_high], [yi, yi], color=P.RANGE_COLOR, lw=7, solid_capstyle='round', zorder=1)
        st = VERDICT_STYLE[r.verdict]
        ax.plot(r.dH_298_mech, yi, ls='', marker=st['marker'], ms=9, color=st['color'], mec='white', mew=0.8, zorder=3)
        if r.verdict != 'OK':                                # selective direct labels, in ink
            left = r.dH_298_mech < r.lit_low
            ax.annotate(r.verdict.lower(), (r.dH_298_mech, yi), xytext=(-8 if left else 8, 0),
                        textcoords='offset points', ha='right' if left else 'left', va='center',
                        fontsize=8, color=P.SECONDARY_INK)
    ax.set_yticks(y, [a.replace('->', '→') for a in df['adsorption']])
    ax.set_xlim(lo_all - 1.5 * pad, hi_all + 2.0 * pad)       # room for the 'too weak' labels
    ax.set_xlabel('adsorption enthalpy at 298 K (kJ/mol)   ← stronger binding', fontsize=10, color=P.INK)
    P.style_axis(ax)
    ax.grid(axis='y', visible=False)
    ax.grid(axis='x', color=P.GRID, linewidth=0.6)
    handles = [plt.Line2D([], [], color=P.RANGE_COLOR, lw=7, solid_capstyle='round', label='literature range')]
    handles += [plt.Line2D([], [], ls='', marker=st['marker'], ms=9, color=st['color'], label=st['label'])
                for v, st in VERDICT_STYLE.items() if v in set(df['verdict'])]
    fig.legend(handles=handles, loc='upper center', bbox_to_anchor=(0.5, 1.02), ncol=len(handles),
               frameon=False, fontsize=9)
    fig.suptitle('Adsorption enthalpies: mechanism vs literature', x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    return fig


def switch_off_items(cfg):
    """[(label, parts)] of the what-ifs in diagnostics['switch_off']. Each entry is a species
    name (every surface reaction that makes or uses it), a reaction number, or a list of
    those switched off together. None or [] = no what-if."""
    raw = cfg.diagnostics.get('switch_off')
    if raw is None:
        return []
    items = list(raw) if isinstance(raw, (list, tuple)) else [raw]
    out = []
    for item in items:
        parts = list(item) if isinstance(item, (list, tuple)) else [item]
        out.append(('+'.join(f'R{p}' if isinstance(p, (int, np.integer)) else str(p) for p in parts), parts))
    return out


def file_label(label):
    """'H(19)' -> 'H_19', 'R47' -> 'R47': the what-if's part of mechanism_switch_off_<label>.csv."""
    s = re.sub(r'_+', '_', re.sub(r'[^A-Za-z0-9.+-]+', '_', label))
    return re.sub(r'_*\+_*', '+', s).strip('_')


def switch_off_reactions(model, parts):
    """Surface reaction numbers of one what-if; a clear error for an unknown species or number."""
    rxns = set()
    for p in parts:
        if isinstance(p, (int, np.integer)) and not isinstance(p, bool):
            if not 0 <= p < model.surf.n_reactions:
                raise ValueError(f"diagnostics['switch_off']: reaction #{p} does not exist "
                                 f"(the surface has reactions 0 to {model.surf.n_reactions - 1})")
            rxns.add(int(p))
            continue
        found = [i for i, r in enumerate(model.surf.reactions()) if p in r.reactants or p in r.products]
        if not found:
            where = ('is not a species of the mechanism' if p not in model.surf.species_names + model.gas.species_names
                     else 'is in no surface reaction')
            raise ValueError(f"diagnostics['switch_off']: '{p}' {where}")
        rxns.update(found)
    return sorted(rxns)


def switch_off(ctx):
    """Every what-if of diagnostics['switch_off']: the listed reactions get rate multiplier 0
    and the main cases are solved again from a clean surface. Returns [(label, table)]."""
    from .equilibrium import K_thermo
    model, s = ctx.model, ctx.cfg.diagnostics
    items = switch_off_items(ctx.cfg)
    if not items:
        return []
    cases = s['switch_off_cases'] or C.case_labels(ctx.cfg)
    main = ctx.main()['summary']
    results = ctx.main()['results']
    wanted = {(pt['T'], round(pt['P'] / 1e5, 12)) for pt in common.points(ctx, 'diagnostics')}
    at_points = [(T, round(Pb, 12)) in wanted for T, Pb in zip(main['T_K'], main['P_ref_bar'].astype(float))]
    base = main[(main['basis'] == 'fixed_probe') & main['case'].isin(cases) & np.array(at_points, dtype=bool)]
    G, n = model.gamma, model.probe_sites
    K = {T: K_thermo(model, T) for T in sorted(set(base['T_K']))}

    def q(theta_H, theta_v, T, P_bar):
        return theta_H / theta_v / (K[T] * P_bar) ** (1.0 / n) if theta_v > 0 else np.nan     # p/P0 with P0 = 1 bar

    def job_func(job):
        theta, info = S.solve_steady_state(model, job['row']['T_K'], job['row']['P_total_bar'] * 1e5, job['row']['X'],
                                           S.clean_surface(model), multipliers=job['multipliers'])
        return {'theta': theta, 'info': info}

    tables = []
    for label, parts in items:
        rxns = switch_off_reactions(model, parts)
        multipliers = {i: 0.0 for i in rxns}
        listed = '; '.join(f'#{i}: {model.surf.reaction(i).equation}' for i in rxns)
        ctx.log(f'  switch-off what-if {label}: {len(rxns)} reaction(s) switched off')
        jobs = [common.job(f"switch_off {label} | {r.case} | {r.T_K:.0f} K | {r.P_ref_bar:g} bar",
                           results[r.key]['row'], multipliers=multipliers) for r in base.itertuples()]
        out = ctx.run_jobs(job_func, jobs, 'cache_diagnostics', describe=lambda res: S.describe(res['info']))
        rows = []
        for r, job in zip(base.itertuples(), jobs):
            res = out[job['key']]
            row = {'switched_off': label, 'reactions': listed, 'series': r.series, 'case': r.case, 'T_K': r.T_K,
                   'P_ref_bar': r.P_ref_bar, 'tof_normal': r.tof_fwd, 'tof_net_normal': r.tof_net,
                   'probe_consumed_normal': r.tof_probe_consumed, 'theta_free_normal': r.theta_free,
                   'theta_H_normal': r.theta_probe, 'q_normal': q(r.theta_probe, r.theta_free, r.T_K, r.P_probe_bar)}
            ok = res['status'] == 'done' and res['info']['converged']
            row['converged'] = ok
            if ok:                                       # rates with the switch-off still applied
                surf = common.load_state(model, job['row'], res['theta'], multipliers)
                fwd, _, net = model.probe_rates(surf)
                th = dict(zip(model.species_names, res['theta']))
                row.update({'tof_switched_off': fwd / G, 'tof_net_switched_off': net / G,
                            'probe_consumed_switched_off': -surf.net_production_rates[model.probe_kin_idx] / G,
                            'theta_free_switched_off': th[model.site], 'theta_H_switched_off': th[model.probe_adsorbate],
                            'q_switched_off': q(th[model.probe_adsorbate], th[model.site], r.T_K, r.P_probe_bar),
                            'top_switched_off': common.top_species(model, res['theta'])})
            rows.append(row)
        df = pd.DataFrame(rows)
        df['tof_ratio'] = df.get('tof_switched_off', np.nan) / df['tof_normal']
        ref = df[df['case'] == ctx.cfg.reference_label].set_index(['T_K', 'P_ref_bar'])
        for col in ('tof_normal', 'tof_switched_off'):
            if col in df and len(ref):
                ref_val = [ref[col].get((T, Pb), np.nan) for T, Pb in zip(df['T_K'], df['P_ref_bar'])]
                df[f'slowdown_{col[4:]}'] = np.array(ref_val, dtype=float) / df[col]
        tables.append((label, df))
    return tables


SPAN_LOG = 30.0           # a switch-off panel turns logarithmic when its values span more than this factor


def fit_y(ax, values, headroom=0.08):
    """The y range of a switch-off panel: the values with 1 in view, logarithmic when they span
    more than SPAN_LOG, else linear without an offset; `headroom` (fraction of the range) on top."""
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    lo, hi = (min(v.min(), 1.0), max(v.max(), 1.0)) if len(v) else (1.0, 1.0)
    if lo > 0 and hi / lo > SPAN_LOG:
        ax.set_yscale('log')
        span = np.log10(hi / lo)
        ax.set_ylim(lo / 10 ** (0.08 * span), hi * 10 ** (max(0.08, headroom) * span))
    else:
        pad = max(0.004, 0.08 * (hi - lo))
        ax.set_ylim(lo - pad, hi + max(pad, headroom * (hi - lo)))
        ax.ticklabel_format(axis='y', style='plain', useOffset=False)


def plot_switch_off(cfg, label, df, probe='H2'):
    """One what-if, per case (color) at each main point (x axis: the temperature series, then
    the pressure series): (a) dissociation TOF switched off / normal; (b) the slowdown of each
    impurity case, normal (dot) and switched off (ring); (c) the equilibrium check q of the
    probe step, normal and switched off. Needs only the settings and the what-if's table."""
    known = C.case_labels(cfg)
    cases = [c for c in known if c in set(df['case'])] + sorted(set(df['case']) - set(known))
    colors = P.colors_for(cases, first=cfg.reference_label)
    present = set(zip(df['T_K'].astype(float), df['P_ref_bar'].astype(float)))

    def col(d, name):
        return d[name].astype(float) if name in d else pd.Series(np.nan, index=d.index)

    fig, axes = plt.subplots(3, 1, figsize=(11, 10))
    pts = [P.point_axis(ax, cfg, present, label_series=(i == 0), labels=(i == 2)) for i, ax in enumerate(axes)][0]
    dot = dict(marker='o', ls='', ms=5, mec='white', mew=0.6)
    ring = dict(marker='o', ls='', ms=10, mfc='none', mew=1.3)
    for j, case in enumerate(cases):
        d = df[df['case'] == case]
        x = [pts.index((float(T), float(Pb))) + (j - (len(cases) - 1) / 2) * 0.08 for T, Pb in zip(d['T_K'], d['P_ref_bar'])]
        axes[0].plot(x, col(d, 'tof_ratio'), color=colors[case], **{**dot, 'ms': 6})
        if case != cfg.reference_label:                       # the reference's slowdown is 1 by definition
            axes[1].plot(x, col(d, 'slowdown_normal'), color=colors[case], **dot)
            axes[1].plot(x, col(d, 'slowdown_switched_off'), mec=colors[case], **ring)
        axes[2].plot(x, col(d, 'q_normal'), color=colors[case], **dot)
        axes[2].plot(x, col(d, 'q_switched_off'), mec=colors[case], **ring)
    others = df[df['case'] != cfg.reference_label]
    fit_y(axes[0], col(df, 'tof_ratio'), headroom=0.3)                   # room for the series names
    fit_y(axes[1], np.r_[col(others, 'slowdown_normal'), col(others, 'slowdown_switched_off')])
    fit_y(axes[2], np.r_[col(df, 'q_normal'), col(df, 'q_switched_off')])
    for ax in axes:
        ax.axhline(1.0, color=P.MUTED, lw=1, ls='--', zorder=1)
    P.style_axis(axes[0], 'a) Dissociation rate: switched off / normal', 'rate ratio')
    P.style_axis(axes[1], 'b) Slowdown by each impurity', 'slowdown $1/\\Pi$')
    P.style_axis(axes[2], f'c) Equilibrium check of the {probe} step',
                 '$q = (\\theta_\\mathrm{H}/\\theta_v)\\,/\\,(K\\,p/P^\\circ)^{1/n}$')
    handles = [plt.Line2D([], [], color=colors[c], label=c, **{**dot, 'ms': 7}) for c in cases]
    handles += [plt.Line2D([], [], color=P.INK, label='normal (b, c)', **dot),
                plt.Line2D([], [], mec=P.INK, label='switched off (b, c)', **ring)]
    fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=9)
    n_rxn = str(df['reactions'].iloc[0]).count('#') if 'reactions' in df and len(df) else 0
    failed = int((~df['converged'].astype(bool)).sum()) if 'converged' in df else 0
    fig.suptitle(f'What-if: {label} switched off' + (f' ({n_rxn} reaction{"s" if n_rxn != 1 else ""})' if n_rxn else '')
                 + (f'; {failed} of {len(df)} points not converged (no marker)' if failed else ''),
                 x=0.01, ha='left', fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    return fig


def top_reactions(ctx, case, T, species, n=5, P=None):
    """Largest sources and sinks of `species` at a main steady state (P in Pa; None = the
    reference pressure)."""
    model = ctx.model
    summary = ctx.main()['summary']
    Pb = (ctx.cfg.reference_pressure if P is None else float(P)) / 1e5
    hit = summary[(summary['basis'] == 'fixed_probe') & (summary['case'] == case) & (summary['T_K'] == float(T))
                  & np.isclose(summary['P_ref_bar'].astype(float), Pb)]
    if hit.empty or not hit['converged'].iloc[0]:
        return []
    res = ctx.main()['results'][hit['key'].iloc[0]]
    surf = common.load_state(model, res['row'], res['theta_a'])
    q = surf.net_rates_of_progress / model.gamma
    out = []
    for i, r in enumerate(surf.reactions()):
        nu = r.products.get(species, 0) - r.reactants.get(species, 0)
        if nu and q[i] != 0:
            note = r.input_data.get('note', '').splitlines()
            out.append({'case': case, 'T_K': float(T), 'P_ref_bar': Pb, 'species': species, 'label': species_label(model, species),
                        'coverage': res['theta_a'][model.species_names.index(species)],
                        'rate_1_per_s': q[i] * nu, 'reaction': i, 'equation': r.equation,
                        'source': note[1] if len(note) > 1 else ''})
    out.sort(key=lambda d: -abs(d['rate_1_per_s']))
    return out[:n]


def run(ctx):
    s = ctx.cfg.diagnostics
    ctx.log('\n== diagnostics ==')
    out = {}
    if s['adsorption_checks']:
        out['binding'] = binding_check(ctx)
        ctx.save_table(out['binding'], 'mechanism_binding_check.csv')
        fig = plot_binding_check(out['binding'])
        P.save(fig, ctx.path('binding_check.png'), ctx.cfg)
        plt.close(fig)
        ctx.log(out['binding'].drop(columns='literature_note').to_string(index=False, float_format=lambda v: f'{v:8.1f}'))
    tables = switch_off(ctx)
    if tables:
        for label, df in tables:
            ctx.save_table(df, f'mechanism_switch_off_{file_label(label)}.csv')
            fig = plot_switch_off(ctx.cfg, label, df, probe_label(ctx.model))
            P.save(fig, ctx.path(f'switch_off_{file_label(label)}.png'), ctx.cfg)
            plt.close(fig)
            ctx.log(df[['case', 'T_K', 'P_ref_bar', 'tof_ratio', 'slowdown_normal', 'slowdown_switched_off', 'q_normal',
                        'q_switched_off']].to_string(index=False, float_format=lambda v: f'{v:.3g}'))
        out['switch_off'] = pd.concat([df for _, df in tables], ignore_index=True)
        ctx.save_table(out['switch_off'], 'mechanism_switch_off.csv')
    else:
        ctx.log("  no switch-off what-if (diagnostics['switch_off'] = [])")
    if s['top_reactions']:
        rows = []
        for entry in s['top_reactions']:                 # (case, T, species) or (case, T, P, species)
            case, T, *rest = entry
            P_pa, sp = (rest[0], rest[1]) if len(rest) == 2 else (None, rest[0])
            rows += top_reactions(ctx, case, T, sp, s['top_n'], P_pa)
        out['top_reactions'] = pd.DataFrame(rows)
        ctx.save_table(out['top_reactions'], 'top_reactions.csv')
        for r in rows:
            ctx.log(f"  {r['label']} at {r['case']}, {r['T_K']:.0f} K: {'makes' if r['rate_1_per_s'] > 0 else 'uses '} "
                    f"{abs(r['rate_1_per_s']):9.3e} /s  #{r['reaction']:5d}  {r['equation']}")
    ctx.outputs['diagnostics'] = out
    return out
