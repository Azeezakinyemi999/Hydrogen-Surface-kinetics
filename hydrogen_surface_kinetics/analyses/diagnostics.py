"""Is the mechanism believable? (documentation: Stage 10)

1. adsorption_checks: reaction enthalpies from the mechanism's own thermo vs literature.
2. switch-off what-if: every reaction that makes or uses `switch_off` (default: the
   molecular adsorbate) gets rate multiplier 0, and the main cases are solved again.
3. top_reactions: largest sources and sinks of chosen species at chosen conditions.

Outputs: mechanism_binding_check.csv + binding_check.png, mechanism_switch_off.csv,
top_reactions.csv.
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .. import conditions as C
from .. import plotting as P
from .. import solver as S
from ..labels import species_label
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


def switch_off(ctx):
    model, s = ctx.model, ctx.cfg.diagnostics
    species = s['switch_off'] or model.molecular_adsorbate
    if not species:
        return None
    rxns = [i for i, r in enumerate(model.surf.reactions()) if species in r.reactants or species in r.products]
    multipliers = {i: 0.0 for i in rxns}
    cases = s['switch_off_cases'] or C.case_labels(ctx.cfg)
    main = ctx.main()['summary']
    wanted = {(pt['T'], round(pt['P'] / 1e5, 12)) for pt in common.points(ctx, 'diagnostics')}
    at_points = [(T, round(Pb, 12)) in wanted for T, Pb in zip(main['T_K'], main['P_ref_bar'].astype(float))]
    base = main[(main['basis'] == 'fixed_probe') & main['case'].isin(cases) & np.array(at_points, dtype=bool)]
    ctx.log(f'  switch-off what-if: {len(rxns)} reactions involve {species}')

    def job_func(job):
        theta, info = S.solve_steady_state(model, job['row']['T_K'], job['row']['P_total_bar'] * 1e5, job['row']['X'],
                                           S.clean_surface(model), multipliers=job['multipliers'])
        return {'theta': theta, 'info': info}

    results = ctx.main()['results']
    jobs = [common.job(f"switch_off {species} | {r.case} | {r.T_K:.0f} K | {r.P_ref_bar:g} bar",
                       results[r.key]['row'], multipliers=multipliers) for r in base.itertuples()]
    out = ctx.run_jobs(job_func, jobs, 'cache_diagnostics', describe=lambda res: S.describe(res['info']))
    rows = []
    for r, job in zip(base.itertuples(), jobs):
        res = out[job['key']]
        row = {'series': r.series, 'case': r.case, 'T_K': r.T_K, 'P_ref_bar': r.P_ref_bar,
               'tof_normal': r.tof_fwd, 'theta_free_normal': r.theta_free}
        ok = res['status'] == 'done' and res['info']['converged']
        row['converged'] = ok
        if ok:
            # rates with the switch-off still applied
            surf = common.load_state(model, job['row'], res['theta'], multipliers)
            row.update({'tof_switched_off': model.probe_rates(surf)[0] / model.gamma,
                        'theta_free_switched_off': res['theta'][model.species_names.index(model.site)],
                        'top_switched_off': common.top_species(model, res['theta'])})
        rows.append(row)
    df = pd.DataFrame(rows)
    df['tof_ratio'] = df.get('tof_switched_off', np.nan) / df['tof_normal']
    ref = df[df['case'] == ctx.cfg.reference_label].set_index(['T_K', 'P_ref_bar'])
    for col in ('tof_normal', 'tof_switched_off'):
        if col in df and len(ref):
            ref_val = [ref[col].get((T, Pb), np.nan) for T, Pb in zip(df['T_K'], df['P_ref_bar'])]
            df[f'slowdown_{col[4:]}'] = np.array(ref_val, dtype=float) / df[col]
    df.insert(0, 'switched_off', species)
    return df


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
    sw = switch_off(ctx)
    if sw is not None:
        out['switch_off'] = sw
        ctx.save_table(sw, 'mechanism_switch_off.csv')
        ctx.log(sw.to_string(index=False, float_format=lambda v: f'{v:.3g}'))
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
