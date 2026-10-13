"""Reaction path diagrams: element fluxes between species at a solved steady state
(documentation: Stage 11e). No new solving. Near equilibrium the net flows cancel, so the
default is 'OneWayFlow' (both directions) with a small threshold.
Outputs: reaction_paths/path_<element>_<case>_<T>K_<P>bar_<flow>.dot (+ .png/.pdf if Graphviz is installed).
"""
import os
import shutil
import subprocess

import cantera as ct
import numpy as np

from ..conditions import case_labels
from ..labels import species_label
from . import common


def default_conditions(ctx):
    """Every case at every point of the paths section's series."""
    return [(case, pt['T'], pt['P']) for pt in common.points(ctx, 'paths') for case in case_labels(ctx.cfg)]


def normalize(cond, ctx):
    """(case, T) or (case, T, P in Pa) -> (case, T, P)."""
    case, T = cond[0], float(cond[1])
    P = float(cond[2]) if len(cond) > 2 else float(ctx.cfg.reference_pressure)
    return case, T, P


def diagram(ctx, case, T, Pa, element):
    model, s = ctx.model, ctx.cfg.paths
    summary = ctx.main()['summary']
    Pb = Pa / 1e5
    hit = summary[(summary['basis'] == 'fixed_probe') & (summary['case'] == case) & (summary['T_K'] == T)
                  & np.isclose(summary['P_ref_bar'].astype(float), Pb)]
    if hit.empty or not hit['converged'].iloc[0]:
        ctx.log(f'  {case}, {T:.0f} K, {Pb:g} bar: no converged main steady state, skipped')
        return None
    res = ctx.main()['results'][hit['key'].iloc[0]]
    surf = common.load_state(model, res['row'], res['theta_a'])
    d = ct.ReactionPathDiagram(surf, element)
    d.threshold = s['threshold']
    d.flow_type = s['flow_type']
    d.title = (f"{element} flux, {case}, {T:.0f} K, {Pb:g} bar  ({s['flow_type']}, width ~ flux, "
               f"threshold {s['threshold']:g} of max; labels = flux / max flux)")
    dot = d.get_dot()
    for name in surf.species_names:
        dot = dot.replace(f'label="{name}"', f'label="{species_label(model, name)}"')
    stem = ctx.path('reaction_paths',
                    f"path_{element}_{case.replace(' ', '')}_{T:.0f}K_{Pb:g}bar_{s['flow_type']}")
    with open(stem + '.dot', 'w') as f:
        f.write(dot)
    files = [stem + '.dot']
    if shutil.which('dot'):
        for fmt in s['formats']:
            extra = ['-Gdpi=150'] if fmt == 'png' else []
            subprocess.run(['dot', f'-T{fmt}', *extra, stem + '.dot', '-o', f'{stem}.{fmt}'], check=True)
            files.append(f'{stem}.{fmt}')
    return files


def run(ctx):
    s = ctx.cfg.paths
    ctx.log('\n== paths: reaction path diagrams ==')
    if not shutil.which('dot'):
        ctx.log('  Graphviz "dot" not found: writing .dot files only')
    elements = s['elements'] or [el for el in ctx.model.gas.species(ctx.model.probe_gas).composition]
    conds = [normalize(c, ctx) for c in (s['conditions'] or default_conditions(ctx))]
    files = {}
    for element in elements:
        for case, T, Pa in conds:
            out = diagram(ctx, case, T, Pa, element)
            if out:
                files[(element, case, T, Pa / 1e5)] = out
                ctx.log(f'  {os.path.relpath(out[-1], ctx.cfg.run_dir)}')
    ctx.outputs['paths'] = files
    return files
