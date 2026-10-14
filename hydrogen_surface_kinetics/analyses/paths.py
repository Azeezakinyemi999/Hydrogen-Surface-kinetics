"""Reaction path diagrams: element fluxes between species at a solved steady state
(documentation: Stage 11e). No new solving. Near equilibrium the net flows cancel, so the
default is 'OneWayFlow' (both directions) with a small threshold.
Nodes can show structure images (paths['species_images']) with the readable name under each
(image_names); orientation 'TB' or 'LR'; the relative-flux numbers on the arrows can be hidden
(edge_labels = False).
Outputs: reaction_paths/path_<element>_<case>_<T>K_<P>bar_<flow>.dot (+ .png/.pdf if Graphviz is installed).
"""
import os
import re
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


IMAGE_INCH_PER_PX = 0.02            # structure images: inches per image pixel (RMG draws small PNGs)


def images_folder(ctx):
    """Where structure images are: paths['species_images'] (a folder; False = none), or by
    default species/ next to settings.py (RMG's own images), else output/flux_species/ (drawn
    by the flux analysis). Images are matched by the species' Cantera name: <name>.png."""
    setting = ctx.cfg.paths['species_images']
    if setting is False:
        return None
    if setting not in (None, True):
        p = os.path.expanduser(str(setting))
        return p if os.path.isabs(p) else os.path.normpath(os.path.join(ctx.cfg.run_dir, p))
    for cand in (os.path.join(ctx.cfg.run_dir, 'species'), os.path.join(ctx.cfg.output_dir, 'flux_species')):
        if os.path.isdir(cand):
            return cand
    return None


def display_name(model, name):
    """Readable name: surface species by species_label, gas species without the '(index)'."""
    if name in model.surf.species_names:
        return species_label(model, name)
    return re.sub(r'\(\d+\)$', '', name)


NODE_RE = re.compile(r'^(s\d+) \[ fontname="Helvetica", label="(.*?)"\];$', re.M)


EDGE_RE = re.compile(r'(s\d+ -> s\d+)\[fontname="Helvetica", penwidth=[\d.]+, arrowsize=[\d.]+, '
                     r'color="[^"]*"\n, label = " ?([^"]*)"\];')
RMG_MAX_EDGE, RMG_MAX_NODE, RMG_AMOUNT_TOL = 8.0, 6.0, 1e-6      # widths (pt) at the largest flux / amount


def log_width(x, tol, w_max, w_min=0.0):
    """RMG's line width: w_max at x = 1, falling linearly in log10(x) to 0 at x = tol."""
    if not x or x < tol:
        return 0.0
    return max(w_min, w_max * (1.0 - np.log10(x) / np.log10(tol)))


def style_dot(model, dot, images=None, image_names=True, orientation='TB', edge_labels=True,
              style='cantera', amounts=None, threshold=1e-8):
    """Restyle Cantera's diagram (nodes still labeled with Cantera names):
      images      : folder of <Cantera name>.png structure images (None = text nodes)
      image_names : print the readable name under each image
      orientation : 'TB' (top to bottom) or 'LR' (left to right)
      edge_labels : keep the relative-flux numbers on the arrows (widths always show flux)
      style       : 'cantera' (Cantera's colors and widths) or 'rmg' (RMG's flux-diagram look:
                    black arrows with width ~ log of the relative flux down to `threshold`, and
                    each species in an outline whose thickness ~ log of its amount, `amounts`:
                    coverage for surface species, mole fraction for gas species)"""
    rmg = style == 'rmg'

    def outline(name):
        if not rmg:
            return ''
        w = log_width((amounts or {}).get(name, 0.0), RMG_AMOUNT_TOL, RMG_MAX_NODE)
        return f'shape=ellipse, penwidth={w:.2f}, color="black", ' if w > 0 else 'shape=ellipse, color="gray80", '

    def node(m):
        nid, name = m.group(1), m.group(2)
        label = display_name(model, name)
        path = os.path.join(images, f'{name}.png') if images else None
        font = 'sans' if rmg else 'Helvetica'
        if not (path and os.path.exists(path)):
            return f'{nid} [ {outline(name)}fontname="{font}", label="{label}"];'
        from PIL import Image
        with Image.open(path) as im:
            w, h = (v * IMAGE_INCH_PER_PX for v in im.size)
        if not image_names and not rmg:
            return (f'{nid} [ shape=none, label="", image="{path}", imagescale=true, fixedsize=true, '
                    f'width={w:.3f}, height={h:.3f} ];')
        text = label.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
        name_row = f'<TR><TD><FONT FACE="{font}" POINT-SIZE="11">{text}</FONT></TD></TR>' if image_names else ''
        shape = outline(name) if rmg else 'shape=none, '
        return (f'{nid} [ {shape}margin=0, label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="0" '
                f'CELLPADDING="1"><TR><TD FIXEDSIZE="TRUE" WIDTH="{w * 72:.0f}" HEIGHT="{h * 72:.0f}">'
                f'<IMG SRC="{path}" SCALE="TRUE"/></TD></TR>{name_row}</TABLE>> ];')
    dot = NODE_RE.sub(node, dot)
    if rmg:                                                          # black arrows, width ~ log relative flux
        def edge(m):
            x = float(m.group(2))
            w = log_width(x, threshold, RMG_MAX_EDGE, w_min=0.6)
            label = f', label=" {m.group(2).strip()}", fontname="sans"' if edge_labels else ''
            return f'{m.group(1)}[penwidth={w:.2f}, arrowsize=0.8, color="black"{label}];'
        dot = EDGE_RE.sub(edge, dot)
        dot = dot.replace('fontname = "Helvetica";', 'fontname = "sans";')
    if orientation not in ('TB', 'LR'):
        raise ValueError("paths['orientation'] must be 'TB' or 'LR'")
    extra = ' overlap=false; fontname="sans"; fontsize=10;' if rmg else ''
    dot = dot.replace('digraph reaction_paths {', f'digraph reaction_paths {{\nrankdir="{orientation}"; '
                      f'nodesep=0.5; ranksep=0.7;{extra}', 1)
    # Cantera scales arrowheads up to 3x with the flux: cap them so the largest pair stays readable
    dot = re.sub(r'arrowsize=([\d.]+)', lambda m: f'arrowsize={min(float(m.group(1)), 1.4):.2f}', dot)
    if not edge_labels:
        dot = re.sub(r'\n, label = "[^"]*"', '', dot)               # edge labels only; the title stays
        dot = dot.replace('; labels = flux / max flux', '')
    if rmg:
        dot = dot.replace('width ~ flux', 'arrow width ~ log(flux / max flux), outline ~ log(coverage or mole fraction)')
    return dot


def raw_dot(ctx, case, T, Pa, element):
    """(dot text, amounts) for a solved main steady state, or None if it did not converge: Cantera's
    element-flux diagram (nodes named by Cantera name), and {species: coverage (surface) or mole
    fraction (gas)} at that state."""
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
    amounts = dict(zip(surf.species_names, surf.coverages))
    amounts.update(model.gas.mole_fraction_dict())
    return d.get_dot(), amounts


def render(dot, stem, formats, dpi=200):
    """Write stem.dot and, if Graphviz is installed, stem.<format> for each format."""
    with open(stem + '.dot', 'w') as f:
        f.write(dot)
    files = [stem + '.dot']
    if shutil.which('dot'):
        for fmt in formats:
            extra = [f'-Gdpi={dpi}'] if fmt == 'png' else []
            subprocess.run(['dot', f'-T{fmt}', *extra, stem + '.dot', '-o', f'{stem}.{fmt}'], check=True)
            files.append(f'{stem}.{fmt}')
    return files


def diagram(ctx, case, T, Pa, element, images=None):
    s = ctx.cfg.paths
    out = raw_dot(ctx, case, T, Pa, element)
    if out is None:
        return None
    dot = style_dot(ctx.model, out[0], images, s['image_names'], s['orientation'], s['edge_labels'],
                    s['style'], out[1], s['threshold'])
    stem = ctx.path('reaction_paths',
                    f"path_{element}_{case.replace(' ', '')}_{T:.0f}K_{Pa / 1e5:g}bar_{s['flow_type']}")
    return render(dot, stem, s['formats'])


def run(ctx):
    s = ctx.cfg.paths
    ctx.log('\n== paths: reaction path diagrams ==')
    if not shutil.which('dot'):
        ctx.log('  Graphviz "dot" not found: writing .dot files only')
    elements = s['elements'] or [el for el in ctx.model.gas.species(ctx.model.probe_gas).composition]
    conds = [normalize(c, ctx) for c in (s['conditions'] or default_conditions(ctx))]
    images = images_folder(ctx)
    ctx.log(f'  structure images: {images}' if images else '  structure images: none found (text labels)')
    files = {}
    for element in elements:
        for case, T, Pa in conds:
            out = diagram(ctx, case, T, Pa, element, images)
            if out:
                files[(element, case, T, Pa / 1e5)] = out
                ctx.log(f'  {os.path.relpath(out[-1], ctx.cfg.run_dir)}')
    ctx.outputs['paths'] = files
    return files
