"""Flux-diagram videos of the RMG input file's own batch reactors (documentation: Stage 2).

Optional: needs RMG-Py (with RDKit), Graphviz and ffmpeg. Every surfaceReactor block in
the RMG input file is run as RMG does (closed, isothermal, constant-volume batch reactor)
with Cantera, then drawn with rmgpy.tools.fluxdiagram.generate_flux_diagram.

Species are matched between the seed mechanism and the Cantera file by label (without
the "(index)" suffix), falling back to structure isomorphism where the labels differ:
the Cantera name is a SMILES with [Pt] for the site, rebuilt as an RMG Molecule.
Output: flux/<reactor number>/ (frames, .dot files, flux_diagram.avi).
"""
import glob
import os
import re
import shutil
import subprocess
import sys
import time

import cantera as ct
import numpy as np


STRIP_INDEX = re.compile(r'\(\d+\)$')
LABEL_RE = re.compile(r'label="t = ([^"]*)"')


def strip_index(name):
    return STRIP_INDEX.sub('', name)


def _path(ctx, p):
    if p is None:
        return None
    p = os.path.expanduser(str(p))
    return p if os.path.isabs(p) else os.path.normpath(os.path.join(ctx.cfg.run_dir, p))


def import_rmg(rmg_path):
    """Import the RMG-Py pieces used here; a clear message if they are missing."""
    if rmg_path and rmg_path not in sys.path:
        sys.path.insert(0, rmg_path)
    try:
        from rdkit import Chem
        from rmgpy.data.kinetics.database import KineticsDatabase
        from rmgpy.molecule import Atom, Bond, Molecule
        from rmgpy.molecule.element import get_element
        from rmgpy.rmg.main import RMG
        from rmgpy.rmg.model import CoreEdgeReactionModel
        from rmgpy.tools.fluxdiagram import generate_flux_diagram
    except ImportError as e:
        raise RuntimeError(f"the 'flux' analysis needs RMG-Py and RDKit ({e}); set flux['rmg_path'] "
                           'or run with the environment that has the compiled RMG-Py') from e
    return dict(Chem=Chem, KineticsDatabase=KineticsDatabase, Atom=Atom, Bond=Bond, Molecule=Molecule,
                get_element=get_element, RMG=RMG, CoreEdgeReactionModel=CoreEdgeReactionModel,
                generate_flux_diagram=generate_flux_diagram)


# --- species and reaction matching -------------------------------------------------------
def load_seed(rmg, seed_dir):
    kdb = rmg['KineticsDatabase']()
    name = os.path.basename(os.path.normpath(seed_dir))
    kdb.load_libraries(path=os.path.dirname(os.path.normpath(seed_dir)), libraries=[seed_dir])
    lib = kdb.libraries[name]
    species_dict = lib.get_species(os.path.join(seed_dir, 'dictionary.txt'))
    return species_dict, [entry.item for entry in lib.entries.values()]


def mol_from_pt_smiles(rmg, smiles):
    """RMG writes Cantera names as SMILES with [Pt] for the site; rebuild an RMG Molecule
    (neither from_smiles nor from_rdkit_mol maps Pt back to the X site atom)."""
    rd = rmg['Chem'].MolFromSmiles(smiles)
    if rd is None:
        return None
    rd = rmg['Chem'].AddHs(rd)
    mol, atoms = rmg['Molecule'](), []
    for a in rd.GetAtoms():
        symbol = 'X' if a.GetSymbol() == 'Pt' else a.GetSymbol()
        radicals = 0 if symbol == 'X' else a.GetNumRadicalElectrons()   # bare [Pt] gets spurious radicals
        atom = rmg['Atom'](element=rmg['get_element'](symbol), radical_electrons=radicals,
                           charge=a.GetFormalCharge(), lone_pairs=0)
        mol.add_atom(atom)
        atoms.append(atom)
    for b in rd.GetBonds():
        mol.add_bond(rmg['Bond'](atoms[b.GetBeginAtomIdx()], atoms[b.GetEndAtomIdx()], order=b.GetBondTypeAsDouble()))
    mol.update_lone_pairs()
    mol.update_multiplicity()
    mol.update_atomtypes(raise_exception=False)
    return mol


def build_matched_model(rmg, species_dict, all_reactions, gas, surf, log=print):
    """Seed species <-> Cantera species (label, then structure). Matched seed species are
    relabeled to the Cantera name; reactions are re-pointed at the canonical species
    objects (generate_flux_diagram looks them up by identity)."""
    ct_name_by_base = {strip_index(n): ('gas', n) for n in gas.species_names}
    ct_name_by_base.update({strip_index(n): ('surface', n) for n in surf.species_names})
    seed_to_ct = {label: ct_name_by_base[label] for label in species_dict if label in ct_name_by_base}
    unmatched_seed = [label for label in species_dict if label not in seed_to_ct]
    taken = {strip_index(n) for _, n in seed_to_ct.values()}
    for base, (phase, ctname) in ct_name_by_base.items():
        if base in taken or not unmatched_seed:
            continue
        try:
            mol = mol_from_pt_smiles(rmg, base)
        except Exception:
            mol = None
        if mol is None:
            continue
        for label in unmatched_seed:
            if species_dict[label].is_isomorphic(mol):
                seed_to_ct[label] = (phase, ctname)
                unmatched_seed.remove(label)
                log(f'  structure match: seed {label} -> cantera {ctname}')
                break

    reactions = []
    for rxn in all_reactions:
        if all(s.label in seed_to_ct for s in rxn.reactants + rxn.products):
            rxn.reactants = [species_dict[s.label] for s in rxn.reactants]
            rxn.products = [species_dict[s.label] for s in rxn.products]
            rxn.pairs = None
            rxn.generate_pairs()
            reactions.append(rxn)

    matched, species = {}, []
    for label, (phase, ctname) in seed_to_ct.items():
        spec = species_dict[label]
        spec.label = strip_index(ctname)          # str(spec) == the Cantera name == image file name
        m = re.search(r'\((\d+)\)$', ctname)
        spec.index = int(m.group(1)) if m else -1
        matched[spec.label] = (phase, ctname)
        species.append(spec)
    log(f'  matched species: {len(species)} / {len(species_dict)} seed species; '
        f'seed species with no Cantera match: {unmatched_seed}')
    log(f'  Cantera species not in the seed (left out): {len(set(ct_name_by_base) - set(matched))}; '
        f'matched reactions: {len(reactions)} / {len(all_reactions)}')
    return species, reactions, matched, ct_name_by_base


def reaction_index(gas, surf):
    """(reactant labels, product labels) -> (phase, Cantera index, sign)."""
    index = {}
    for phase, rxns in (('gas', gas.reactions()), ('surface', surf.reactions())):
        for i, r in enumerate(rxns):
            reac = frozenset(strip_index(n) for n in r.reactants)
            prod = frozenset(strip_index(n) for n in r.products)
            index[(reac, prod)] = (phase, i, 1.0)
            index.setdefault((prod, reac), (phase, i, -1.0))
    return index


def draw_species_images(species, species_dir, log=print):
    os.makedirs(species_dir, exist_ok=True)
    n = 0
    for spec in species:
        path = os.path.join(species_dir, f'{spec}.png')
        if not os.path.exists(path):
            try:
                spec.molecule[0].draw(path)
                n += 1
            except Exception as e:
                log(f'  could not draw {spec}: {e}')
    log(f'  drew {n} species images into {species_dir}')


# --- one reactor -------------------------------------------------------------------------
def simulate(rs, gas, surf, ct_name_by_base, n_points, rtol, atol):
    """RMG surfaceReactor semantics: closed, isothermal, constant-volume batch reactor with
    V from 1 mol of gas (RMG's own convention) and A = V * surface_volume_ratio."""
    T, P = rs.T.value_si, rs.P_initial.value_si
    V = ct.gas_constant / 1e3 * T / P
    term_time = next((t.time.value_si for t in rs.termination if hasattr(t, 'time')), 1000.0)
    X = {ct_name_by_base[strip_index(s.label)][1]: x for s, x in rs.initial_gas_mole_fractions.items() if x}
    cov = {ct_name_by_base[strip_index(s.label)][1]: x for s, x in rs.initial_surface_coverages.items() if x}
    gas.TPX = T, P, X
    reactor = ct.IdealGasReactor(gas, energy='off')
    reactor.volume = V
    surf.TP = T, P
    surf.coverages = cov
    ct.ReactorSurface(surf, reactor, A=V * rs.surface_volume_ratio.value_si)
    net = ct.ReactorNet([reactor])
    net.rtol, net.atol = rtol, atol
    hist = {k: [] for k in ('t', 'gas_c', 'surf_c', 'gas_q', 'surf_q')}
    for t in np.logspace(-8, np.log10(term_time), n_points):
        net.advance(t)
        hist['t'].append(net.time)
        hist['gas_c'].append(reactor.thermo.concentrations.copy())
        hist['surf_c'].append(surf.concentrations.copy())
        hist['gas_q'].append(gas.net_rates_of_progress.copy())
        hist['surf_q'].append(surf.net_rates_of_progress.copy())
    return {k: np.array(v) for k, v in hist.items()}


def condition_label(rs):
    comp = ' '.join(f'{strip_index(s.label)}={x:.2f}' for s, x in rs.initial_gas_mole_fractions.items() if x > 0)
    return f'T={rs.T.value_si:.0f}K, P={rs.P_initial.value_si / 1e5:.1f}bar\\n{comp}'


def bake_condition_label(folder, label):
    """generate_flux_diagram labels each frame with the time only: add the condition,
    re-render the frames and re-stitch the video the same way it does."""
    for dot_path in sorted(glob.glob(os.path.join(folder, '*.dot'))):
        with open(dot_path) as f:
            content = f.read()
        content, n = LABEL_RE.subn(lambda m: f'label="{label}\\nt = {m.group(1)}"', content, count=1)
        if n:
            with open(dot_path, 'w') as f:
                f.write(content)
            subprocess.run(['dot', dot_path, '-Tpng', '-o', dot_path[:-4] + '.png'], check=True)
    avi = os.path.join(folder, 'flux_diagram.avi')
    if os.path.exists(avi):
        os.remove(avi)
    subprocess.run(['ffmpeg', '-framerate', '6', '-i', 'flux_diagram_%04d.png', '-c:v', 'mpeg4', '-r', '30',
                    '-pix_fmt', 'yuv420p', 'flux_diagram.avi'], check=True, cwd=folder, capture_output=True)


# --- the analysis ------------------------------------------------------------------------
def run(ctx):
    s, cfg, log = ctx.cfg.flux, ctx.cfg, ctx.log
    log('\n== flux: flux-diagram videos of the RMG reactors ==')
    missing = [k for k in ('rmg_input', 'seed_dir') if not s[k]]
    if missing:
        raise ValueError(f"flux: set {missing} in the flux section of the input file")
    tools = [t for t in ('dot', 'ffmpeg') if not shutil.which(t)]
    if tools:
        raise RuntimeError(f'flux: {tools} not found on PATH')
    rmg = import_rmg(_path(ctx, s['rmg_path']))
    rmg_input, seed_dir = _path(ctx, s['rmg_input']), _path(ctx, s['seed_dir'])
    species_dir = _path(ctx, s['species_dir']) or ctx.path('flux_species', '')
    out_root = ctx.path('flux', '')

    species_dict, all_reactions = load_seed(rmg, seed_dir)
    gas, surf = ctx.model.gas, ctx.model.new_surface()
    species, reactions, matched, ct_name_by_base = build_matched_model(rmg, species_dict, all_reactions, gas, surf, log)
    draw_species_images(species, species_dir, log)
    index = reaction_index(gas, surf)
    mapping = [index.get((frozenset(strip_index(x.label) for x in r.reactants),
                          frozenset(strip_index(x.label) for x in r.products))) for r in reactions]
    log(f'  reaction rates mapped: {sum(m is not None for m in mapping)} / {len(reactions)}')
    model = rmg['CoreEdgeReactionModel']()
    model.core.species, model.core.reactions = species, reactions
    pos = {spec.label: i for i, spec in enumerate(species)}
    rmg_to_ct = {pos[label]: v for label, v in matched.items()}

    job = rmg['RMG'](input_file=rmg_input)
    job.load_input(rmg_input)
    systems = job.reaction_systems
    indices = s['reactor_indices'] if s['reactor_indices'] is not None else list(range(len(systems)))
    log(f'  {len(systems)} reactors in {rmg_input}; running {len(indices)}')

    gas_idx = {n: i for i, n in enumerate(gas.species_names)}
    surf_idx = {n: i for i, n in enumerate(surf.species_names)}

    def one(jb):
        idx = jb['index']
        rs = systems[idx]
        out_dir = os.path.join(out_root, str(idx + 1))
        if s['skip_existing'] and os.path.exists(os.path.join(out_dir, 'flux_diagram.avi')):
            return {'skipped': True, 'folder': out_dir}
        h = simulate(rs, gas, ctx.model.new_surface(), ct_name_by_base, s['n_points'], s['rtol'], s['atol'])
        conc = np.zeros((len(h['t']), len(species)))
        for p, (phase, name) in rmg_to_ct.items():
            conc[:, p] = h['gas_c'][:, gas_idx[name]] if phase == 'gas' else h['surf_c'][:, surf_idx[name]]
        rates = np.zeros((len(h['t']), len(reactions)))
        for i, m in enumerate(mapping):
            if m is not None:
                phase, ci, sign = m
                rates[:, i] = sign * (h['gas_q'][:, ci] if phase == 'gas' else h['surf_q'][:, ci])
        os.makedirs(out_dir, exist_ok=True)
        rmg['generate_flux_diagram'](model, h['t'], conc, rates, out_dir, species_directory=species_dir)
        label = condition_label(rs)
        bake_condition_label(out_dir, label)
        return {'skipped': False, 'folder': out_dir, 'label': label}

    # One reactor at a time in this process: forking after RDKit/Cairo have drawn images
    # crashes the child on macOS, and each reactor is a single long Cantera run anyway.
    res = {}
    for i in indices:
        key, t0 = f'flux reactor {i + 1}', time.time()
        try:
            res[key] = {'status': 'done', **one({'index': i})}
            note = 'skipped (video exists)' if res[key]['skipped'] else os.path.relpath(res[key]['folder'], cfg.run_dir)
        except Exception as e:                      # keep going: one stiff reactor must not stop the batch
            res[key] = {'status': 'error', 'error': f'{type(e).__name__}: {str(e).strip().splitlines()[0][:200]}'}
            note = 'FAILED - ' + res[key]['error']
        log(f'  {key}: {time.time() - t0:5.0f} s  {note}')
    failed = [k for k, r in res.items() if r['status'] != 'done']
    log(f"  failed: {failed or 'none'}")
    ctx.outputs['flux'] = res
    return res
