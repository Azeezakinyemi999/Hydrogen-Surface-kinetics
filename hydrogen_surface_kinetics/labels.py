"""Readable names for species in tables and plots."""
ORDER = ['C', 'H', 'N', 'O']
SITE_ELEMENTS = {'X', 'Pt', 'Ni', 'Cu', 'Pd', 'Rh', 'Ru', 'Fe', 'Co', 'Ag', 'Au', 'Ir'}


def formula(composition, skip=SITE_ELEMENTS):
    comp = {el: int(round(n)) for el, n in composition.items() if el not in skip and n}
    els = [el for el in ORDER if el in comp] + sorted(set(comp) - set(ORDER))
    return ''.join(f"{el}{comp[el] if comp[el] > 1 else ''}" for el in els)


def probe_label(model):
    return model.cfg.probe_label or formula(model.gas.species(model.probe_gas).composition)


def species_label(model, name):
    """cfg.species_labels first, then role names, then formula + '*' per site + mechanism number."""
    labels = model.cfg.species_labels
    if name in labels:
        return labels[name]
    if name == model.site:
        return 'empty site'
    sp = model.surf.species(name)
    star = '*' * int(round(sp.size))
    f = formula(sp.composition)
    if name == model.probe_adsorbate:
        return f + star
    if name == model.molecular_adsorbate:
        return f'{f}{star} (molecular)'
    number = name.rsplit('(', 1)[-1].rstrip(')') if name.endswith(')') else ''
    return f'{f}{star}' + (f' #{number}' if number else '')
