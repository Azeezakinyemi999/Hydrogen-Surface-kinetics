"""Build workflow_explained.pptx: what each step of the workflow does, and why (method only;
running example runs/try_2/Cu).

    ~/anaconda3/bin/python3 build_deck.py            # build the deck next to this file
    ~/anaconda3/bin/python3 render.py workflow_explained.pptx <png folder>   # check every slide

python-pptx is in the anaconda base Python, not in rmg4.0. One module per section.
"""
import importlib
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from deck_lib import Deck  # noqa: E402

MODULES = ['slides_a', 'slides_b', 'slides_c', 'slides_d', 'slides_e', 'slides_f', 'slides_g', 'slides_h',
           'slides_i', 'slides_j', 'slides_k', 'slides_l', 'slides_m']


def build(path=os.path.join(HERE, 'workflow_explained.pptx')):
    d = Deck()
    for name in MODULES:
        try:
            mod = importlib.import_module(name)
        except ModuleNotFoundError as e:
            if e.name == name:                  # sections not written yet
                continue
            raise
        for f in mod.SLIDES:
            f(d)
    d.save(path)
    print(f'{d.count} slides -> {path}')
    return path


if __name__ == '__main__':
    build()
