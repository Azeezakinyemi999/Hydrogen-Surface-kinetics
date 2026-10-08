"""Shared plot helpers and the two output styles.

Plot functions draw in a plain style and call save(fig, path, cfg), which applies
  'screen'      : bolder lines and text, 300 dpi PNG (slides, quick viewing)
  'publication' : journal width, 7 pt regular text, thin lines, bold panel letters,
                  only condition titles ("300 K"), subscripted formulas,
                  vector PDF + 600 dpi PNG
"""
import os
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

SERIES_COLORS = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100', '#e87ba4', '#008300', '#4a3aa7', '#e34948']
EMPTY_COLOR, OTHER_COLOR = '#e4e3dd', '#8a8983'
INK, MUTED, GRID = '#33322f', '#8a8983', '#e4e3dd'
MARKERS = ['o', 's', 'D', '^', 'v', 'P', 'X', '*']
PUB_FONT_PT = 7
PUB_MAX_HEIGHT_CM = 22.0          # printable page height

mpl.rcParams.update({'figure.dpi': 150, 'savefig.dpi': 300, 'pdf.fonttype': 42, 'ps.fonttype': 42,
                     'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
                     'mathtext.default': 'regular'})
_FORMULA = re.compile(r'(?<=[A-Z])(\d+)|(?<=[A-Z][a-z])(\d+)')   # digits after an element symbol


def colors_for(items, first=None):
    """Fixed palette order (never cycled past its length without warning)."""
    out = {} if first is None else {first: INK}
    rest = [i for i in items if i != first]
    for i, item in enumerate(rest):
        out[item] = SERIES_COLORS[i % len(SERIES_COLORS)]
    return out


def style_axis(ax, title='', ylabel=''):
    ax.set_title(title, fontsize=11, loc='left', color=INK)
    ax.set_ylabel(ylabel, fontsize=10, color=INK)
    ax.spines[['top', 'right']].set_visible(False)
    ax.spines[['left', 'bottom']].set_color(MUTED)
    ax.tick_params(colors='#55544f', labelsize=9)
    ax.grid(axis='y', color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def two_states_handle():
    return plt.Line2D([], [], marker='o', ls='', ms=9, mfc='none', mec=INK)


# --- Styles ----------------------------------------------------------------------------
def _screen(fig):
    for ax in fig.axes:
        for line in ax.get_lines():
            line.set_linewidth(line.get_linewidth() * 1.6)
            line.set_markersize(line.get_markersize() * 1.35)
        for coll in ax.collections:
            if hasattr(coll, 'set_sizes') and len(coll.get_sizes()):
                coll.set_sizes(coll.get_sizes() * 1.8)
        for spine in ax.spines.values():
            spine.set_linewidth(1.4)
        ax.tick_params(width=1.4, length=5)
    for text in fig.findobj(mpl.text.Text):
        if text.get_text():
            text.set_fontsize(text.get_fontsize() * 1.12)
            text.set_fontweight('bold')


def sub(s):
    return s if '$' in s else _FORMULA.sub(lambda m: f'$_{{{m.group(0)}}}$', s)   # H2 -> H$_2$


def _publication(fig, width_cm):
    w, h = fig.get_size_inches()
    W = width_cm / 2.54
    H = h * W / w * 1.3                                      # final print size; a little taller for labels
    if H > PUB_MAX_HEIGHT_CM / 2.54:                         # tall, narrow figures: fit the page, keep aspect
        W, H = W * PUB_MAX_HEIGHT_CM / 2.54 / H, PUB_MAX_HEIGHT_CM / 2.54
    fig.set_size_inches(W, H)
    if fig._suptitle is not None:
        fig._suptitle.set_text('')                          # titles belong in the caption
    axes = [ax for ax in fig.axes if ax.get_visible()]
    for i, ax in enumerate(axes):
        title = re.sub(r'^\s*[a-z]\)\s*', '', ax.get_title(loc='left'))
        title = title.split(':')[0]                         # keep condition titles ("300 K", "CO, 300 K")
        title = title if re.search(r'\d\s*K\b', title) else ''
        letter = rf'$\bf{{({chr(97 + i)})}}$ ' if len(axes) > 1 else ''
        ax.set_title(letter + sub(title), loc='left', fontsize=PUB_FONT_PT + 1)
        for axis, set_ticks in ((ax.xaxis, ax.set_xticks), (ax.yaxis, ax.set_yticks)):   # category labels
            labels = axis.get_ticklabels()
            texts = [t.get_text() for t in labels]
            if texts and any(re.search('[A-Za-z]', t) for t in texts) and not any('$' in t for t in texts):
                set_ticks(axis.get_ticklocs(), [sub(t) for t in texts], rotation=labels[0].get_rotation(),
                          ha=labels[0].get_horizontalalignment())
        for line in ax.get_lines():
            line.set_linewidth(max(0.6, line.get_linewidth() * 0.6))
            line.set_markersize(line.get_markersize() * 0.65)
            line.set_markeredgewidth(min(line.get_markeredgewidth(), 0.6))
        for coll in ax.collections:
            if hasattr(coll, 'set_sizes') and len(coll.get_sizes()):
                coll.set_sizes(coll.get_sizes() * 0.4)
                coll.set_linewidths(0.5)
        for spine in ax.spines.values():
            spine.set_linewidth(0.6)
        ax.tick_params(width=0.6, length=3, labelsize=PUB_FONT_PT - 0.5)
    # Figure legends: re-create them inside the layout (above, or to the right) so nothing overlaps
    for leg in list(fig.legends):
        anchor = leg.get_bbox_to_anchor().transformed(fig.transFigure.inverted())
        loc = 'outside right center' if anchor.x0 >= 0.99 else 'outside upper center'
        handles, labels = leg.legend_handles, [t.get_text() for t in leg.get_texts()]
        title = leg.get_title().get_text()
        ncol = 1 if loc.endswith('right center') else max(1, len(labels))
        leg.remove()
        fig.legend(handles, [sub(l) for l in labels], loc=loc, ncol=ncol, frameon=False,
                   title=sub(title) or None, fontsize=PUB_FONT_PT, title_fontsize=PUB_FONT_PT)
    for text in fig.findobj(mpl.text.Text):
        t = text.get_text()
        if t:
            text.set_text(sub(t))
            if '$\\bf' not in t:
                text.set_fontsize(min(text.get_fontsize(), PUB_FONT_PT))
    fig.set_layout_engine('constrained')


def apply_style(fig, style='screen', width_cm=17.8):
    if getattr(fig, '_hsk_styled', False):
        return fig
    if style == 'publication':
        _publication(fig, width_cm)
    else:
        _screen(fig)
    fig._hsk_styled = True
    return fig


def save(fig, path, cfg=None, style=None):
    """Apply the run's plot style and save to `path` (.png); publication also writes a .pdf.
    Returns the list of files written."""
    style = style or (cfg.plot_style if cfg is not None else 'screen')
    width = cfg.publication_width_cm if cfg is not None else 17.8
    apply_style(fig, style, width)
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    stem = os.path.splitext(path)[0]
    written = []
    if style == 'publication':
        fig.savefig(stem + '.pdf', pad_inches=0.02, bbox_inches='tight')
        fig.savefig(stem + '.png', dpi=600, pad_inches=0.02, bbox_inches='tight')
        written += [stem + '.pdf', stem + '.png']
    else:
        fig.savefig(stem + '.png', dpi=300, bbox_inches='tight')
        written.append(stem + '.png')
    return written


def log_floor(values, floor=1e-12):
    return np.clip(np.asarray(values, dtype=float), floor, None)
