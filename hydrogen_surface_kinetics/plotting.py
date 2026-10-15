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
SECONDARY_INK, RANGE_COLOR = '#52514e', '#c3c2b7'
DIVERGING_LOW, DIVERGING_HIGH = '#2a78d6', '#e34948'   # blue <-> red poles; neutral midpoint = ink
MARKERS = ['o', 's', 'D', '^', 'v', 'P', 'X', '*']
PUB_FONT_PT = 7
PUB_MAX_HEIGHT_CM = 22.0          # printable page height
MAX_PANELS, MAX_ROWS, MAX_COLS = 9, 3, 4   # per figure; more panels are split into name_1, name_2, ...

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


def pages(items, size):
    """`items` in the fewest groups of at most `size`, as even as possible (9 items, size 4:
    3 + 3 + 3)."""
    items = list(items)
    if not items:
        return []
    n = -(-len(items) // max(1, size))
    k = -(-len(items) // n)
    return [items[i:i + k] for i in range(0, len(items), k)]


def rows_per_page(ncol):
    """How many rows of `ncol` panels one figure holds."""
    return max(1, min(MAX_ROWS, MAX_PANELS // ncol))


def grid_shape(n, max_cols=MAX_COLS):
    """(rows, columns) for n <= MAX_PANELS panels: one row up to max_cols, then as square as
    the limits allow (5 -> 2 x 3, 7 -> 2 x 4, 9 -> 3 x 3)."""
    nrow = -(-n // max_cols)
    return nrow, -(-n // nrow)


def keep_title(ax, text):
    """A panel title that names the panel's condition (case, temperature): the publication
    style keeps it, while descriptive titles move to the caption."""
    ax.set_title(text, fontsize=11, loc='left', color=INK)
    ax._hsk_keep_title = True


HALF_MARK = dict(ls='', marker='D', ms=7, mfc='white', mew=1.5, zorder=4)   # 'rate halved here' marker


def half_mark_handle(label):
    return plt.Line2D([], [], **{**HALF_MARK, 'mec': INK}, label=label)


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


def panel_letter(i):
    """The label of panel i (from 0): a, b, ..., z, aa, ab, ..."""
    s = ''
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(97 + r) + s
    return s


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
        if not (getattr(ax, '_hsk_keep_title', False) or re.search(r'\d\s*K\b', title)):
            title = ''
        letter = rf'$\bf{{({panel_letter(i)})}}$ ' if len(axes) > 1 else ''
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
        entry_cm = 1.0 + 0.13 * max((len(re.sub(r'[${}_\\]', '', l)) for l in labels), default=0)  # marker + 7 pt text
        ncol = 1 if loc.endswith('right center') else max(1, min(len(labels), int(W * 2.54 // entry_cm)))
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


def save_parts(figs, path, cfg=None, style=None):
    """Save the parts of one figure and close them: a single part as `path`, several as
    path_1, path_2, ... Parts left from an earlier run with a different count are removed.
    Returns the list of files written."""
    figs = [f for f in figs if f is not None]
    if not figs:
        return []
    stem, ext = os.path.splitext(path)
    folder, base = os.path.split(stem)
    names = [path] if len(figs) == 1 else [f'{stem}_{i}{ext}' for i in range(1, len(figs) + 1)]
    keep = {os.path.splitext(os.path.basename(n))[0] for n in names}
    part = re.compile(re.escape(base) + r'(_\d+)?\.(png|pdf)')
    if os.path.isdir(folder or '.'):
        for f in os.listdir(folder or '.'):
            if part.fullmatch(f) and os.path.splitext(f)[0] not in keep:
                os.remove(os.path.join(folder, f))
    written = []
    for fig, name in zip(figs, names):
        written += save(fig, name, cfg, style)
        plt.close(fig)
    return written


def signed(x, digits=2):
    """'+1.00' / '−2.00' with a true minus sign."""
    return f'{x:+.{digits}f}'.replace('-', '\u2212')


def log_floor(values, floor=1e-12):
    return np.clip(np.asarray(values, dtype=float), floor, None)
