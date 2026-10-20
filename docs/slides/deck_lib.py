"""Drawing helpers for the workflow deck: native PowerPoint shapes only (python-pptx).
(From Surface-Impurity-Effect/docs/slides/deck_lib.py; here with O2 and N2 colours and a Cu(111) slab.)

Every picture is built from text boxes, rectangles, circles, freeform lines and connectors,
so each box, arrow and label stays editable in PowerPoint.

Text markup used by text(), shape() and actions():
    _{x}   subscript          ^{x}   superscript
    **x**  bold               `x`    code font (Courier New)
    [[KEY|x]]  bold text in colour KEY (a key of COL, e.g. [[CO|CO]])
"""
import math
import re

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_LINE
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, MSO_AUTO_SIZE, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

W, H = 13.333, 7.5          # slide size, inches (16:9)
LM = 0.55                   # side margin
TITLE, BODY, CODE, MATH = 'Cambria', 'Calibri', 'Courier New', 'Cambria'

COL = dict(
    # neutrals
    ink='1F2328', muted='4E5661', faint='8C939D', rule='D0D4DA', panel='F2F3F5', panel2='E6E9ED',
    white='FFFFFF', dark='3A4048',
    # one colour per gas species (and everything made from it)
    H='2B6CB0', CO='B83246', H2O='6B4FA0', NH3='C0508F', CH4='8A5A2B', CO2='4A5568', O2='00838F', N2='3949AB',
    # meaning
    ok='2F855A', okbg='E7F4EC', amber='B7791F', amberbg='FFF4DB',
    # surface and code
    ni='CDD2D9', niline='A7AEB8', codebg='F6F8FA', codeline='D0D7DE',
)
SPECIES_TEXT = dict(H='white', CO='white', H2O='white', NH3='white', CH4='white', CO2='white', O2='white', N2='white')


def rgb(c):
    return RGBColor.from_string(COL.get(c, c).lstrip('#').upper())


def E(v):
    """inches -> EMU (int)"""
    return int(Inches(v))


# --- text ------------------------------------------------------------------------------
_TOK = re.compile(r'(\*\*.+?\*\*|`[^`]+`|_\{[^}]*\}|\^\{[^}]*\}|\[\[[A-Za-z0-9]+\|[^\]]+\]\])')


def _runs(p, s, size, color, font, bold=False, italic=False, base=None):
    for part in _TOK.split(s):
        if not part:
            continue
        if part.startswith('**') and part.endswith('**') and len(part) > 4:
            _runs(p, part[2:-2], size, color, font, True, italic, base)
            continue
        if part.startswith('[['):
            key, inner = part[2:-2].split('|', 1)
            _runs(p, inner, size, key, font, True, italic, base)
            continue
        txt, f, b, bs = part, font, bold, base
        if part.startswith('`') and part.endswith('`') and len(part) > 2:
            txt, f = part[1:-1], CODE
        elif part.startswith('_{'):
            txt, bs = part[2:-1], -25000
        elif part.startswith('^{'):
            txt, bs = part[2:-1], 30000
        r = p.add_run()
        r.text = txt
        r.font.size = Pt(size)
        r.font.name = f
        r.font.bold = b
        r.font.italic = italic
        r.font.color.rgb = rgb(color)
        if bs:
            r._r.get_or_add_rPr().set('baseline', str(bs))


def _fill(tf, content, size=14, color='ink', font=BODY, bold=False, italic=False, align='l',
          anchor='t', wrap=True, spacing=None, margin=0.03, space_after=0):
    tf.word_wrap = wrap
    tf.auto_size = MSO_AUTO_SIZE.NONE
    tf.margin_left = tf.margin_right = Inches(margin)
    tf.margin_top = tf.margin_bottom = Inches(min(margin, 0.03))
    tf.vertical_anchor = {'t': MSO_ANCHOR.TOP, 'm': MSO_ANCHOR.MIDDLE, 'b': MSO_ANCHOR.BOTTOM}[anchor]
    paras = content if isinstance(content, (list, tuple)) else [content]
    for i, para in enumerate(paras):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        markup, ov = (para if isinstance(para, tuple) else (para, {}))
        p.alignment = {'l': PP_ALIGN.LEFT, 'c': PP_ALIGN.CENTER, 'r': PP_ALIGN.RIGHT}[ov.get('align', align)]
        if spacing:
            p.line_spacing = spacing
        if space_after or ov.get('space_after'):
            p.space_after = Pt(ov.get('space_after', space_after))
        if ov.get('space_before'):
            p.space_before = Pt(ov['space_before'])
        _runs(p, markup, ov.get('size', size), ov.get('color', color), ov.get('font', font),
              ov.get('bold', bold), ov.get('italic', italic))


def text(sl, x, y, w, h, content, size=14, color='ink', font=BODY, bold=False, italic=False,
         align='l', anchor='t', wrap=True, spacing=None, margin=0.03, space_after=0):
    tb = sl.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    _fill(tb.text_frame, content, size, color, font, bold, italic, align, anchor, wrap, spacing,
          margin, space_after)
    return tb


# --- shapes ----------------------------------------------------------------------------
_KINDS = {'rect': MSO_SHAPE.RECTANGLE, 'rrect': MSO_SHAPE.ROUNDED_RECTANGLE, 'oval': MSO_SHAPE.OVAL,
          'diamond': MSO_SHAPE.DIAMOND, 'hex': MSO_SHAPE.HEXAGON, 'chevron': MSO_SHAPE.CHEVRON,
          'pentagon': MSO_SHAPE.PENTAGON, 'doc': MSO_SHAPE.FLOWCHART_DOCUMENT,
          'can': MSO_SHAPE.CAN, 'tri': MSO_SHAPE.ISOSCELES_TRIANGLE}


def shape(sl, kind, x, y, w, h, fill='white', line='faint', lw=1.0, dash=False, radius=None,
          content=None, size=14, color='ink', font=BODY, bold=False, italic=False, align='c',
          anchor='m', margin=0.05, spacing=None, rotation=None, wrap=True):
    s = sl.shapes.add_shape(_KINDS[kind], Inches(x), Inches(y), Inches(w), Inches(h))
    s.shadow.inherit = False
    if fill is None:
        s.fill.background()
    else:
        s.fill.solid()
        s.fill.fore_color.rgb = rgb(fill)
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = rgb(line)
        s.line.width = Pt(lw)
        if dash:
            s.line.dash_style = MSO_LINE.DASH
    if kind == 'rrect':
        s.adjustments[0] = 0.5 if radius == 'pill' else min(0.5, (radius if radius is not None else 0.1) / max(min(w, h), 1e-6))
    if rotation is not None:
        s.rotation = rotation
    _fill(s.text_frame, content if content is not None else '', size, color, font, bold, italic,
          align, anchor, wrap, spacing, margin)
    return s


def circle(sl, cx, cy, d, fill='white', line='faint', lw=1.0, content=None, size=14, color='ink',
           font=BODY, bold=False, dash=False):
    return shape(sl, 'oval', cx - d / 2, cy - d / 2, d, d, fill=fill, line=line, lw=lw, dash=dash,
                 content=content, size=size, color=color, font=font, bold=bold, margin=0.0, wrap=False)


def _line_end(ln, which, kind='triangle', w='med', length='med'):
    tag = qn('a:headEnd') if which == 'head' else qn('a:tailEnd')
    el = ln.find(tag)
    if el is None:
        el = etree.SubElement(ln, tag)
    el.set('type', kind)
    el.set('w', w)
    el.set('len', length)


def arrow(sl, x1, y1, x2, y2, color='muted', lw=1.5, end=True, start=False, dash=False, head='med'):
    c = sl.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1), Inches(x2), Inches(y2))
    c.shadow.inherit = False            # the template's connector style adds a faint shadow (effectRef 1)
    c.line.color.rgb = rgb(color)
    c.line.width = Pt(lw)
    if dash:
        c.line.dash_style = MSO_LINE.DASH
    ln = c.line._get_or_add_ln()
    if start:
        _line_end(ln, 'head', w=head, length=head)
    if end:
        _line_end(ln, 'tail', w=head, length=head)
    return c


def line(sl, x1, y1, x2, y2, color='rule', lw=1.0, dash=False):
    return arrow(sl, x1, y1, x2, y2, color=color, lw=lw, end=False, dash=dash)


def poly(sl, pts, color='muted', lw=1.5, end=False, start=False, dash=False, closed=False, fill=None,
         head='med'):
    fb = sl.shapes.build_freeform(E(pts[0][0]), E(pts[0][1]), scale=1.0)
    fb.add_line_segments([(E(x), E(y)) for x, y in pts[1:]], close=closed)
    s = fb.convert_to_shape()
    s.shadow.inherit = False
    if fill:
        s.fill.solid()
        s.fill.fore_color.rgb = rgb(fill)
    else:
        s.fill.background()
    if color is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = rgb(color)
        s.line.width = Pt(lw)
        if dash:
            s.line.dash_style = MSO_LINE.DASH
        ln = s.line._get_or_add_ln()
        if start:
            _line_end(ln, 'head', w=head, length=head)
        if end:
            _line_end(ln, 'tail', w=head, length=head)
    return s


def curve(sl, f, x0, x1, X, Y, Wd, Ht, n=60, **kw):
    """Plot y = f(x) for x in [x0, x1] into the box (X, Y, Wd, Ht); f returns 0..1 (0 = bottom)."""
    pts = []
    for i in range(n + 1):
        t = x0 + (x1 - x0) * i / n
        pts.append((X + Wd * (t - x0) / (x1 - x0), Y + Ht * (1 - f(t))))
    return poly(sl, pts, **kw)


# --- deck ------------------------------------------------------------------------------
class Deck:
    def __init__(self):
        self.prs = Presentation()
        self.prs.slide_width, self.prs.slide_height = Inches(W), Inches(H)
        self.layout = self.prs.slide_layouts[6]
        self.count = 0

    def slide(self, title, subtitle=None, part=None, tag=None, number=True):
        sl = self.prs.slides.add_slide(self.layout)
        self.count += 1
        bg = sl.background.fill
        bg.solid()
        bg.fore_color.rgb = rgb('white')
        if part:
            text(sl, LM, 0.13, 8.0, 0.3, part.upper(), size=12, color='faint', bold=True)
        if tag:
            tag_box(sl, tag)
        size = 26 if len(title) <= 64 else (24 if len(title) <= 70 else 22)
        text(sl, LM, 0.40, W - 2 * LM, 0.62, title, size=size, color='ink', font=TITLE, anchor='m')
        if subtitle:
            text(sl, LM, 1.02, W - 2 * LM, 0.4, subtitle, size=17, color='muted')
        shape(sl, 'rect', LM, 1.48, W - 2 * LM, 0.012, fill='rule', line=None)
        if number:
            text(sl, W - LM - 0.8, 7.17, 0.8, 0.26, str(self.count), size=10, color='faint', align='r')
        return sl

    def save(self, path):
        self.prs.save(path)


def tag_box(sl, label):
    w = 0.22 + 0.083 * len(label)
    shape(sl, 'rrect', W - LM - w, 0.11, w, 0.3, fill='white', line='rule', lw=0.75, radius='pill',
          content=label, size=12, color='muted', margin=0.02)


def notes(sl, lines):
    sl.notes_slide.notes_text_frame.text = '\n'.join(lines)


def actions(sl, items, y=6.04, h=1.0):
    n = len(items)
    gap = 0.14
    bw = (W - 2 * LM - gap * (n - 1)) / n
    for i, t in enumerate(items):
        x = LM + i * (bw + gap)
        shape(sl, 'rrect', x, y, bw, h, fill='panel', line=None, radius=0.08)
        circle(sl, x + 0.3, y + h / 2, 0.38, fill='dark', line=None, content=str(i + 1), size=14,
               color='white', bold=True)
        text(sl, x + 0.56, y + 0.03, bw - 0.62, h - 0.06, t, size=14, color='ink', anchor='m')


def caveat(sl, x, y, w, h, content, size=14):
    shape(sl, 'rrect', x, y, w, h, fill='amberbg', line='amber', lw=1.0, radius=0.08)
    circle(sl, x + 0.24, y + 0.24, 0.28, fill='amber', line=None, content='!', size=14, color='white', bold=True)
    text(sl, x + 0.46, y + 0.04, w - 0.52, h - 0.08, content, size=size, color='ink', anchor='m')


def good(sl, x, y, w, h, content, size=14):
    shape(sl, 'rrect', x, y, w, h, fill='okbg', line='ok', lw=1.0, radius=0.08,
          content=content, size=size, color='ink')


def schematic(sl, x, y, w=1.4, label='schematic'):
    text(sl, x, y, w, 0.26, label, size=12, color='faint', italic=True)


def label(sl, x, y, w, content, size=14, color='muted', align='l', h=0.3, bold=False, italic=False,
          font=BODY, anchor='t'):
    return text(sl, x, y, w, h, content, size=size, color=color, align=align, bold=bold, italic=italic,
                font=font, anchor=anchor)


def ref(sl, x, y, w, content, align='l', color='faint'):
    """A code location label, e.g. solver.py:146"""
    return text(sl, x, y, w, 0.26, content, size=12, color=color, font=CODE, align=align)


def code(sl, x, y, w, h, lines, size=14, first=None):
    """Code box. lines: list of strings; first: line number of the first line (shown faint)."""
    shape(sl, 'rrect', x, y, w, h, fill='codebg', line='codeline', lw=0.75, radius=0.06)
    paras = []
    for i, s in enumerate(lines):
        if first is not None:
            paras.append((f'{first + i:>3}  ' + s, {}))
        else:
            paras.append((s, {}))
    tb = text(sl, x + 0.1, y + 0.06, w - 0.2, h - 0.12, [p[0] for p in paras], size=size, font=CODE,
              color='ink', spacing=1.0)
    return tb


def step(sl, x, y, w, h, content, where=None, fill='white', line='faint', size=14, lw=1.0, align='c',
         dash=False):
    """Flow step: a box with optional code location under it."""
    s = shape(sl, 'rrect', x, y, w, h, fill=fill, line=line, lw=lw, radius=0.08, content=content,
              size=size, align=align, dash=dash)
    if where:
        ref(sl, x, y + h + 0.01, w, where, align='r')
    return s


def gate(sl, x, y, w, h, content, size=14, where=None):
    s = shape(sl, 'diamond', x, y, w, h, fill='white', line='dark', lw=1.25, content=content,
              size=size, margin=0.02)
    if where:
        ref(sl, x, y + h + 0.0, w, where, align='c')
    return s


def yes_no(sl, x, y, word, color=None, w=0.5):
    col = color or ('ok' if word in ('yes', 'agree', 'pass') else 'muted')
    text(sl, x, y, w, 0.28, word, size=12, color=col, bold=True)


def badge(sl, x, y, n, d=0.34, fill='dark'):
    circle(sl, x, y, d, fill=fill, line=None, content=str(n), size=13, color='white', bold=True)


# --- the surface -----------------------------------------------------------------------
def site(sl, cx, cy, d=0.44, lw=1.25):
    return circle(sl, cx, cy, d, fill='white', line='faint', lw=lw)


def adsorbate(sl, cx, cy, key, txt=None, d=0.44, n_sites=1, pitch=0.62, size=None):
    txt = key if txt is None else txt
    if size is None:
        size = (14 if len(txt) <= 2 else 13) if d >= 0.4 else (12 if len(txt) <= 1 else 10)
    if n_sites == 1:
        return circle(sl, cx, cy, d, fill=key, line=None, content=txt, size=size,
                      color=SPECIES_TEXT.get(key, 'white'), bold=True)
    w = pitch * (n_sites - 1) + d
    return shape(sl, 'rrect', cx - w / 2, cy - d / 2, w, d, fill=key, line=None, radius='pill',
                 content=txt, size=size, color=SPECIES_TEXT.get(key, 'white'), bold=True, margin=0.0)


def slab(sl, x, y, n, pitch=0.62, d=0.44, occ=None, metal_label='Cu(111)', label_size=12):
    """Side view: n sites in a row on a metal slab. y is the centre line of the sites.
    occ: {first site index: (species key, text, n_sites)}. Returns the site centre x list."""
    occ = occ or {}
    shape(sl, 'rrect', x, y + d / 2 + 0.03, n * pitch, 0.34, fill='ni', line=None, radius=0.06,
          content=metal_label if metal_label else '', size=label_size, color='muted')
    covered = set()
    for start, (_, _, ns) in occ.items():
        covered.update(range(start, start + ns))
    xs = [x + pitch * (i + 0.5) for i in range(n)]
    for i in range(n):
        if i not in covered:
            site(sl, xs[i], y, d)
    for start, (key, txt, ns) in occ.items():
        cx = xs[start] + pitch * (ns - 1) / 2
        adsorbate(sl, cx, y, key, txt, d, ns, pitch)
    return xs


def grid(sl, x, y, cols, rows, pitch=0.46, d=0.38, occ=None, plate=True):
    """Top view: a close-packed (hexagonal) layer of sites. occ: {(row, col): (key, text, n_sites)};
    a two-site species covers (row, col) and (row, col + 1). Returns {(r, c): (cx, cy)}."""
    occ = occ or {}
    pos = {}
    dy = pitch * math.sqrt(3) / 2
    for r in range(rows):
        for c in range(cols):
            pos[(r, c)] = (x + d / 2 + c * pitch + (r % 2) * pitch / 2, y + d / 2 + r * dy)
    if plate:
        pw = (cols - 1) * pitch + pitch / 2 + d + 0.24
        ph = (rows - 1) * dy + d + 0.24
        shape(sl, 'rrect', x - 0.12, y - 0.12, pw, ph, fill='panel', line=None, radius=0.1)
    covered = set()
    for (r, c), (_, _, ns) in occ.items():
        covered.update((r, c + k) for k in range(ns))
    for key in pos:
        if key not in covered:
            site(sl, *pos[key], d=d, lw=1.0)
    for (r, c), (key, txt, ns) in occ.items():
        cx, cy = pos[(r, c)]
        cx += pitch * (ns - 1) / 2
        adsorbate(sl, cx, cy, key, txt, d, ns, pitch)
    return pos


def h2(sl, cx, cy, d=0.36, key='H'):
    """Gas-phase H2: two touching atoms."""
    circle(sl, cx - d * 0.45, cy, d, fill=key, line=None, content='H', size=13, color='white', bold=True)
    circle(sl, cx + d * 0.45, cy, d, fill=key, line=None, content='H', size=13, color='white', bold=True)


def gas(sl, cx, cy, key, txt=None, w=None, h=0.36, size=13):
    """A gas molecule drawn as a labelled capsule in its species colour."""
    txt = key if txt is None else txt
    w = w or max(0.5, 0.16 + 0.11 * len(txt))
    return shape(sl, 'rrect', cx - w / 2, cy - h / 2, w, h, fill=key, line=None, radius='pill',
                 content=txt, size=size, color='white', bold=True, margin=0.0)


def axes(sl, x, y, w, h, xlabel, ylabel, color='muted', size=14, ylabel_w=1.6):
    """L-shaped axes with arrowheads; (x, y) is the top-left of the plotting box."""
    arrow(sl, x, y + h, x + w, y + h, color=color, lw=1.25)
    arrow(sl, x, y + h, x, y - 0.05, color=color, lw=1.25)
    text(sl, x, y + h + 0.04, w, 0.3, xlabel, size=size, color=color, align='c')
    text(sl, x - ylabel_w - 0.06, y + h / 2 - 0.3, ylabel_w, 0.6, ylabel, size=size, color=color, align='r',
         anchor='m')


def level(sl, x, y, w, txt, color='ink', lw=2.5, size=14, label_side='below'):
    """An energy level: a thick horizontal bar with a label."""
    line(sl, x, y, x + w, y, color=color, lw=lw)
    ly = y + 0.04 if label_side == 'below' else y - 0.34
    text(sl, x - 0.3, ly, w + 0.6, 0.3, txt, size=size, color=color, align='c')


# --- panels and chips --------------------------------------------------------------------
def card(sl, x, y, w, h, heading=None, fill='panel', line=None, size=16, heading_color='ink', lw=1.0,
         dash=False):
    shape(sl, 'rrect', x, y, w, h, fill=fill, line=line, lw=lw, radius=0.1, dash=dash)
    if heading:
        text(sl, x + 0.15, y + 0.1, w - 0.3, 0.36, heading, size=size, color=heading_color, bold=True)


def chip(sl, x, y, txt, key=None, w=None, h=0.36, size=14, outline='faint'):
    """A pill label. key = species colour (filled) or None (white with outline). x, y = top-left."""
    w = w or (0.3 + 0.1 * len(re.sub(r'_\{|\^\{|\}', '', txt)))
    if key:
        return shape(sl, 'rrect', x, y, w, h, fill=key, line=None, radius='pill', content=txt, size=size,
                     color='white', bold=True, margin=0.0)
    return shape(sl, 'rrect', x, y, w, h, fill='white', line=outline, radius='pill', content=txt, size=size,
                 color='ink', margin=0.0)


def filecard(sl, x, y, w, h, title, size=14, fill='white'):
    shape(sl, 'doc', x, y, w, h, fill=fill, line='faint', lw=1.0)
    text(sl, x + 0.12, y + 0.08, w - 0.24, 0.32, title, size=size, color='ink', font=CODE, bold=True)


def tick(sl, x, y, ok=True, d=0.32):
    """A round ✓ (green) or ✗ (dark) mark centred at (x, y)."""
    circle(sl, x, y, d, fill='ok' if ok else 'dark', line=None, content='✓' if ok else '✗', size=13,
           color='white', bold=True)
