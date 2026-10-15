"""Run report: one PDF per run, built from the run's output folder and settings.py.

Nothing is solved: the report reads the tables (.csv) and figures (.pdf/.png) the analyses
wrote, so it can be rebuilt at any time from a finished or partly finished run. Every
analysis with outputs gets a section: its goal, what was computed, key numbers, its tables
(all borders), its figures, and how to read each figure panel by panel. The fixed text is in
report_text.py.

    from hydrogen_surface_kinetics import report
    report.build('runs/try_4/Cu/settings.py')      # -> runs/try_4/Cu/output/report/try_4_Cu_report.pdf

LaTeX (latexmk + pdflatex) builds the PDF; without it the .tex is kept and the log says so.
"""
import datetime
import glob
import math
import os
import re
import shutil
import subprocess

import numpy as np
import pandas as pd

from . import __version__
from . import report_text as RT
from .config import ANALYSIS_ORDER, DEFAULTS, REPORT_DEFAULTS, load_config

TEXT_WIDTH_CM = 17.4                # A4 with 1.8 cm margins
LATEX_TIMEOUT = 900                 # s

# The outputs of each analysis (glob patterns in the output folder): a section is written
# for every analysis with at least one of them.
SECTION_OUTPUTS = {
    'steady_state': ['steady_state_results.csv', 'rates_summary_vs_*.p[dn][fg]'],
    'composition': ['steady_state_surface_species.csv', 'steady_state_coverages.csv', 'surface_composition_vs_*.p[dn][fg]'],
    'diagnostics': ['mechanism_binding_check.csv', 'mechanism_switch_off.csv', 'top_reactions.csv', 'binding_check.p[dn][fg]'],
    'convergence': ['convergence.csv', 'convergence*.p[dn][fg]'],
    'curves': ['poisoning_curves.csv', 'poisoning_tolerance.csv', 'poisoning_curves_vs_*.p[dn][fg]', 'site_blocking*.p[dn][fg]'],
    'sweep': ['temperature_sweep.csv', 'release_temperature.csv', 'temperature_sweep*.p[dn][fg]', 'sweep_coverages*.p[dn][fg]'],
    'arrhenius': ['arrhenius_fits.csv', 'arrhenius_intrinsic.csv', 'arrhenius_local.csv', 'arrhenius*.p[dn][fg]'],
    'pressure_sweep': ['pressure_sweep.csv', 'half_rate_pressure.csv', 'pressure_sweep*.p[dn][fg]'],
    'orders': ['reaction_orders.csv', 'reaction_orders_vs_*.p[dn][fg]'],
    'rate_control': ['sensitivity.csv', 'rate_control_*.p[dn][fg]'],
    'paths': ['reaction_paths/path_*.p[dn][fg]'],
    'transients': ['poisoning_over_time*.csv', 'poisoning_over_time_*.p[dn][fg]'],
    'flux': ['flux/*/flux_diagram*.png', 'flux/*/flux_diagram.avi'],
}
PATH_RE = re.compile(r'path_(?P<el>[A-Za-z]+)_(?P<case>.+)_(?P<T>[\d.]+)K_(?P<P>[0-9.e+-]+)bar_(?P<flow>[A-Za-z]+)\.(png|pdf)')


class ReportError(RuntimeError):
    pass


# --- LaTeX text ----------------------------------------------------------------------------
_SPECIAL = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$', '#': r'\#', '_': r'\_',
            '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}', '^': r'\textasciicircum{}',
            '<': r'\textless{}', '>': r'\textgreater{}', '|': r'\textbar{}'}
_UNICODE = {'−': '$-$', '×': r'$\times$', 'θ': r'$\theta$', 'Π': r'$\Pi$', '≈': r'$\approx$', '⇌': r'$\rightleftharpoons$',
            '→': r'$\rightarrow$', '←': r'$\leftarrow$', '°': r'$^\circ$', 'µ': r'$\mu$', 'μ': r'$\mu$', '≤': r'$\le$',
            '≥': r'$\ge$', '·': r'$\cdot$', '±': r'$\pm$', 'Δ': r'$\Delta$', '…': r'\ldots{}', '–': '--', '—': '---',
            '’': "'", '‘': '`', '“': '``', '”': "''", '²': '$^2$', '³': '$^3$', '⁻': '$^-$', '¹': '$^1$',
            '₀': '$_0$', '₁': '$_1$', '₂': '$_2$', '₃': '$_3$', '₄': '$_4$', 'α': r'$\alpha$', 'γ': r'$\gamma$',
            'σ': r'$\sigma$', 'ν': r'$\nu$', 'Γ': r'$\Gamma$', 'ε': r'$\varepsilon$', 'Å': r'\AA{}', '∗': r'$\ast$'}
_ASCII = {'−': '-', '×': 'x', 'θ': 'theta', 'Π': 'Pi', '≈': '~', '→': '->', '←': '<-', '°': 'deg', 'µ': 'u', 'μ': 'u',
          '≤': '<=', '≥': '>=', '·': '*', '±': '+/-', 'Δ': 'd', '…': '...', '–': '-', '—': '-', '’': "'", '‘': "'",
          '“': '"', '”': '"', '²': '2', '³': '3', '⁻': '-', '¹': '1', '₀': '0', '₁': '1', '₂': '2', '₃': '3', '₄': '4'}


def tex_escape(text):
    """Plain text (species names, labels, file names) -> LaTeX. Unknown non-ASCII -> '?'."""
    out = []
    for ch in str(text):
        if ch in _SPECIAL:
            out.append(_SPECIAL[ch])
        elif ch in _UNICODE:
            out.append(_UNICODE[ch])
        elif ord(ch) < 128:
            out.append(ch)
        else:
            out.append('?')
    return ''.join(out)


def tex_formula(text):
    """A species label as LaTeX text: digits after an element symbol become subscripts and
    '*' a site, e.g. 'H2' -> H\\textsubscript{2}, 'H*' -> H$\\ast$."""
    if text is None or (isinstance(text, float) and math.isnan(text)):
        return '--'
    s = tex_escape(text).replace('*', r'$\ast$')
    return re.sub(r'(?<=[A-Za-z)])(\d+)(?![^(]*\))', r'\\textsubscript{\1}', s)


def tex_num(v, sig=3):
    """A number in math mode, 3 significant figures; scientific outside [1e-3, 1e5).
    NaN/None -> '--', booleans -> yes/no, strings escaped."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return '--'
    if isinstance(v, (bool, np.bool_)):
        return 'yes' if v else 'no'
    if isinstance(v, str):
        return tex_escape(v)
    v = float(v)
    if not math.isfinite(v):
        return r'$\infty$' if v > 0 else r'$-\infty$'
    if v == 0:
        return '$0$'
    a = abs(v)
    if v.is_integer() and a < 1e5:
        return f'${int(v)}$'
    e = math.floor(math.log10(a))
    if -3 <= e < 5:
        if e >= sig - 1:
            s = f'{round(v, sig - 1 - e):.0f}'
        else:
            s = f'{v:#.{sig}g}'.rstrip('.')         # '#' keeps trailing zeros: 0.100, 1.00
            if 'e' in s:
                s = f'{v:.0f}'
        return f'${s}$'
    m = round(v / 10 ** e, sig - 1)
    if abs(m) >= 10:
        m, e = m / 10, e + 1
    return f'${m:.{sig - 1}f}{{\\times}}10^{{{e}}}$'


def tex_coverage(v):
    """A coverage: values below 1e-12 (numerical noise) print as 0."""
    return '$0$' if isinstance(v, (int, float)) and abs(v) < 1e-12 else tex_num(v)


def tex_kj(v):
    """An energy in kJ/mol with one decimal and a sign."""
    if v is None or not np.isfinite(v):
        return '--'
    return f'${v:+.1f}$'


def tex_flag(v):
    """True -> yes, False -> no, missing -> --."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return '--'
    if isinstance(v, str):
        return 'yes' if v.strip().lower() == 'true' else 'no' if v.strip().lower() == 'false' else tex_escape(v)
    return 'yes' if bool(v) else 'no'


def tex_arrhenius(A, Ea):
    """A exp(-Ea/RT) with Ea in kJ/mol, e.g. $4.21{\\times}10^{6}\\,\\exp(-55.1/RT)$."""
    if A is None or Ea is None or not (np.isfinite(A) and np.isfinite(Ea)):
        return '--'
    a = tex_num(A)[1:-1]
    if abs(Ea) < 0.05:
        return f'${a}$ (no $T$ dependence)'
    return f'${a}\\,\\exp({"-" if Ea > 0 else "+"}{abs(Ea):.1f}/RT)$'


def tex_path(path):
    """A file path that may break at '/' (url's \\path); escaped typewriter text if it has
    characters \\path cannot take."""
    path = str(path)
    if re.search(r'[{}%#\\]', path):
        return r'\texttt{' + tex_escape(path) + '}'
    return r'\path{' + path + '}'


def tex_file(name):
    return r'\texttt{' + tex_escape(name) + '}'


def plain_width(tex):
    """Rough printed width (characters) of a LaTeX snippet."""
    s = re.sub(r'\\(textsubscript|mathrm|textbf|emph|texttt|bfseries)\{([^}]*)\}', r'\2', tex)
    s = re.sub(r'\\(exp|ln|log)\b', r'\1', s)
    s = re.sub(r'\\[a-zA-Z]+', 'x', s)
    s = s.replace(r'\,', ' ')
    s = re.sub(r'[${}^_\\]', '', s)
    return len(s)


def is_missing(v):
    return v is None or (isinstance(v, float) and math.isnan(v))


def ascii_text(text):
    return ''.join(ch if ord(ch) < 128 else _ASCII.get(ch, '?') for ch in text)


def fig_refs(labels):
    """'Figure 3' / 'Figures 3--5, 8' from labels fig:<n>."""
    nums = sorted(int(lab.split(':')[1]) for lab in labels)
    runs = []
    for n in nums:
        if runs and n == runs[-1][1] + 1:
            runs[-1][1] = n
        else:
            runs.append([n, n])
    parts = [rf'\ref{{fig:{a}}}' if a == b else rf'\ref{{fig:{a}}}--\ref{{fig:{b}}}' for a, b in runs]
    return ('Figure ' if len(nums) == 1 else 'Figures ') + ', '.join(parts)


def bullets(items):
    items = [i for i in items if i]
    return (r'\begin{itemize}[leftmargin=*]' + ''.join(r'\item ' + i for i in items) + r'\end{itemize}') if items else ''


def point_tex(T, Pb):
    return f'{T:g} K, {Pb:g} bar'


def finite(df, col):
    return df[col].replace([np.inf, -np.inf], np.nan).dropna() if col in df else pd.Series(dtype=float)


# --- Report ---------------------------------------------------------------------------------
class Report:
    """The LaTeX of one run's report, built section by section."""

    def __init__(self, cfg, model=None, log=print):
        self.cfg, self.model, self.log = cfg, model, log
        self.out = cfg.output_dir
        self.dir = os.path.join(self.out, 'report')
        self.opts = {**REPORT_DEFAULTS, **(cfg.get('report') or {})}
        self.used = set()                   # output files shown or referenced in the report
        self.n_fig = 0

    # files ------------------------------------------------------------------------------
    def csv(self, name):
        path = os.path.join(self.out, name)
        if not os.path.exists(path):
            return None
        self.used.add(name)
        return pd.read_csv(path)

    def outputs(self, name):
        """The output files of one analysis (paths relative to the output folder)."""
        found = set()
        for pattern in SECTION_OUTPUTS.get(name, []):
            found.update(os.path.relpath(p, self.out) for p in glob.glob(os.path.join(self.out, pattern)))
        return sorted(found)

    def figure_files(self, stem):
        """The files of one saved figure, in order: stem.pdf (else .png), or its parts stem_1, stem_2, ..."""
        def best(base):
            for ext in ('.pdf', '.png'):
                if os.path.exists(os.path.join(self.out, base + ext)):
                    return base + ext
            return None
        one = best(stem)
        if one:
            return [one]
        parts = set()
        for path in glob.glob(os.path.join(self.out, glob.escape(stem) + '_*.*')):
            m = re.fullmatch(re.escape(os.path.basename(stem)) + r'_(\d+)\.(pdf|png)', os.path.basename(path))
            if m:
                parts.add(int(m.group(1)))
        return [f for f in (best(f'{stem}_{i}') for i in sorted(parts)) if f]

    def stems(self, pattern):
        """Stems of the saved figures (parts merged) whose name matches the regex `pattern`."""
        found = set()
        for path in glob.glob(os.path.join(self.out, '*.p[dn][fg]')):
            stem = os.path.splitext(os.path.basename(path))[0]
            m = re.fullmatch(r'(.+)_\d+', stem)
            if m and not any(os.path.exists(os.path.join(self.out, stem.rsplit('_', 1)[0] + e)) for e in ('.pdf', '.png')) \
                    and self.figure_files(m.group(1)):
                stem = m.group(1)
            if re.fullmatch(pattern, stem):
                found.add(stem)
        return sorted(found)

    # points and cases -------------------------------------------------------------------
    def case_order(self):
        from .conditions import case_labels
        return case_labels(self.cfg)

    def sort_cases(self, df, by=()):
        order = {c: i for i, c in enumerate(self.case_order())}
        return df.assign(_o=df['case'].map(order).fillna(len(order))).sort_values(['_o', *by]).drop(columns='_o')

    def point_order(self):
        """(T, P bar) of the series points: the temperature series first, then each pressure series."""
        from .conditions import series_points
        return [(p['T'], p['P'] / 1e5) for p in series_points(self.cfg)]

    def sort_points(self, pts):
        order = {pt: i for i, pt in enumerate(self.point_order())}
        return sorted(set(pts), key=lambda p: (order.get(p, len(order)), p))

    def probe_tex(self):
        if self.model is not None:
            from .labels import probe_label, species_label
            return tex_formula(probe_label(self.model)), tex_formula(species_label(self.model, self.model.probe_adsorbate))
        return tex_formula(self.cfg.probe_label or self.cfg.probe_gas), tex_formula(self.cfg.probe_adsorbate)

    # LaTeX blocks ---------------------------------------------------------------------------
    def table(self, df, columns, caption, source=None, group_by=None, full=False):
        """A bordered longtable. columns: [(column, header LaTeX, align, formatter or None)];
        align 'l', 'r', or 'P<cm>' for a wrapping column of that width. Rows beyond max_rows
        are not shown (unless full) and the caption names the file with all of them.
        group_by: a column; a heavier rule separates the rows where its value changes."""
        if df is None or not len(df):
            return ''
        n_all = len(df)
        if n_all > self.opts['max_rows'] and not full:
            df = df.head(self.opts['max_rows'])
        recs = df.to_dict('records')
        cells = [['--' if is_missing(rec.get(col)) else (fmt or tex_num)(rec.get(col)) for col, _, _, fmt in columns]
                 for rec in recs]
        heads = [h for _, h, _, _ in columns]

        fixed = {}                                  # printed width of the columns that do not wrap, per size
        for size, (cm_per_char, pad) in ((r'\scriptsize', (0.14, 0.23)), (r'\tiny', (0.115, 0.23))):
            fixed[size] = sum(pad + (0 if a.startswith('P') else cm_per_char * max(
                [max(plain_width(x) for x in h.split(r'\\'))] + [plain_width(c[i]) for c in cells]))
                for i, (_, h, a, _) in enumerate(columns))
        wraps = [float(a[1:]) for _, _, a, _ in columns if a.startswith('P')]
        # the first that fits: small font; small font with narrower wrapping columns; tiny; tiny narrower; landscape
        size, scale, landscape = r'\tiny', 1.0, True
        for sz, least in ((r'\scriptsize', None), (r'\scriptsize', 2.6), (r'\tiny', None), (r'\tiny', 2.0)):
            room = TEXT_WIDTH_CM - 0.4 - fixed[sz]              # 0.4 cm to spare: the widths are estimates
            if least is None and room >= sum(wraps):
                size, scale, landscape = sz, 1.0, False
                break
            if least is not None and wraps and room >= least * len(wraps) and room < sum(wraps):
                size, scale, landscape = sz, room / sum(wraps), False
                break
        spec = '|' + '|'.join(r'>{\raggedright\arraybackslash}p{' + f'{float(a[1:]) * scale:.2f}' + 'cm}'
                              if a.startswith('P') else a for _, _, a, _ in columns) + '|'
        head = ' & '.join(r'\makecell{\bfseries ' + h + '}' for h in heads) + r' \\ \thickhline'
        note = ''
        if source:
            self.used.add(source)
            note = (f' First {len(df)} of {n_all} rows; all of them in {tex_file(source)}.'
                    if len(df) < n_all else f' File: {tex_file(source)}.')
        lines = [r'\begingroup' + size + r'\setlength{\tabcolsep}{3pt}',
                 r'\begin{longtable}{' + spec + '}',
                 r'\caption{' + caption + note + r'}\\', r'\hline', head, r'\endfirsthead',
                 r'\hline', head, r'\endhead',
                 r'\multicolumn{' + str(len(columns)) + r'}{r}{\emph{continued on the next page}}\\', r'\endfoot',
                 r'\endlastfoot']
        for k, (rec, row) in enumerate(zip(recs, cells)):
            new_group = (group_by is not None and k + 1 < len(recs)
                         and recs[k + 1].get(group_by) != rec.get(group_by))
            lines.append(' & '.join(row) + (r' \\ \thickhline' if new_group else r' \\ \hline'))
        lines += [r'\end{longtable}', r'\endgroup']
        body = '\n'.join(lines)
        return (r'\begin{landscape}' + '\n' + body + '\n' + r'\end{landscape}') if landscape else body

    def series_tables(self, df, value_cols, what, source, lead=(), sort_extra=(), group_by='case'):
        """A results table split like the series: the temperature series (rows by case, then T)
        and one table per pressure-series temperature (rows by case, then p)."""
        if df is None or not len(df) or 'series' not in df:
            return ''
        case_col = [('case', 'case', 'l', tex_escape)]
        out = []
        P_ref = self.cfg.reference_pressure / 1e5
        t_part = df[df['series'].isin(['temperature', 'both']) & np.isclose(df['P_ref_bar'].astype(float), P_ref)]
        if len(t_part):
            out.append(self.table(self.sort_cases(t_part, (*sort_extra, 'T_K')),
                                  case_col + list(lead) + [('T_K', '$T$\\\\(K)', 'r', None)] + list(value_cols),
                                  f'{what}: temperature series at $p$(\\probe) = {P_ref:g} bar.', source, group_by))
        p_part = df[df['series'].isin(['pressure', 'both'])]
        for T in sorted(set(p_part['T_K'])):
            part = p_part[p_part['T_K'] == T]
            if part['P_ref_bar'].nunique() < 2:
                continue
            out.append(self.table(self.sort_cases(part, (*sort_extra, 'P_ref_bar')),
                                  case_col + list(lead) + [('P_ref_bar', '$p$(\\probe)\\\\(bar)', 'r', None)] + list(value_cols),
                                  f'{what}: pressure series at {T:g} K.', source, group_by))
        return '\n\n'.join(out)

    def pivot(self, df, index, columns, value, caption, index_head, source=None, fmt=None, index_fmt=None,
              index_order=None, col_order=None, col_head=None, flag=None):
        """One value as a grid: a row per `index` value, a column per `columns` value. flag: a
        boolean column; flagged cells get a dagger."""
        if df is None or not len(df) or value not in df:
            return ''
        piv = df.pivot_table(index=index, columns=columns, values=value, aggfunc='first', dropna=False)
        flags = None
        if flag and flag in df:
            flags = df.assign(_f=df[flag].fillna(False).astype(bool)).pivot_table(
                index=index, columns=columns, values='_f', aggfunc='max')
        rows = [r for r in (index_order or list(piv.index)) if r in piv.index]
        cols = [c for c in (col_order or list(piv.columns)) if c in piv.columns]
        recs = []
        for r in rows:
            rec = {'_idx': r}
            for i, c in enumerate(cols):
                s = (fmt or tex_num)(piv.loc[r, c])
                if flags is not None and c in flags.columns and r in flags.index and bool(flags.loc[r, c] or False):
                    s += r'$^\dagger$'
                rec[f'c{i}'] = s
            recs.append(rec)
        spec = [('_idx', index_head, 'l', index_fmt or tex_num)] + \
               [(f'c{i}', (col_head or tex_escape)(c), 'r', lambda s: s) for i, c in enumerate(cols)]
        return self.table(pd.DataFrame(recs), spec, caption, source)

    def figure(self, stem, caption, height=0.72):
        """The figure (all its parts) as floats placed here. Returns (LaTeX, figure labels)."""
        files = self.figure_files(stem)
        out, labels = [], []
        for k, f in enumerate(files):
            self.used.add(f)
            self.n_fig += 1
            label = f'fig:{self.n_fig}'
            labels.append(label)
            part = f' (part {k + 1} of {len(files)})' if len(files) > 1 else ''
            out.append(r'\begin{figure}[H]\centering'
                       rf'\includegraphics[width=\linewidth,height={height}\textheight,keepaspectratio]{{../{f}}}'
                       r'\caption{' + caption + part + r' File: ' + tex_file(f) + r'.}\label{' + label + '}'
                       r'\end{figure}')
        return '\n'.join(out), labels

    @staticmethod
    def how_to_read(spec, labels):
        """Per panel: what it shows, how to read it, and what each trend means."""
        out = [r'\subsubsection*{' + fig_refs(labels) + ': ' + spec['caption'] + '}',
               r'\begin{description}[leftmargin=1.2em,style=nextline]']
        for name, shows, read, trends in spec['panels']:
            out.append(r'\item[' + name + '] ' + shows + ' ' + read)
            if trends:
                out.append(r'\begin{itemize}[leftmargin=1.4em,itemsep=0pt,topsep=2pt]')
                out += [r'\item \emph{' + see + r'} $\Rightarrow$ ' + means + '.' for see, means in trends]
                out.append(r'\end{itemize}')
        out.append(r'\end{description}')
        return '\n'.join(out)

    @staticmethod
    def column_list(columns):
        return (r'\begin{description}[leftmargin=1.2em,itemsep=0pt,font=\normalfont\itshape]' + '\n'
                + '\n'.join(r'\item[' + c + '] ' + text for c, text in columns) + '\n' + r'\end{description}')

    def section(self, name, parts):
        """One analysis: goal, then parts, then how to read the tables and each kind of figure
        shown. parts: (subsection title or None, LaTeX) or ('figure', key, stem[, caption
        addition[, height]]), key naming the figure's text in report_text."""
        text = RT.SECTIONS[name]
        specs = {f['key']: f for f in text.get('figures', [])}
        body = [r'\section{' + text['title'] + '}', r'\label{sec:' + name + '}',
                r'\paragraph{Goal.} ' + text['goal']]
        shown = {}
        for part in parts:
            if part[0] == 'figure':
                key, stem = part[1], part[2]
                extra = part[3] if len(part) > 3 else ''
                height = part[4] if len(part) > 4 else 0.72
                caption = specs[key]['caption'] + (' ' + extra if extra else '')
                tex, labels = self.figure(stem, caption, height)
                if labels:
                    body.append(tex)
                    shown.setdefault(key, []).extend(labels)
                continue
            title, content = part
            if not content:
                continue
            if title:
                body.append(r'\subsection*{' + title + '}')
            body.append(content)
        if text.get('columns'):
            body += [r'\subsection*{How to read the tables}', self.column_list(text['columns'])]
        if shown:
            body.append(r'\subsection*{How to read the figures}')
            body += [self.how_to_read(specs[key], labels) for key, labels in shown.items()]
        return '\n\n'.join(body)

    def series_figures(self, key, base):
        """('figure', ...) parts for <base>_vs_T and <base>_vs_P."""
        return [('figure', key, f'{base}_vs_T', '(Temperature series.)'),
                ('figure', key, f'{base}_vs_P', '(Pressure series.)')]


# --- Run facts ---------------------------------------------------------------------------
def run_name(cfg):
    """The run folder below runs/ (e.g. 'try_4/Cu'), else the folder's name."""
    parts = os.path.normpath(cfg.run_dir).split(os.sep)
    if 'runs' in parts:
        i = len(parts) - 1 - parts[::-1].index('runs')
        if i < len(parts) - 1:
            return '/'.join(parts[i + 1:])
    return parts[-1]


def rmg_metal(cfg):
    """The catalyst of the RMG input.py next to settings.py, e.g. 'Cu111' (None if absent)."""
    path = os.path.join(cfg.run_dir, 'input.py')
    if not os.path.exists(path):
        return None
    with open(path, errors='replace') as f:
        m = re.search(r'catalystProperties\(\s*metal\s*=\s*[\'"]([^\'"]+)[\'"]', f.read())
    return m.group(1) if m else None


def package_commit():
    here = os.path.dirname(os.path.abspath(__file__))
    try:
        sha = subprocess.run(['git', '-C', here, 'rev-parse', '--short', 'HEAD'], capture_output=True,
                             text=True, timeout=10).stdout.strip()
        dirty = subprocess.run(['git', '-C', here, 'status', '--porcelain', '--', '.'], capture_output=True,
                               text=True, timeout=10).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return None
    return (sha + (' (with uncommitted changes)' if dirty else '')) if sha else None


def cantera_version():
    try:
        import cantera as ct
        return ct.__version__
    except ImportError:
        return None


def non_default_settings(cfg):
    """[(setting, value, default)] of the solver and run sections that differ from the defaults."""
    rows = []
    for section in ('solver', 'run'):
        for key, default in DEFAULTS[section].items():
            value = cfg[section].get(key, default)
            if value != default:
                rows.append((f'{section}.{key}', value, default))
    return rows


def status_word(rec):
    if not bool(rec.get('converged')):
        return 'not converged'
    agree = rec.get('starts_agree')
    if agree is None or (isinstance(agree, float) and math.isnan(agree)):
        return 'one start'
    return 'ok' if bool(agree) else 'two states'


def with_status(df):
    return df.assign(status=[status_word(r) for r in df.to_dict('records')]) if df is not None and len(df) else df


def value_or_text(v):
    return tex_num(v) if isinstance(v, (int, float)) and not isinstance(v, bool) else tex_escape(v)


# --- Pages -------------------------------------------------------------------------------
def preamble(rep, title):
    probe, ads = rep.probe_tex()
    return '\n'.join([
        r'\documentclass[10pt,a4paper]{article}',
        r'\usepackage[T1]{fontenc}', r'\usepackage[utf8]{inputenc}', r'\usepackage{lmodern}',
        r'\usepackage[margin=1.8cm]{geometry}', r'\usepackage{amsmath,amssymb}', r'\usepackage{graphicx}',
        r'\usepackage{longtable,array,makecell}', r'\usepackage{xcolor}', r'\usepackage{float}',
        r'\usepackage{caption}', r'\usepackage{enumitem}', r'\usepackage{pdflscape}', r'\usepackage{listings}',
        r'\usepackage[hidelinks]{hyperref}',
        r'\newcommand{\thickhline}{\noalign{\global\arrayrulewidth=0.9pt}\hline\noalign{\global\arrayrulewidth=0.4pt}}',
        r'\renewcommand{\arraystretch}{1.18}', r'\setlength{\parindent}{0pt}', r'\setlength{\parskip}{4pt}',
        r'\setlength{\LTpre}{6pt}', r'\setlength{\LTpost}{8pt}',
        r'\captionsetup{font=small,labelfont=bf}', r'\captionsetup[longtable]{font=small,labelfont=bf}',
        r'\lstset{basicstyle=\ttfamily\scriptsize,breaklines=true,columns=fullflexible,frame=single,'
        r'commentstyle=\color{gray},keywordstyle=\bfseries,language=Python}',
        r'\newcommand{\probe}{' + probe + '}', r'\newcommand{\adsorbate}{' + ads + '}',
        r'\hypersetup{pdftitle={' + title + '}}',
    ])


def title_page(rep):
    cfg, model = rep.cfg, rep.model
    name, metal = run_name(cfg), rmg_metal(cfg)
    probe, _ = rep.probe_tex()
    title = (f'{probe} on {tex_escape(metal)}' if metal else probe) + ': surface kinetics report'
    facts = [('run', r'\texttt{' + tex_escape(name) + '}'), ('settings', tex_path(cfg.settings_file)),
             ('output folder', tex_path(cfg.output_dir)), ('mechanism', tex_path(cfg.mechanism))]
    if model is not None:
        facts.append(('mechanism sha256', r'\texttt{' + model.mechanism_sha256[:16] + r'\ldots}'))
        if model.mechanism_file != cfg.mechanism:
            facts.append(('loaded from', tex_path(model.mechanism_file) + ' (duplicate flags repaired)'))
    if metal:
        facts.append(('catalyst (RMG input.py)', tex_escape(metal)))
    commit = package_commit()
    facts += [('report built', datetime.datetime.now().strftime('%Y-%m-%d %H:%M')),
              ('package', f'hydrogen\\_surface\\_kinetics {__version__}' + (f', commit {tex_escape(commit)}' if commit else '')),
              ('Cantera', tex_escape(cantera_version() or 'not available'))]
    rows = '\n'.join(r'\textbf{' + k + '} & ' + v + r' \\ \hline' for k, v in facts)
    body = [r'\begin{titlepage}\centering', r'\vspace*{3cm}',
            r'{\LARGE\bfseries ' + title + r'\par}', r'\vspace{0.6cm}',
            r'{\large run \texttt{' + tex_escape(name) + r'}\par}', r'\vspace{1.5cm}',
            r'\small\renewcommand{\arraystretch}{1.4}',
            r'\begin{tabular}{|>{\raggedright\arraybackslash}p{4cm}|>{\raggedright\arraybackslash}p{11.5cm}|}\hline',
            rows, r'\end{tabular}', r'\vfill',
            r'{\small Built from the run\textquotesingle s output files; nothing was solved to make it. '
            r'Rebuild with \texttt{bash run.sh runs/' + tex_escape(name) + r' report}.\par}',
            r'\end{titlepage}']
    return title, '\n'.join(body)


def run_summary(rep, built):
    from .conditions import ratio_for, series_points
    cfg, model = rep.cfg, rep.model
    out = [r'\section{Run summary}', r'\label{sec:summary}']

    # the system
    if model is not None:
        rxn = model.surf.reaction(model.i_probe)
        sys_rows = [('gas phase', f'{model.gas.n_species} species, {model.gas.n_reactions} reactions'),
                    ('surface phase', f'{model.surf.n_species} species, {model.surf.n_reactions} reactions'),
                    ('site density', tex_num(model.gamma * 1e-1) + r' mol/cm$^2$ (' + tex_num(model.gamma) + r' kmol/m$^2$)'),
                    ('probe gas', tex_file(cfg.probe_gas)), ('empty site', tex_file(cfg.empty_site)),
                    ('probe adsorbate', tex_file(cfg.probe_adsorbate)),
                    ('molecular adsorbate', tex_file(cfg.molecular_adsorbate) if cfg.molecular_adsorbate else 'none'),
                    ('probe reaction', f'\\#{model.i_probe}: ' + tex_escape(rxn.equation)
                     + (' (written as desorption; read in the adsorption direction)' if model.sign_probe == -1 else ''))]
        out += [r'\subsection{The system}', r'\begin{tabular}{|l|>{\raggedright\arraybackslash}p{12cm}|}\hline',
                '\n'.join(r'\textbf{' + k + '} & ' + v + r' \\ \hline' for k, v in sys_rows), r'\end{tabular}']

    # conditions
    pts = series_points(cfg)
    t_pts = sorted(p['T'] for p in pts if p['series'] in ('temperature', 'both'))
    P_ref = cfg.reference_pressure
    cond = [r'\item \textbf{Temperature series:} ' + ', '.join(f'{T:g}' for T in t_pts)
            + r' K at $p$(\probe) = ' + f'{P_ref / 1e5:g}' + ' bar.']
    p_pts = [p for p in pts if p['series'] in ('pressure', 'both')]
    for T in sorted({p['T'] for p in p_pts}):
        Ps = ', '.join(f"{p['P'] / 1e5:g}" for p in sorted(p_pts, key=lambda p: p['P']) if p['T'] == T)
        cond.append(r'\textbf{Pressure series at ' + f'{T:g}' + r' K:} $p$(\probe) = ' + Ps + ' bar.')
    if not p_pts:
        cond.append(r'\textbf{Pressure series:} none.')
    n_cases = len(rep.case_order())
    cond.append(f'{len(pts)} distinct points $\\times$ {n_cases} cases = {len(pts) * n_cases} main conditions'
                + (' on each of two bases (fixed probe pressure, and fixed mole fractions).' if cfg.include_input_basis
                   else ' (impurities added at a fixed probe-gas pressure).'))
    if cfg.mixture_label and len(cfg.impurities) > 1:
        cond.append(r'Mixture case \emph{' + tex_escape(cfg.mixture_label) + '}: all impurities together, '
                    + ('each at its own ratio.' if cfg.mixture_mode == 'same_each' else 'their ratios divided by the number of impurities.'))
    out += [r'\subsection{Conditions}', bullets([c.replace(r'\item ', '') for c in cond])]
    if cfg.impurities:
        rows = []
        for label in cfg.impurities:
            r = ratio_for(cfg.impurity_ratio, label)
            species = model.impurities[label] if model is not None else str(cfg.impurities[label])
            rows.append({'label': label, 'species': species, 'r': r, 'ppm': r * 1e6, 'p_bar': r * P_ref / 1e5})
        out.append(rep.table(pd.DataFrame(rows), [
            ('label', 'case', 'l', tex_escape), ('species', 'mechanism\\\\species', 'l', tex_file),
            ('r', r'ratio $r$\\$p_\mathrm{imp}/p$(\probe)', 'r', None), ('ppm', r'ppm of\\the probe gas', 'r', None),
            ('p_bar', r'$p_\mathrm{imp}$ (bar) at\\' + f'{P_ref / 1e5:g}' + r' bar \probe', 'r', None)],
            r'Impurities and their levels. The same ratio applies at every point of both series.'))

    # analyses
    rows = [{'analysis': name, 'in settings': name in cfg.analyses, 'outputs': bool(rep.outputs(name)),
             'section': (r'\ref{sec:' + name + '}') if name in built else '--'} for name in ANALYSIS_ORDER]
    out += [r'\subsection{Analyses}',
            rep.table(pd.DataFrame(rows), [('analysis', 'analysis', 'l', tex_file),
                                           ('in settings', r'in\\settings.py', 'l', None), ('outputs', r'outputs\\found', 'l', None),
                                           ('section', r'section', 'l', lambda s: s)],
                      'The analyses of the package, in the order they run. A section is shown for each one with outputs.')]

    # solver settings
    nd = non_default_settings(cfg)
    if nd:
        out += [r'\subsection{Solver and run settings that differ from the defaults}',
                rep.table(pd.DataFrame(nd, columns=['setting', 'value', 'default']),
                          [('setting', 'setting', 'l', tex_file), ('value', 'value', 'r', value_or_text),
                           ('default', 'default', 'r', value_or_text)],
                          'Settings changed from the package defaults.')]

    # warnings
    warns = []
    if model is not None:
        warns += [tex_escape(w) for w in model.warnings]
        if model.sticking_warnings:
            warns.append(f'{len(model.sticking_warnings)} sticking coefficients exceed 1 at some temperature '
                         '(Cantera warns when the mechanism is loaded); first: ' + tex_escape(model.sticking_warnings[0][:300]))
    else:
        warns.append('the mechanism was not loaded for this report, so its load warnings are not listed')
    out += [r'\subsection{Warnings from the mechanism}', bullets(warns) or 'None.']

    # convergence
    rows = []
    for path in sorted(glob.glob(os.path.join(rep.out, '*.csv'))):
        f = os.path.basename(path)
        df = pd.read_csv(path)
        if 'converged' not in df:
            continue
        conv = df['converged'].astype(bool)
        two = int((df['starts_agree'] == False).sum()) if 'starts_agree' in df else np.nan    # noqa: E712
        rows.append({'file': f, 'rows': len(df), 'converged': int(conv.sum()), 'not': int((~conv).sum()), 'two': two})
    if rows:
        out += [r'\subsection{Convergence}',
                rep.table(pd.DataFrame(rows), [('file', 'table', 'l', tex_file), ('rows', 'rows', 'r', None),
                                               ('converged', 'converged', 'r', None),
                                               ('not', r'not\\converged', 'r', None), ('two', r'two steady\\states', 'r', None)],
                          'Steady states in each results table. Two steady states: the clean and the covered start '
                          'ended in different states (history dependence).')]
        ss = rep.csv('steady_state_results.csv')
        if ss is not None:
            ss = with_status(ss)
            bad = ss[ss['status'] != 'ok']
            if len(bad):
                out.append(rep.table(rep.sort_cases(bad, ('T_K', 'P_ref_bar')), [
                    ('case', 'case', 'l', tex_escape), ('T_K', '$T$ (K)', 'r', None), ('P_ref_bar', r'$p$(\probe) (bar)', 'r', None),
                    ('series', 'series', 'l', tex_escape), ('status', 'status', 'l', tex_escape),
                    ('start_A_solved_by', r'start A\\solved by', 'l', tex_escape), ('start_B_solved_by', r'start B\\solved by', 'l', tex_escape)],
                    'Main conditions that are not plain single steady states.'))
            else:
                out.append(r'All main conditions converged to a single steady state from both starts.')
    out += [r'\subsection{Conventions used in this report}', RT.CONVENTIONS]
    return '\n\n'.join(out)


# --- Sections ----------------------------------------------------------------------------
STEADY_COLS = [('P_impurity_bar', r'$p_\mathrm{imp}$\\(bar)', 'r', None),
               ('tof_fwd', r'diss. TOF\\(s$^{-1}$)', 'r', None), ('tof_rev', r'rec. TOF\\(s$^{-1}$)', 'r', None),
               ('tof_net', r'net TOF\\(s$^{-1}$)', 'r', None), ('theta_free', r'$\theta_\ast$', 'r', tex_coverage),
               ('theta_probe', r'$\theta_\mathrm{H}$', 'r', tex_coverage)]


def section_steady_state(rep):
    cfg = rep.cfg
    ss = with_status(rep.csv('steady_state_results.csv'))
    parts = []
    if ss is not None:
        main = ss[ss['basis'] == 'fixed_probe']
        ref = cfg.reference_label
        n_pts = len(main[['T_K', 'P_ref_bar']].drop_duplicates())
        parts.append(('What was computed',
                      f'{len(main)} steady states: {n_pts} points (temperature and pressure series) $\\times$ '
                      f'{main["case"].nunique()} cases ({tex_escape(", ".join(rep.case_order()))}). '
                      'Each was solved from a clean surface and from a covered one; the rates, coverages and '
                      'poisoning factors are those of the clean start.'))
        keys = []
        counts = main['status'].value_counts()
        keys.append(f"{counts.get('ok', 0)} of {len(main)} conditions are single, converged steady states"
                    + ''.join(f"; {n} {tex_escape(s)}" for s, n in counts.items() if s != 'ok') + '.')
        ok = main[main['converged'].astype(bool)].assign(two=lambda d: d['status'] == 'two states')
        r_ok = ok[ok['case'] == ref]
        if len(r_ok):
            lo, hi = r_ok['tof_fwd'].idxmin(), r_ok['tof_fwd'].idxmax()
            keys.append(f"Reference ({tex_escape(ref)}) dissociation TOF from {tex_num(r_ok.loc[lo, 'tof_fwd'])} s$^{{-1}}$ "
                        f"({point_tex(r_ok.loc[lo, 'T_K'], r_ok.loc[lo, 'P_ref_bar'])}) to {tex_num(r_ok.loc[hi, 'tof_fwd'])} s$^{{-1}}$ "
                        f"({point_tex(r_ok.loc[hi, 'T_K'], r_ok.loc[hi, 'P_ref_bar'])}); empty sites "
                        f"$\\theta_\\ast$ = {tex_num(r_ok['theta_free'].min())} to {tex_num(r_ok['theta_free'].max())}.")
            eq = (ok['tof_net'].abs() / ok['tof_fwd']).replace([np.inf], np.nan)
            if eq.notna().any():
                keys.append(r'Distance from equilibrium of the probe step, $|$net$|$/dissociation: at most '
                            + tex_num(eq.max()) + (r' --- the step is at equilibrium everywhere (dissociation $\approx$ recombination).'
                                                   if eq.max() < 1e-3 else '.'))
        worst = ok[ok['case'].isin(list(cfg.impurities))]
        if len(worst) and worst['slowdown_x'].notna().any():
            top = worst.loc[worst.groupby(['T_K', 'P_ref_bar'])['slowdown_x'].idxmax().dropna()]
            for case, grp in top.groupby('case', sort=False):
                pts = ', '.join(f"{point_tex(r['T_K'], r['P_ref_bar'])} ({tex_num(r['slowdown_x'])}$\\times$)"
                                for r in grp.sort_values(['T_K', 'P_ref_bar']).to_dict('records'))
                keys.append(f'Strongest single impurity: \\textbf{{{tex_escape(case)}}} at {pts}.')
        mix = ok[ok['case'] == cfg.mixture_label] if cfg.mixture_label else ok.iloc[0:0]
        if len(mix):
            keys.append(f"All impurities together ({tex_escape(cfg.mixture_label)}): slowdown "
                        f"{tex_num(mix['slowdown_x'].min())} to {tex_num(mix['slowdown_x'].max())}$\\times$.")
        parts.append(('Key numbers', bullets(keys)))
        if len(ok):
            pts = rep.sort_points(map(tuple, ok[['T_K', 'P_ref_bar']].to_numpy()))
            parts.append(('Slowdown at a glance', rep.pivot(
                ok, 'case', ['T_K', 'P_ref_bar'], 'slowdown_x',
                r'Slowdown $1/\Pi$ of dissociation (times slower than the probe gas alone) per case and point. 1 = no effect; $\dagger$ = two steady states.',
                'case', index_fmt=tex_escape, index_order=rep.case_order(), col_order=pts,
                col_head=lambda c: f'{c[0]:g} K\\\\{c[1]:g} bar', flag='two')))
        cols = list(STEADY_COLS)
        if cfg.molecular_adsorbate:
            cols.append(('theta_molec', r'$\theta_\mathrm{H_2}$', 'r', tex_coverage))
        cols += [('theta_impurity', r'$\theta_\mathrm{imp}$', 'r', tex_coverage), ('site_count_factor', r'site\\factor', 'r', None),
                 ('poisoning_factor', r'$\Pi$', 'r', None), ('slowdown_x', r'$1/\Pi$', 'r', None)]
        if 'poisoning_factor_B' in main:
            cols.append(('poisoning_factor_B', r'$\Pi_\mathrm{B}$', 'r', None))
        cols.append(('status', 'status', 'l', tex_escape))
        tables = [rep.series_tables(bdf, cols, 'Main steady states' + ('' if basis == 'fixed_probe' else ' (fixed mole-fraction basis)'),
                                    'steady_state_results.csv') for basis, bdf in ss.groupby('basis', sort=False)]
        parts.append(('All main steady states', '\n\n'.join(tables)))
    parts += [('figure', 'rates_summary_vs_T', 'rates_summary_vs_T'), ('figure', 'rates_summary_vs_P', 'rates_summary_vs_P')]
    return rep.section('steady_state', parts)


def section_composition(rep):
    sp = rep.csv('steady_state_surface_species.csv')
    rep.csv('steady_state_coverages.csv')
    parts = []
    if sp is not None and len(sp):
        sp = sp[sp['basis'] == 'fixed_probe'] if 'basis' in sp else sp
        n_cond = len(sp[['case', 'T_K', 'P_ref_bar']].drop_duplicates())
        parts.append(('What was computed',
                      f'The coverage of every surface species at the {n_cond} converged main conditions; species above '
                      f'{tex_num(rep.cfg.composition["min_listed"])} are listed in {tex_file("steady_state_surface_species.csv")}, '
                      f'all of them in {tex_file("steady_state_coverages.csv")}.'))
        ads = sp[sp['species'] != 'empty site'].sort_values('coverage', ascending=False)
        top = ads.groupby(['case', 'T_K', 'P_ref_bar', 'series'], sort=False).head(2)
        cells = (top.assign(txt=[f"{tex_formula(s)} {tex_num(c)}" for s, c in zip(top['species'], top['coverage'])])
                 .groupby(['case', 'T_K', 'P_ref_bar', 'series'], sort=False)['txt'].agg('; '.join).reset_index())
        keys = []
        first = ads.groupby(['case', 'T_K', 'P_ref_bar'], sort=False).head(1)
        for case in rep.case_order():
            d = first[first['case'] == case]
            if not len(d):
                continue
            main = d['species'].value_counts()
            s = main.index[0]
            dd = d[d['species'] == s]
            keys.append(f"\\textbf{{{tex_escape(case)}}}: {tex_formula(s)} holds the most sites at {len(dd)} of {len(d)} points "
                        f"($\\theta$ = {tex_num(dd['coverage'].min())} to {tex_num(dd['coverage'].max())})"
                        + (', otherwise ' + ', '.join(tex_formula(x) for x in main.index[1:]) if len(main) > 1 else '') + '.')
        parts.append(('Key numbers', bullets(keys)))
        parts.append(('The two most abundant adsorbates', rep.series_tables(
            cells, [('txt', 'most abundant adsorbates and their coverage', 'l', lambda s: s)],
            'The two most abundant adsorbates (empty sites not counted)', 'steady_state_surface_species.csv')))
        if len(sp) <= rep.opts['max_rows']:
            parts.append(('All listed species', rep.table(
                rep.sort_cases(sp, ('T_K', 'P_ref_bar', 'coverage')),
                [('case', 'case', 'l', tex_escape), ('T_K', '$T$ (K)', 'r', None), ('P_ref_bar', '$p$ (bar)', 'r', None),
                 ('species', 'species', 'l', tex_formula), ('mechanism_name', 'mechanism name', 'l', tex_file),
                 ('coverage', 'coverage', 'r', tex_coverage)],
                'Surface species above the listing threshold.', 'steady_state_surface_species.csv', group_by='case')))
    parts += rep.series_figures('surface_composition', 'surface_composition')
    return rep.section('composition', parts)


def section_diagnostics(rep):
    parts, keys = [], []
    bc = rep.csv('mechanism_binding_check.csv')
    if bc is not None and len(bc):
        dg = next((c for c in bc.columns if c.startswith('dG_')), None)
        verdicts = bc['verdict'].value_counts()
        keys.append('Binding check: ' + ', '.join(f'{n} {tex_escape(v)}' for v, n in verdicts.items()) + '.'
                    + ''.join(f" {tex_escape(r['adsorption'])}: {tex_escape(r['verdict'])} ({tex_kj(r['dH_298_mech'])} kJ/mol "
                              f"against {tex_escape(r['dH_literature'])})." for r in bc.to_dict('records') if r['verdict'] != 'OK'))
        cols = [('adsorption', 'adsorption', 'l', tex_escape), ('dH_298_mech', r'$\Delta H_{298}$\\mechanism', 'r', tex_kj),
                ('dH_literature', r'literature\\(kJ/mol)', 'l', tex_escape), ('verdict', 'verdict', 'l', tex_escape)]
        if dg:
            cols.append((dg, r'$\Delta G$\\at ' + dg[3:] + ' K', 'r', tex_kj))
        cols.append(('literature_note', 'source', 'P4.5', tex_escape))
        parts.append(('Binding energies against literature', rep.table(bc, cols, 'Adsorption enthalpies (kJ/mol): mechanism against literature.',
                                                                       'mechanism_binding_check.csv')))
    parts.append(('figure', 'binding_check', 'binding_check'))
    sw = rep.csv('mechanism_switch_off.csv')
    if sw is not None and len(sw):
        species = ', '.join(sorted(set(sw['switched_off'].astype(str))))
        ok = sw[sw['converged'].astype(bool)]
        if len(ok):
            i = (np.log(ok['tof_ratio'].abs())).abs().idxmax()
            r = ok.loc[i]
            keys.append(f'Switching off {tex_formula(species)} changes the dissociation rate by up to a factor '
                        f"{tex_num(r['tof_ratio'])} ({tex_escape(r['case'])}, {point_tex(r['T_K'], r['P_ref_bar'])}); the slowdown goes from "
                        f"{tex_num(r['slowdown_normal'])} to {tex_num(r['slowdown_switched_off'])}."
                        + (' The conclusions do not depend on this species.' if (ok['tof_ratio'] - 1).abs().max() < 0.01 else ''))
        parts.append((f'What-if: reactions of {tex_formula(species)} switched off', rep.series_tables(
            sw, [('tof_ratio', r'rate\\ratio', 'r', None), ('slowdown_normal', r'slowdown\\normal', 'r', None),
                 ('slowdown_switched_off', r'slowdown\\switched off', 'r', None), ('theta_free_normal', r'$\theta_\ast$\\normal', 'r', tex_coverage),
                 ('theta_free_switched_off', r'$\theta_\ast$\\switched off', 'r', tex_coverage),
                 ('top_switched_off', r'most abundant species\\switched off', 'P4.2', tex_escape), ('converged', 'converged', 'l', tex_flag)],
            'Switch-off what-if', 'mechanism_switch_off.csv')))
    top = rep.csv('top_reactions.csv')
    if top is not None and len(top):
        parts.append(('Largest sources and sinks of chosen species', rep.table(
            top, [('case', 'case', 'l', tex_escape), ('T_K', '$T$ (K)', 'r', None), ('P_ref_bar', '$p$ (bar)', 'r', None),
                  ('label', 'species', 'l', tex_formula), ('coverage', 'coverage', 'r', tex_coverage),
                  ('reaction', r'reaction\\index', 'r', None), ('equation', 'equation', 'P5.5', tex_escape),
                  ('rate_1_per_s', r'rate\\(s$^{-1}$)', 'r', None)],
            r'Reactions that make (+) or consume ($-$) chosen species, per site.', 'top_reactions.csv', group_by='label')))
    parts.insert(0, ('Key numbers', bullets(keys)))
    return rep.section('diagnostics', parts)


def section_convergence(rep):
    cv = rep.csv('convergence.csv')
    parts = []
    if cv is not None and len(cv):
        rows = []
        for (case, T, Pb), g in cv.groupby(['case', 'T_K', 'P_ref_bar'], sort=False):
            g = g.sort_values('t_s')
            acc = g['accepted_s'].dropna()
            acc = float(acc.iloc[0]) if len(acc) else np.nan
            at = g[g['accepted_here'].astype(bool)] if 'accepted_here' in g else g.iloc[0:0]
            after = g[g['t_s'] > acc] if np.isfinite(acc) else g.iloc[0:0]
            tof_acc = at['tof_fwd'].iloc[0] if len(at) else np.nan
            err = g['error'].dropna() if 'error' in g else pd.Series(dtype=object)
            rows.append({'case': case, 'T_K': T, 'P_ref_bar': Pb, 'accepted_s': acc,
                         'res_acc': at['norm_residual'].iloc[0] if len(at) else np.nan,
                         'res_end': g['norm_residual'].iloc[-1], 'res_after': after['norm_residual'].max() if len(after) else np.nan,
                         'tof_change': (g['tof_fwd'].iloc[-1] - tof_acc) / tof_acc if np.isfinite(tof_acc) and tof_acc else np.nan,
                         't_end': g['t_s'].iloc[-1], 'error': err.iloc[0] if len(err) else ''})
        summ = pd.DataFrame(rows)
        n_acc = int(summ['accepted_s'].notna().sum())
        n_bad = int((summ['res_after'] > 1).sum())
        parts.append(('What was computed', f'{len(summ)} conditions integrated from a clean surface to '
                                           f'{tex_num(summ["t_end"].max())} s with every checkpoint recorded.'))
        parts.append(('Key numbers', bullets([
            f'{n_acc} of {len(summ)} traced conditions were accepted by the reactor method.',
            (f'{n_bad} of them {"has" if n_bad == 1 else "have"} a residual above 1 after acceptance: '
             f'{"its" if n_bad == 1 else "their"} accepted state was not final.' if n_bad
             else 'After acceptance the residual stays below 1 in every trace: the accepted states are final.'),
            f"Largest change of the dissociation TOF after acceptance: {tex_num(summ['tof_change'].abs().max())} (relative)."])))
        parts.append(('Traced conditions', rep.table(
            summ, [('case', 'case', 'l', tex_escape), ('T_K', '$T$ (K)', 'r', None), ('P_ref_bar', '$p$ (bar)', 'r', None),
                   ('accepted_s', r'accepted\\at (s)', 'r', None), ('res_acc', r'residual at\\acceptance', 'r', None),
                   ('res_after', r'largest residual\\after acceptance', 'r', None), ('res_end', r'residual\\at the end', 'r', None),
                   ('tof_change', r'TOF change\\after acceptance', 'r', None), ('error', 'error', 'P3.5', tex_escape)],
            'One row per traced condition (summary of the checkpoints).', 'convergence.csv', full=True)))
    parts.append(('figure', 'convergence', 'convergence'))
    return rep.section('convergence', parts)


def tolerance_text(rec):
    note = str(rec.get('note') or '')
    if note == 'nan':
        note = ''
    v = rec.get('ppm_of_probe')
    if 'already' in note:
        return r'$<$' + tex_num(v)
    if 'never' in note:
        return r'$>$' + tex_num(v) if v is not None and np.isfinite(v) else 'never'
    return tex_num(v)


def section_curves(rep):
    cfg = rep.cfg
    tol = rep.csv('poisoning_tolerance.csv')
    cur = with_status(rep.csv('poisoning_curves.csv'))
    parts, keys = [], []
    if cur is not None and len(cur):
        imp = cur[cur['case'] != cfg.reference_label]
        levels = sorted(set(imp['ratio'].dropna()))
        n_pts = len(imp[['T_K', 'P_ref_bar']].drop_duplicates())
        parts.append(('What was computed',
                      f'{len(levels)} impurity levels $r$ from {tex_num(min(levels))} to {tex_num(max(levels))} '
                      f'$\\times$ {imp["case"].nunique()} impurities $\\times$ {n_pts} points = {len(imp)} steady states'
                      f' ({int(imp["converged"].astype(bool).sum())} converged, {int((imp["status"] == "two states").sum())} with two steady states).'))
    if tol is not None and len(tol):
        tol = tol.assign(cell=[tolerance_text(r) for r in tol.to_dict('records')])
        for case in [c for c in rep.case_order() if c in set(tol['case'])]:
            d = tol[(tol['case'] == case) & tol['ppm_of_probe'].notna() & (tol['note'].fillna('') == '')]
            never = tol[(tol['case'] == case) & tol['note'].fillna('').str.contains('never')]
            if len(d):
                lo = d.loc[d['ppm_of_probe'].idxmin()]
                keys.append(f"\\textbf{{{tex_escape(case)}}} halves the rate at {tex_num(d['ppm_of_probe'].min())} to "
                            f"{tex_num(d['ppm_of_probe'].max())} ppm of the probe gas (most sensitive point: "
                            f"{point_tex(lo['T_K'], lo['P_ref_bar'])})" + (f'; never halved at {len(never)} points.' if len(never) else '.'))
            elif len(never):
                keys.append(f'\\textbf{{{tex_escape(case)}}} never halves the rate in the range tested.')
        pts = rep.sort_points(map(tuple, tol[['T_K', 'P_ref_bar']].to_numpy()))
        parts.append(('Key numbers', bullets(keys)))
        parts.append(('Tolerated impurity level', rep.pivot(
            tol, 'case', ['T_K', 'P_ref_bar'], 'cell',
            r'Impurity level that halves the dissociation rate, in ppm of the probe gas. $<$: already halved at the lowest level '
            r'(an upper bound); $>$: never halved up to that level.', 'impurity', fmt=lambda s: s if isinstance(s, str) else '--',
            index_fmt=tex_escape, index_order=rep.case_order(), col_order=pts, col_head=lambda c: f'{c[0]:g} K\\\\{c[1]:g} bar',
            source='poisoning_tolerance.csv')))
        parts.append(('', rep.table(rep.sort_cases(tol, ('T_K', 'P_ref_bar')), [
            ('case', 'impurity', 'l', tex_escape), ('series', 'series', 'l', tex_escape), ('T_K', '$T$ (K)', 'r', None),
            ('P_ref_bar', '$p$ (bar)', 'r', None), ('ratio_at_half_rate', r'$r$ at\\half rate', 'r', None),
            ('ppm_of_probe', 'ppm', 'r', None), ('note', 'note', 'P5', tex_escape)],
            'Half-rate levels, every point.', 'poisoning_tolerance.csv', group_by='case')))
    if cur is not None and len(cur):
        imp = cur[(cur['case'] != cfg.reference_label)]
        if len(imp) <= rep.opts['max_rows']:
            parts.append(('Poisoning curves', rep.series_tables(
                imp, [('poisoning_factor', r'$\Pi$', 'r', None), ('theta_impurity', r'$\theta_\mathrm{imp}$', 'r', tex_coverage),
                      ('theta_free', r'$\theta_\ast$', 'r', tex_coverage), ('status', 'status', 'l', tex_escape)],
                'Poisoning curves', 'poisoning_curves.csv', lead=[('ratio', '$r$', 'r', None)], sort_extra=('ratio',))))
        else:
            rows = []
            for (case, T, Pb, series), g in imp.groupby(['case', 'T_K', 'P_ref_bar', 'series'], sort=False):
                g = g.sort_values('ratio')
                c = g[g['converged'].astype(bool)]
                rows.append({'case': case, 'T_K': T, 'P_ref_bar': Pb, 'series': series, 'levels': len(g),
                             'pi_lo': c['poisoning_factor'].iloc[0] if len(c) else np.nan,
                             'pi_hi': c['poisoning_factor'].iloc[-1] if len(c) else np.nan,
                             'th_hi': c['theta_impurity'].iloc[-1] if len(c) else np.nan,
                             'two': int((g['status'] == 'two states').sum()), 'notconv': int((g['status'] == 'not converged').sum())})
            parts.append(('Poisoning curves (summary)', rep.series_tables(
                pd.DataFrame(rows), [('levels', 'levels', 'r', None), ('pi_lo', r'$\Pi$ at the\\lowest level', 'r', None),
                                     ('pi_hi', r'$\Pi$ at the\\highest level', 'r', None),
                                     ('th_hi', r'$\theta_\mathrm{imp}$ at the\\highest level', 'r', tex_coverage),
                                     ('two', r'two-state\\levels', 'r', None), ('notconv', r'not\\converged', 'r', None)],
                'Poisoning curves (summary per impurity and point)', 'poisoning_curves.csv')))
    parts += rep.series_figures('poisoning_curves', 'poisoning_curves')
    parts.append(('figure', 'site_blocking', 'site_blocking'))
    return rep.section('curves', parts)


def sweep_pivots(rep, df, index, label, values, source, flag='two'):
    """Pivot tables (index x case) of several columns of a sweep."""
    out = []
    cases = [c for c in rep.case_order() if c in set(df['case'])]
    for col, what, fmt in values:
        sub = df if col != 'poisoning_factor' else df[df['case'] != rep.cfg.reference_label]
        out.append(rep.pivot(sub, index, 'case', col, what + r' $\dagger$ = two steady states; -- = no steady state.', label,
                             source=source, fmt=fmt, col_order=cases, flag=flag))
    return '\n\n'.join(out)


def section_sweep(rep):
    cfg = rep.cfg
    sw = with_status(rep.csv('temperature_sweep.csv'))
    rel = rep.csv('release_temperature.csv')
    parts, keys = [], []
    if sw is not None and len(sw):
        sw = sw.assign(two=sw['status'] == 'two states')
        sw.loc[~sw['converged'].astype(bool), ['tof_fwd', 'poisoning_factor', 'E_app_kJ_per_mol']] = np.nan
        temps = sorted(sw['T_K'].unique())
        parts.append(('What was computed', f'{len(temps)} temperatures from {temps[0]:g} to {temps[-1]:g} K at '
                                           f'$p$(\\probe) = {cfg.reference_pressure / 1e5:g} bar, {sw["case"].nunique()} cases '
                                           f'({len(sw)} steady states, {int(sw["converged"].astype(bool).sum())} converged).'))
        r = sw[(sw['case'] == cfg.reference_label)]
        if finite(r, 'E_app_kJ_per_mol').size:
            keys.append(f"Reference: apparent activation energy {tex_kj(finite(r, 'E_app_kJ_per_mol').min())} to "
                        f"{tex_kj(finite(r, 'E_app_kJ_per_mol').max())} kJ/mol over the sweep.")
    if rel is not None and len(rel):
        for rec in rel.to_dict('records'):
            note = '' if pd.isna(rec.get('note')) else str(rec['note'])
            keys.append(f"\\textbf{{{tex_escape(rec['case'])}}}: " + (f"rate back to half above {tex_num(rec['release_T_K'])} K"
                                                                      if pd.notna(rec.get('release_T_K')) else tex_escape(note)) + '.')
        parts.append(('Key numbers', bullets(keys)))
        parts.append(('Release temperatures', rep.table(
            rep.sort_cases(rel), [('case', 'impurity case', 'l', tex_escape), ('release_T_K', r'release\\$T$ (K)', 'r', None),
                                  ('note', 'note', 'P7', tex_escape)],
            r'Release temperature: the lowest temperature at which the rate is back to half ($\Pi \geq 0.5$).', 'release_temperature.csv')))
    elif keys:
        parts.append(('Key numbers', bullets(keys)))
    if sw is not None and len(sw):
        parts.append(('Rates against temperature', sweep_pivots(rep, sw, 'T_K', '$T$ (K)', [
            ('tof_fwd', r'Dissociation TOF (s$^{-1}$) per temperature and case.', None),
            ('poisoning_factor', r'Poisoning factor $\Pi$ per temperature and impurity case.', None),
            ('E_app_kJ_per_mol', r'Apparent activation energy $E_\mathrm{app}$ (kJ/mol).', tex_kj)], 'temperature_sweep.csv')))
    parts += [('figure', 'temperature_sweep', 'temperature_sweep'), ('figure', 'sweep_coverages', 'sweep_coverages')]
    return rep.section('sweep', parts)


def section_arrhenius(rep):
    fits = rep.csv('arrhenius_fits.csv')
    intr = rep.csv('arrhenius_intrinsic.csv')
    loc = rep.csv('arrhenius_local.csv')
    cfg = rep.cfg
    parts, keys = [], []
    if rep.model is not None:
        from .analyses.arrhenius import mechanism_parameters
        keys.append(r'The probe reaction in the mechanism: \texttt{' + tex_escape(mechanism_parameters(rep.model)) + '}.')
    if fits is not None and len(fits):
        fits = fits.assign(expr=[tex_arrhenius(a, e) for a, e in zip(fits['A'], fits['Ea_kJ_per_mol'])],
                           units=['--' if ('K =' in q or '/' in q) else 's$^{-1}$' for q in fits['quantity']],
                           T_range=[f'{a:g}--{b:g}' for a, b in zip(fits['T_min_K'], fits['T_max_K'])])
        intrinsic = fits[fits['level'] == 'intrinsic']
        for rec in intrinsic.to_dict('records'):
            if 'K =' in rec['quantity']:
                keys.append(f"Equilibrium of the \\probe{{}} step: $K = $ {rec['expr']}, i.e.\\ "
                            f"$\\Delta H^\\circ$ = {tex_kj(rec['dH_kJ_per_mol'])} kJ/mol, $\\Delta S^\\circ$ = "
                            f"{tex_num(rec['dS_J_per_mol_K'])} J/(mol K) ({rec['T_range']} K).")
            else:
                name = 'dissociation on the clean metal, $k_\\mathrm{diss}$' if 'diss' in rec['quantity'] else \
                    'recombination from a full surface, $k_\\mathrm{rec}$'
                keys.append(f"Intrinsic {name} = {rec['expr']} s$^{{-1}}$ per site ({rec['T_range']} K; "
                            f"{tex_escape(rec['conditions'])}).")
        obs = fits[fits['level'] == 'observed']
        diss = obs[obs['quantity'].str.startswith('dissociation')]
        ref = diss[diss['case'] == cfg.reference_label]
        if len(ref):
            keys.append(f"Observed dissociation on the steady-state surface, probe gas alone: {ref['expr'].iloc[0]} s$^{{-1}}$ "
                        f"($E_a$ = {tex_kj(ref['Ea_kJ_per_mol'].iloc[0])} kJ/mol).")
        imp = diss[diss['case'] != cfg.reference_label].dropna(subset=['dEa_vs_reference_kJ_per_mol'])
        if len(imp):
            imp = imp.sort_values('dEa_vs_reference_kJ_per_mol', ascending=False)
            keys.append('Impurity effect on the dissociation $E_a$ ($\\Delta E_a$, kJ/mol; $A/A_\\mathrm{ref}$): '
                        + '; '.join(f"{tex_escape(r['case'])} {tex_kj(r['dEa_vs_reference_kJ_per_mol'])} ({tex_num(r['A_over_reference'])})"
                                    for r in imp.to_dict('records')) + '.')
        curved = obs[obs['curved'].astype(str).str.lower() == 'true']
        if len(curved):
            keys.append('Curved (one $E_a$ does not fit the whole range; see the local values): '
                        + ', '.join(sorted({tex_escape(c) for c in curved['case']})) + '.')
        parts.append(('Key numbers', bullets(keys)))
        cols = [('quantity', 'quantity', 'P3.3', tex_escape), ('case', 'case', 'l', tex_escape), ('expr', r'$A\,\exp(-E_a/RT)$', 'l', lambda s: s),
                ('units', 'units', 'l', lambda s: s), ('T_range', '$T$ (K)', 'l', tex_escape), ('conditions', 'conditions', 'P4.2', tex_escape),
                ('R2', '$R^2$', 'r', lambda v: tex_num(v, 4)), ('dEa_vs_reference_kJ_per_mol', r'$\Delta E_a$\\(kJ/mol)', 'r', tex_kj),
                ('A_over_reference', r'$A/A_\mathrm{ref}$', 'r', None), ('curved', 'curved', 'l', tex_flag)]
        parts.append(('Arrhenius expressions', rep.table(
            pd.concat([intrinsic, rep.sort_cases(obs)]), cols,
            r'Rates of hydrogen in Arrhenius form, $E_a$ in kJ/mol, $R = 8.314\times10^{-3}$ kJ/(mol K). Intrinsic rows: the '
            r'\probe{} step alone; observed rows: the steady-state surface of each case.', 'arrhenius_fits.csv', group_by='case')))
    elif keys:
        parts.append(('Key numbers', bullets(keys)))
    if intr is not None and len(intr):
        parts.append(('Intrinsic rate constants', rep.table(
            intr, [('T_K', '$T$ (K)', 'r', None), ('k_diss_clean_1_per_s', r'$k_\mathrm{diss}$ (s$^{-1}$)\\empty surface', 'r', None),
                   ('k_rec_full_1_per_s', r'$k_\mathrm{rec}$ (s$^{-1}$)\\full of \adsorbate{}', 'r', None),
                   ('K_eq', r'$K$\\(dimensionless)', 'r', None), ('dH_thermo_kJ_per_mol', r'$\Delta H^\circ$\\(kJ/mol)', 'r', tex_kj),
                   ('dS_thermo_J_per_mol_K', r'$\Delta S^\circ$\\(J/(mol K))', 'r', None)],
            r'Intrinsic constants of the \probe{} step at each temperature; $\Delta H^\circ$, $\Delta S^\circ$ from the mechanism thermochemistry.',
            'arrhenius_intrinsic.csv')))
    if loc is not None and len(loc):
        cases = [c for c in rep.case_order() if c in set(loc['case'])]
        parts.append(('Local activation energies', '\n\n'.join([
            rep.pivot(loc, 'T_K', 'case', 'Ea_local_fwd_kJ_per_mol', r'Local activation energy of dissociation (kJ/mol), '
                      r'$-R\,d\ln\mathrm{TOF}/d(1/T)$ by central differences.', '$T$ (K)', 'arrhenius_local.csv', fmt=tex_kj, col_order=cases),
            rep.pivot(loc, 'T_K', 'case', 'fwd_over_rev', r'Dissociation / recombination TOF: 1 = the \probe{} step at equilibrium.',
                      '$T$ (K)', 'arrhenius_local.csv', fmt=lambda v: tex_num(v, 4), col_order=cases)])))
    parts.append(('figure', 'arrhenius', 'arrhenius'))
    return rep.section('arrhenius', parts)


def section_pressure_sweep(rep):
    ps = with_status(rep.csv('pressure_sweep.csv'))
    half = rep.csv('half_rate_pressure.csv')
    parts, keys, tables = [], [], []
    if ps is not None and len(ps):
        ps = ps.assign(two=ps['status'] == 'two states')
        ps.loc[~ps['converged'].astype(bool), ['tof_fwd', 'poisoning_factor', 'slope_dlnTOF_dlnp']] = np.nan
        temps = sorted(ps['T_K'].unique())
        Ps = sorted(ps['P_ref_bar'].unique())
        parts.append(('What was computed', f'{len(Ps)} probe-gas pressures from {tex_num(Ps[0])} to {tex_num(Ps[-1])} bar at '
                                           f'{", ".join(f"{T:g}" for T in temps)} K, {ps["case"].nunique()} cases, every impurity at '
                                           f'its fixed ratio to the probe gas ({len(ps)} steady states).'))
        for T in temps:
            r = ps[(ps['T_K'] == T) & (ps['case'] == rep.cfg.reference_label)]
            s = finite(r, 'slope_dlnTOF_dlnp')
            if s.size:
                keys.append(f'{T:g} K, probe gas alone: order in \\probe{{}} from {tex_num(s.iloc[0])} (lowest pressure) to '
                            f'{tex_num(s.iloc[-1])} (highest).')
            d = ps[ps['T_K'] == T]
            tables.append(sweep_pivots(rep, d, 'P_ref_bar', '$p$(\\probe) (bar)', [
                ('tof_fwd', f'Dissociation TOF (s$^{{-1}}$) at {T:g} K.', None),
                ('poisoning_factor', f'Poisoning factor $\\Pi$ at {T:g} K.', None),
                ('slope_dlnTOF_dlnp', f'Local slope $d\\ln\\mathrm{{TOF}}/d\\ln p$ at {T:g} K.', None)], 'pressure_sweep.csv'))
    if half is not None and len(half):
        for rec in half.to_dict('records'):
            if pd.notna(rec.get('half_rate_P_bar')):
                keys.append(f"\\textbf{{{tex_escape(rec['case'])}}} at {rec['T_K']:g} K: $\\Pi$ crosses 0.5 at "
                            f"{tex_num(rec['half_rate_P_bar'])} bar ({tex_escape(rec['direction'])}).")
    parts.append(('Key numbers', bullets(keys)))
    if half is not None and len(half):
        parts.append(('Half-rate pressures', rep.table(
            rep.sort_cases(half, ('T_K',)), [('T_K', '$T$ (K)', 'r', None), ('case', 'case', 'l', tex_escape),
                                             ('half_rate_P_bar', r'half-rate\\pressure (bar)', 'r', None),
                                             ('direction', 'direction', 'l', tex_escape), ('note', 'note', 'P6', tex_escape)],
            r'Probe-gas pressure where $\Pi$ crosses 0.5.', 'half_rate_pressure.csv')))
    parts.append(('Rates against pressure', '\n\n'.join(tables)))
    parts.append(('figure', 'pressure_sweep', 'pressure_sweep'))
    for stem in rep.stems(r'pressure_sweep_coverages_.+'):
        parts.append(('figure', 'pressure_sweep_coverages', stem, f'({stem.rsplit("_", 1)[-1].replace("K", " K")}.)'))
    return rep.section('pressure_sweep', parts)


def section_orders(rep):
    od = rep.csv('reaction_orders.csv')
    parts, keys = [], []
    if od is not None and len(od):
        ref = od[od['case'] == rep.cfg.reference_label]
        if finite(ref, 'order_probe').size:
            keys.append(f"Probe gas alone: order in \\probe{{}} {tex_num(ref['order_probe'].min())} to {tex_num(ref['order_probe'].max())}.")
        for case in [c for c in rep.case_order() if c in set(od['case']) and c != rep.cfg.reference_label]:
            d = od[od['case'] == case]
            if finite(d, 'order_impurity').size:
                keys.append(f"\\textbf{{{tex_escape(case)}}}: order in the impurity {tex_num(d['order_impurity'].min())} to "
                            f"{tex_num(d['order_impurity'].max())}; in \\probe{{}} {tex_num(d['order_probe'].min())} to "
                            f"{tex_num(d['order_probe'].max())}.")
        parts.append(('Key numbers', bullets(keys)))
        parts.append(('Orders at every point', rep.series_tables(
            od, [('tof_fwd_base', r'diss. TOF\\(s$^{-1}$)', 'r', None), ('order_probe', r'order in\\\probe{}', 'r', None),
                 ('order_impurity', r'order in the\\impurity', 'r', None), ('all_starts_agree', r'starts\\agree', 'l', tex_flag)],
            'Reaction orders of the dissociation rate', 'reaction_orders.csv')))
    parts += rep.series_figures('reaction_orders', 'reaction_orders')
    return rep.section('orders', parts)


def is_probe_step(rep, item, label):
    """Is a DRC row the probe reaction? By index when the mechanism is loaded, else by its species."""
    if rep.model is not None:
        return str(item) == str(rep.model.i_probe)
    sides = re.split(r'<?=>', str(label))
    return len(sides) == 2 and any(rep.cfg.probe_gas in side.split() for side in sides) and rep.cfg.probe_adsorbate in str(label)


def section_rate_control(rep):
    sens = rep.csv('sensitivity.csv')
    parts, keys = [], []
    if sens is not None and len(sens):
        rows = []
        name = lambda r: tex_formula(str(r['label']).split('  [')[0])                              # noqa: E731
        step = lambda d, k: (f"\\#{d['item'].iloc[k]}: {tex_escape(d['label'].iloc[k])} ({tex_num(d['X_tof_fwd'].iloc[k])})"
                             if len(d) > k else '--')                                              # noqa: E731
        for (series, T, Pb, case), g in sens.groupby(['series', 'T_K', 'P_ref_bar', 'case'], sort=False):
            trc = g[g['kind'] == 'TRC'].sort_values('X_tof_fwd')
            drc = g[g['kind'] == 'DRC'].assign(_a=lambda d: d['X_tof_fwd'].abs()).sort_values('_a', ascending=False)
            probe = drc[[is_probe_step(rep, i, lab) for i, lab in zip(drc['item'], drc['label'])]]
            other = drc.drop(probe.index)
            rows.append({'series': series, 'T_K': T, 'P_ref_bar': Pb, 'case': case,
                         'trc1': f"{name(trc.iloc[0])} ({tex_num(trc['X_tof_fwd'].iloc[0])})" if len(trc) else '--',
                         'trc2': f"{name(trc.iloc[1])} ({tex_num(trc['X_tof_fwd'].iloc[1])})" if len(trc) > 1 else '--',
                         'trc_sum': trc['X_tof_fwd'].sum() if len(trc) else np.nan,
                         'drc_probe': probe['X_tof_fwd'].iloc[0] if len(probe) else np.nan,
                         'drc1': step(other, 0), 'drc2': step(other, 1)})
        summ = pd.DataFrame(rows)
        for case in [c for c in rep.case_order() if c in set(summ['case'])]:
            holders = summ[summ['case'] == case]['trc1'].str.split(' \\(').str[0].value_counts()
            keys.append(f"\\textbf{{{tex_escape(case)}}}: main site holder (most negative TRC) "
                        + ', '.join(f'{h} at {n} point{"s" if n > 1 else ""}' for h, n in holders.items()) + '.')
        dp = finite(summ, 'drc_probe')
        if dp.size:
            keys.append('DRC of the probe step: '
                        + (f'{tex_num(dp.min())} at every point' if dp.max() - dp.min() < 0.005 else
                           f'{tex_num(dp.min())} to {tex_num(dp.max())}')
                        + (' (it limits the rate everywhere).' if dp.min() > 0.9 else '.'))
        parts.append(('Key numbers', bullets(keys)))
        parts.append(('Rate control at every point', rep.series_tables(
            summ, [('trc1', r'top TRC species', 'l', lambda s: s), ('trc2', r'second', 'l', lambda s: s),
                   ('trc_sum', r'TRC\\sum', 'r', None), ('drc_probe', r'DRC of the\\probe step', 'r', None),
                   ('drc1', r'largest other DRC step', 'P4.4', lambda s: s), ('drc2', 'second other DRC step', 'P4.4', lambda s: s)],
            r'Rate control of the dissociation rate (values in parentheses)', 'sensitivity.csv')))
    pts = []
    for stem in rep.stems(r'rate_control_[\d.]+K_[0-9.e+-]+bar(_DRC)?'):
        m = re.fullmatch(r'rate_control_([\d.]+)K_([0-9.e+-]+)bar(_DRC)?', stem)
        pts.append((float(m.group(1)), float(m.group(2))))
    for T, Pb in rep.sort_points(pts):
        stem = f'rate_control_{T:.0f}K_{Pb:g}bar'
        parts.append(('figure', 'rate_control', stem, f'({point_tex(T, Pb)}.)'))
        parts.append(('figure', 'rate_control_DRC', stem + '_DRC', f'({point_tex(T, Pb)}.)'))
    return rep.section('rate_control', parts)


def path_diagrams(rep):
    """[(element, case, T, P bar, file)] of the reaction-path diagrams, .pdf preferred over .png."""
    found = {}
    slug = {c.replace(' ', ''): c for c in rep.case_order()}
    for path in sorted(glob.glob(os.path.join(rep.out, 'reaction_paths', 'path_*.p[dn][fg]'))):
        m = PATH_RE.fullmatch(os.path.basename(path))
        if not m:
            continue
        key = (m['el'], slug.get(m['case'], m['case']), float(m['T']), float(m['P']))
        if key not in found or path.endswith('.pdf'):
            found[key] = os.path.relpath(path, rep.out)
    order = {c: i for i, c in enumerate(rep.case_order())}
    pts = {pt: i for i, pt in enumerate(rep.point_order())}
    return sorted(((*k, f) for k, f in found.items()), key=lambda d: (pts.get((d[2], d[3]), 99), d[0], order.get(d[1], 99)))


def paths_reference_point(rep, diagrams):
    P_ref = rep.cfg.reference_pressure / 1e5
    at_ref = sorted({(d[2], d[3]) for d in diagrams if np.isclose(d[3], P_ref)})
    return at_ref[0] if at_ref else (diagrams[0][2], diagrams[0][3])


def section_paths(rep):
    diagrams = path_diagrams(rep)
    if not diagrams:
        return rep.section('paths', [])
    mode = rep.opts['path_diagrams']
    ref = paths_reference_point(rep, diagrams)
    here = diagrams if mode == 'all' else [d for d in diagrams if (d[2], d[3]) == ref]
    elsewhere = [d for d in diagrams if d not in here]
    count = pd.DataFrame([{'T_K': d[2], 'P_ref_bar': d[3], 'element': d[0]} for d in diagrams]) \
        .groupby(['T_K', 'P_ref_bar', 'element'], sort=False).size().reset_index(name='n')
    count['where'] = ['this section' if mode == 'all' or (t, p) == ref else
                      ('Appendix~\\ref{app:paths}' if mode == 'reference+appendix' else 'not shown')
                      for t, p in zip(count['T_K'], count['P_ref_bar'])]
    parts = [('What was computed', f'{len(diagrams)} diagrams: every case at {len(count)} point(s). '
                                   + ('All are shown here.' if mode == 'all' else
                                      f'The diagrams at {point_tex(*ref)} are shown here'
                                      + (', the others in Appendix~\\ref{app:paths}.' if elsewhere and mode == 'reference+appendix' else '.'))),
             ('', rep.table(count, [('T_K', '$T$ (K)', 'r', None), ('P_ref_bar', '$p$ (bar)', 'r', None),
                                    ('element', 'element', 'l', tex_escape), ('n', 'diagrams', 'r', None), ('where', 'shown in', 'l', lambda s: s)],
                            'Reaction-path diagrams of this run (folder \\texttt{reaction\\_paths}).', full=True))]
    for el, case, T, Pb, f in here:
        parts.append(('figure', 'reaction_paths', os.path.splitext(f)[0], f'{tex_escape(el)}, {tex_escape(case)}, {point_tex(T, Pb)}.', 0.4))
    rep.paths_elsewhere = elsewhere if mode == 'reference+appendix' else []
    return rep.section('paths', parts)


def section_transients(rep):
    pot = rep.csv('poisoning_over_time.csv')
    rep.csv('poisoning_over_time_trajectories.csv')
    parts, keys = [], []
    if pot is not None and len(pot):
        parts.append(('What was computed', f'{len(pot)} transients: each starts from the steady state in the probe gas alone; '
                                           f'the impurity is added at $t = 0$ at the levels {", ".join(tex_num(r) for r in sorted(set(pot["ratio"])))}. '
                                           f'The full curves are in {tex_file("poisoning_over_time_trajectories.csv")}.'))
        for case in [c for c in rep.case_order() if c in set(pot['case'])]:
            d = pot[(pot['case'] == case) & pot['t_half_s'].notna()]
            if len(d):
                r = d.loc[d['t_half_s'].idxmin()]
                keys.append(f"\\textbf{{{tex_escape(case)}}}: fastest half time {tex_num(r['t_half_s'])} s "
                            f"($r$ = {tex_num(r['ratio'])}, {point_tex(r['T_K'], r['P_ref_bar'])}); deepest drop to "
                            f"{tex_num(pot[pot['case'] == case]['min_over_t0'].min())} of the starting rate.")
            else:
                keys.append(f'\\textbf{{{tex_escape(case)}}}: the rate never halves in the transients.')
        parts.append(('Key numbers', bullets(keys)))
        cols = [('t_half_s', r'half time\\(s)', 'r', None), ('t_tenth_s', r'tenth time\\(s)', 'r', None),
                ('min_over_t0', r'deepest\\(of $t=0$)', 'r', None), ('end_over_t0', r'final\\(of $t=0$)', 'r', None),
                ('t_settle_s', r'settle\\time (s)', 'r', None), ('top_species_end', r'most abundant\\species at the end', 'P4', tex_escape)]
        if pot['error'].notna().any():
            cols.append(('error', 'error', 'P3', tex_escape))
        parts.append(('Poisoning times', rep.series_tables(pot, cols, 'Poisoning over time', 'poisoning_over_time.csv',
                                                           lead=[('ratio', '$r$', 'r', None)], sort_extra=('ratio',))))
    for stem in rep.stems(r'poisoning_over_time_.+'):
        parts.append(('figure', 'poisoning_over_time', stem, f'({tex_escape(stem.replace("poisoning_over_time_", ""))}.)'))
    return rep.section('transients', parts)


def section_flux(rep):
    parts, rows = [], []
    for folder in sorted(glob.glob(os.path.join(rep.out, 'flux', '*')), key=lambda p: (len(p), p)):
        if not os.path.isdir(folder):
            continue
        frames = sorted(glob.glob(os.path.join(folder, 'flux_diagram_*.png')))
        avi = os.path.join(folder, 'flux_diagram.avi')
        name = os.path.basename(folder)
        rows.append({'reactor': name, 'frames': len(frames),
                     'video': os.path.relpath(avi, rep.out) if os.path.exists(avi) else '--'})
        if frames:
            mid = frames[len(frames) // 2]
            parts.append(('figure', 'flux', os.path.splitext(os.path.relpath(mid, rep.out))[0],
                          f'(reactor {tex_escape(name)}, frame {len(frames) // 2 + 1} of {len(frames)}).', 0.42))
    parts.insert(0, ('Reactors', rep.table(pd.DataFrame(rows), [('reactor', 'reactor', 'l', tex_escape), ('frames', 'frames', 'r', None),
                                                               ('video', 'video', 'l', tex_file)],
                                           'Flux-diagram folders (one per RMG reactor).', full=True) if rows else ''))
    return rep.section('flux', parts)


SECTION_BUILDERS = {
    'steady_state': section_steady_state, 'composition': section_composition, 'diagnostics': section_diagnostics,
    'convergence': section_convergence, 'curves': section_curves, 'sweep': section_sweep, 'arrhenius': section_arrhenius,
    'pressure_sweep': section_pressure_sweep, 'orders': section_orders, 'rate_control': section_rate_control,
    'paths': section_paths, 'transients': section_transients, 'flux': section_flux,
}


# --- Appendices --------------------------------------------------------------------------
def appendix_paths(rep):
    diagrams = getattr(rep, 'paths_elsewhere', [])
    if not diagrams:
        return ''
    out = [r'\section{Reaction-path diagrams at the other points}', r'\label{app:paths}',
           'How to read them: see the end of Section~\\ref{sec:paths}.']
    last = None
    for el, case, T, Pb, f in diagrams:
        if (T, Pb) != last:
            out.append(r'\subsection*{' + point_tex(T, Pb) + '}')
            last = (T, Pb)
        tex, _ = rep.figure(os.path.splitext(f)[0], RT.SECTIONS['paths']['figures'][0]['caption']
                            + f' {tex_escape(el)}, {tex_escape(case)}, {point_tex(T, Pb)}.', 0.4)
        out.append(tex)
    return '\n\n'.join(out)


def appendix_files(rep):
    """Every file of the output folder: tables with their size, figures, folders by count."""
    owner = {}
    for name in SECTION_OUTPUTS:
        for f in rep.outputs(name):
            owner[f] = name
    rows = []
    for entry in sorted(os.listdir(rep.out)):
        path = os.path.join(rep.out, entry)
        if entry in ('report',) or entry.startswith('.'):
            continue
        if os.path.isdir(path):
            n = sum(len(fs) for _, _, fs in os.walk(path))
            what = 'cached results' if entry.startswith('cache') else 'folder'
            sec = next((s for f, s in owner.items() if f.startswith(entry + os.sep)), None)
            rows.append({'file': entry + '/', 'kind': what, 'size': f'{n} files', 'section': sec, 'shown': '--'})
            continue
        ext = os.path.splitext(entry)[1].lower()
        size = ''
        if ext == '.csv':
            try:
                df = pd.read_csv(path)
                size = f'{len(df)} rows $\\times$ {len(df.columns)} columns'
            except (pd.errors.ParserError, pd.errors.EmptyDataError, UnicodeDecodeError):
                size = 'unreadable'
        stem = os.path.splitext(entry)[0]
        if ext == '.png' and os.path.exists(os.path.join(rep.out, stem + '.pdf')):
            continue                                            # listed with its .pdf
        if ext == '.pdf' and os.path.exists(os.path.join(rep.out, stem + '.png')):
            rows.append({'file': entry + ', .png', 'kind': 'figure', 'size': '', 'section': owner.get(entry),
                         'shown': 'yes' if {entry, stem + '.png'} & rep.used else 'no'})
            continue
        kind = {'.csv': 'table', '.png': 'figure', '.pdf': 'figure', '.yaml': 'mechanism', '.log': 'log'}.get(ext, 'other')
        rows.append({'file': entry, 'kind': kind, 'size': size, 'section': owner.get(entry),
                     'shown': 'yes' if entry in rep.used else 'no'})
    df = pd.DataFrame(rows)
    df['section'] = [(r'\ref{sec:' + s + '}') if isinstance(s, str) and s in rep.built else '--' for s in df['section']]
    return '\n\n'.join([r'\section{Output files}', r'\label{app:files}',
                        'Every file in the output folder, ' + tex_path(rep.out) + '. A figure is shown in the report as its '
                        'PDF when the run wrote one (publication style), else as its PNG.',
                        rep.table(df, [('file', 'file', 'P6.5', tex_escape), ('kind', 'kind', 'l', tex_escape),
                                       ('size', 'size', 'l', lambda s: s), ('section', 'section', 'l', lambda s: s),
                                       ('shown', 'shown', 'l', tex_escape)],
                                  'Contents of the output folder.', full=True)])


def appendix_settings(rep):
    os.makedirs(rep.dir, exist_ok=True)
    with open(rep.cfg.settings_file, errors='replace') as f:
        text = ascii_text(f.read())
    with open(os.path.join(rep.dir, 'settings_listing.py'), 'w') as f:
        f.write(text)
    return '\n\n'.join([r'\section{The settings of this run}', r'\label{app:settings}',
                        r'The run\textquotesingle s \texttt{settings.py} as it is now (non-ASCII characters transliterated). '
                        'Settings not in this file take the package defaults.',
                        r'\lstinputlisting{settings_listing.py}'])


# --- Build -------------------------------------------------------------------------------
def latexmk_path(opts):
    for p in (opts.get('latexmk'), shutil.which('latexmk'), '/Library/TeX/texbin/latexmk'):
        if p and os.path.exists(p):
            return p
    return None


def latex_errors(log_path, n=12):
    if not os.path.exists(log_path):
        return 'no LaTeX log'
    with open(log_path, errors='replace') as f:
        lines = f.read().splitlines()
    hits = [i for i, line in enumerate(lines) if line.startswith('!') or ':' in line[:200] and 'Error' in line]
    return '\n'.join('\n'.join(lines[i:i + 3]) for i in hits[:n]) or '\n'.join(lines[-n:])


def write_tex(rep):
    title, title_tex = title_page(rep)
    rep.built = [name for name in ANALYSIS_ORDER if name in SECTION_BUILDERS and rep.outputs(name)]
    sections = []
    for name in rep.built:
        sections += [r'\clearpage', SECTION_BUILDERS[name](rep)]
    appendices = [a for a in (appendix_paths(rep), appendix_files(rep), appendix_settings(rep)) if a]
    parts = [preamble(rep, re.sub(r'\\[a-z]+\{([^}]*)\}', r'\1', title)), r'\begin{document}', title_tex,
             r'\tableofcontents', r'\clearpage', run_summary(rep, rep.built), *sections,
             r'\clearpage', r'\appendix', *[r'\clearpage' + '\n' + a for a in appendices], r'\end{document}']
    os.makedirs(rep.dir, exist_ok=True)
    stem = run_name(rep.cfg).replace('/', '_') + '_report'
    path = os.path.join(rep.dir, stem + '.tex')
    with open(path, 'w') as f:
        f.write('\n\n'.join(parts) + '\n')
    return path


def build(cfg_or_path, model=None, log=print, load_model=True, latex=True):
    """Write output/report/<run>_report.tex and build the PDF from it. model: a loaded Model
    (for the mechanism facts and warnings); None loads one when load_model is True (that only
    reads the mechanism: nothing is solved). latex=False writes the .tex only.
    Returns the PDF path, or the .tex path when the PDF is not built."""
    cfg = load_config(cfg_or_path) if isinstance(cfg_or_path, str) else cfg_or_path
    if model is None and load_model:
        from .model import Model
        model = Model(cfg)
    rep = Report(cfg, model, log)
    tex = write_tex(rep)
    exe = latexmk_path(rep.opts) if latex else None
    if exe is None:
        log(f'report: wrote {tex}' + ('' if not latex else
                                       f'; LaTeX (latexmk) not found. Build it with: latexmk -pdf {os.path.basename(tex)}'))
        return tex
    stem = os.path.splitext(os.path.basename(tex))[0]
    env = dict(os.environ, PATH=os.path.dirname(exe) + os.pathsep + os.environ.get('PATH', ''))
    cmd = [exe, '-pdf', '-interaction=nonstopmode', '-halt-on-error', '-file-line-error', os.path.basename(tex)]
    log(f'report: building {stem}.pdf with LaTeX ...')
    try:
        proc = subprocess.run(cmd, cwd=rep.dir, capture_output=True, text=True, timeout=LATEX_TIMEOUT, env=env)
    except subprocess.TimeoutExpired:
        raise ReportError(f'LaTeX took longer than {LATEX_TIMEOUT} s on {tex}') from None
    pdf = os.path.join(rep.dir, stem + '.pdf')
    if proc.returncode != 0 or not os.path.exists(pdf):
        raise ReportError(f'LaTeX failed on {tex}:\n' + latex_errors(os.path.join(rep.dir, stem + '.log')))
    subprocess.run([exe, '-c', os.path.basename(tex)], cwd=rep.dir, capture_output=True, env=env)   # keep .tex and .pdf
    log(f'report: {pdf}')
    return pdf
