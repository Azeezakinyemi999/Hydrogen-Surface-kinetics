r"""Fixed text of the run report (LaTeX snippets, not escaped).

\probe{} and \adsorbate{} are defined by the report from the run (e.g. H$_2$ and H*).
Each section has a 'title', a 'goal', 'columns' (what each table column means) and
'figures': one entry per kind of figure ('key', used by report.py) with a 'caption' and its
'panels': (panel, what it shows, how to read it, [(what you see, what it
means), ...]). Wording follows docs/simulation_workflow.md, "Reading the results".
"""

CONVENTIONS = r"""
\begin{itemize}[leftmargin=*]
\item \textbf{Rates} are turnover frequencies (TOF): events per site per second (s$^{-1}$).
  \emph{Dissociation} (forward) is the probe step in the adsorption direction,
  \probe{} + 2$\ast$ $\rightarrow$ 2\,\adsorbate{}; \emph{recombination} (reverse) is
  2\,\adsorbate{} $\rightarrow$ \probe{} + 2$\ast$; \emph{net} = dissociation $-$ recombination,
  close to 0 when the step is at equilibrium.
\item \textbf{Poisoning factor} $\Pi$ = dissociation TOF with the impurity / dissociation TOF of
  the probe gas alone, at the same temperature and probe-gas pressure. $\Pi = 1$: no effect;
  $\Pi = 0.5$: the rate is halved; $\Pi = 10^{-6}$: a million times slower.
  The \emph{slowdown} is $1/\Pi$ (10 = ten times slower).
\item \textbf{Impurity level} $r = p_\mathrm{impurity}/p_\mathrm{probe}$; $r = 10^{-6}$ is 1\,ppm
  of the probe gas. The impurities are \emph{added} to the probe gas: the probe-gas pressure is
  the same with and without them.
\item \textbf{Coverages} $\theta$ are fractions of the surface sites: $\theta_\ast$ empty,
  $\theta_\mathrm{H}$ the probe adsorbate (\adsorbate{}), $\theta_\mathrm{H_2}$ the molecularly
  adsorbed probe gas (if the mechanism has it) and $\theta_\mathrm{imp}$ everything else, i.e.\
  sites held by the impurities and what they turn into.
\item \textbf{Two series.} Temperature and pressure are never varied together. The
  \emph{temperature series} varies $T$ at the reference probe-gas pressure; each
  \emph{pressure series} varies the probe-gas pressure at one fixed $T$ (the impurities keep
  their ratio to it). Figures and tables marked ``vs T'' belong to the temperature series,
  ``vs P'' to the pressure series. A point in both is labeled ``both''.
\item \textbf{Two starting surfaces.} Every steady state is solved from a clean surface (start A)
  and from a covered one (start B). If they end in different steady states the point has
  \emph{two steady states} (history dependence): an open ring in the figures, ``two states'' in
  the tables. Values are those of start A unless stated.
\item \textbf{Markers.} An open diamond on the 0.5 line marks where the rate is halved; a gap in a
  line is a point without a converged steady state; ``--'' in a table is a missing value.
\item \textbf{Colors} follow the case or the species, never the rank: an impurity or an
  adsorbate has the same color in every figure of this report.
\end{itemize}
"""

SECTIONS = {
    'steady_state': {
        'title': 'Main steady states: dissociation rate and poisoning',
        'goal': r"""How fast does \probe{} dissociate on the metal, per site, at each point of the
temperature and pressure series, and how much slower is it when each impurity (alone, or all
together) is in the gas? Each condition is the steady state of the surface under a gas of fixed
composition: the probe gas at its pressure plus the impurities at their ratios.""",
        'columns': [
            ('case', r'the gas: the probe gas alone (the reference), one impurity added, or all impurities added'),
            (r'$p_\mathrm{imp}$', r'total impurity partial pressure (bar), $r \times p$(\probe{}) summed over the impurities in the case'),
            ('diss. / rec. / net TOF', r'dissociation, recombination and net rate of the probe step (s$^{-1}$); net $\ll$ diss.\ means the step is at equilibrium'),
            (r'$\theta_\ast$, $\theta_\mathrm{H}$, $\theta_\mathrm{H_2}$, $\theta_\mathrm{imp}$', 'fractions of sites: empty, probe adsorbate, molecular probe gas, impurity-derived'),
            ('site factor', r'the number of sites Cantera implies, per nominal site (1 with one-site adsorbates only, up to 2 with two-site ones); used for the fixed-site-count estimate'),
            (r'$\Pi$, $1/\Pi$', r'poisoning factor and slowdown relative to the reference at the same point'),
            (r'$\Pi_\mathrm{B}$', r'the poisoning factor of the steady state reached from start B (equal to $\Pi$ unless there are two steady states)'),
            ('status', r'ok; two states (starts A and B disagree); not converged'),
        ],
        'figures': [
            {'key': 'rates_summary_vs_T', 'caption': 'Main steady states along the temperature series: one marker per temperature.',
             'panels': [
                 ('(a) dissociation rate',
                  r'The forward TOF of \probe{} dissociation for each case (x axis); one marker shape and color per temperature.',
                  r'Log scale: one decade up is ten times faster. Compare a case with the reference case at the same temperature (same marker).',
                  [('a marker far below the reference marker', 'that impurity slows dissociation at that temperature'),
                   ('markers going up with temperature', 'dissociation is activated, or an impurity desorbs and frees sites as T rises'),
                   ('markers going down with temperature', r'fewer empty sites at high T (an impurity product forms and binds), or a sticking coefficient that falls with T'),
                   ('markers that hardly change with temperature', r'the rate is set by how often \probe{} hits empty sites, not by a barrier')]),
                 ('(b) slowdown',
                  r'How many times slower than the probe gas alone (filled markers). Open markers: the same if the number of sites were fixed.',
                  r"The dashed line is ``no effect''; 10 means ten times slower, 1000 a thousand times.",
                  [('a filled marker above the line', 'the impurity poisons at that temperature'),
                   ('the slowdown coming down as T rises', 'the poison desorbs or reacts away at high temperature: poisoning lifts'),
                   ('the slowdown going up as T rises', 'the impurity (or a product of it) binds more strongly at high temperature, e.g.\\ after it decomposes'),
                   ('an open marker far from the filled one', 'multi-site adsorbates change the site count; the size of the effect depends on that convention')]),
                 ('(c) empty sites',
                  r'The fraction of sites that are empty, $\theta_\ast$.',
                  r'Dissociation needs two empty sites, so its rate scales with $\theta_\ast^2$: halving $\theta_\ast$ cuts the rate four times.',
                  [(r'$\theta_\ast$ lower than for the reference', 'the impurity (or what it turns into) takes the sites'),
                   (r'$\theta_\ast$ going up with T', r'adsorbates (\adsorbate{} or impurity) desorb; more room for \probe{}'),
                   (r'$\theta_\ast$ near 0', r'a nearly full surface: the rate is limited by finding free sites')]),
             ]},
            {'key': 'rates_summary_vs_P', 'caption': 'Main steady states along the pressure series: one row of panels per series temperature, one marker per probe-gas pressure.',
             'panels': [
                 ('(a) dissociation rate',
                  r'The forward TOF per case; one marker per probe-gas pressure.',
                  r'The impurities keep their ratio to the probe gas, so all partial pressures rise together.',
                  [('the TOF going up with pressure', r'a positive order in \probe{}: more \probe{} reaches the empty sites'),
                   ('the TOF flat with pressure', r'the surface is saturated with \adsorbate{}; more \probe{} does not help'),
                   ('the TOF coming down with pressure', r'the impurity gains more sites than \probe{} as everything is compressed')]),
                 ('(b) slowdown', r'As in the temperature series, at each pressure.', r"The dashed line is ``no effect''.",
                  [('the slowdown going up with pressure', 'the impurity competes better for sites at high pressure'),
                   ('the slowdown coming down with pressure', r'\probe{} wins the sites back at high pressure')]),
                 ('(c) empty sites', r'The fraction of empty sites at each pressure.', r'It falls as the surface fills.',
                  [(r'$\theta_\ast$ coming down with pressure', r'adsorbates (\adsorbate{} or impurity) fill the surface')]),
             ]},
        ],
    },

    'composition': {
        'title': 'What covers the surface',
        'goal': r"""Which species hold the sites at each steady state, and what does each impurity
displace? This explains the rates of the previous section: dissociation needs empty sites, so
whatever holds them sets the rate.""",
        'columns': [
            ('dominant species', r'the most abundant surface species (not counting empty sites) and its coverage; then the second one'),
            ('species, coverage', r'every surface species above the listing threshold, with its readable name (formula, $\ast$ per site, mechanism number) and fraction of sites'),
        ],
        'figures': [
            {'key': 'surface_composition', 'caption': 'What covers the surface at steady state.',
             'panels': [
                 ('each panel (one per point)',
                  r"Stacked bars, one per case (x axis). Each colored segment is the fraction of sites held by one species; pale gray: empty sites; dark gray ``other'': every species below the plotting threshold together.",
                  r'Each bar adds up to 1. Compare a case with the reference bar in the same panel to see what the impurity displaced.',
                  [('the tallest segment of a bar', 'the species that holds the surface in that case'),
                   (r'the \adsorbate{} segment shrinking in an impurity bar', 'the impurity (or its products) takes the sites from hydrogen'),
                   ('a segment growing from panel to panel', 'that adsorbate gains sites as the temperature (or pressure) changes'),
                   ('pale gray growing from panel to panel', 'adsorbates desorb and the surface empties'),
                   ('a large dark-gray segment', r'many minor species; they are listed in \texttt{steady\_state\_surface\_species.csv}')]),
             ]},
        ],
    },

    'diagnostics': {
        'title': 'Is the mechanism believable?',
        'goal': r"""RMG estimates binding energies; an adsorbate that binds too strongly is predicted to
cover too much of the surface and to poison too much. This section compares the mechanism's
adsorption enthalpies with literature ranges, switches off a suspect species to see how much
the conclusions depend on it, and lists the reactions that make and consume chosen species.""",
        'columns': [
            (r'$\Delta H_{298}$', r'adsorption enthalpy at 298 K from the mechanism\textquotesingle s own thermochemistry (kJ/mol); more negative = stronger binding'),
            ('literature', r'the range entered in settings.py (kJ/mol) and its source note'),
            ('verdict', 'OK (inside the range), TOO STRONG (below it), TOO WEAK (above it)'),
            (r'$\Delta G$', 'adsorption free energy at the stated temperature (kJ/mol)'),
            ('rate ratio', r'dissociation TOF with the species switched off / with it (1 = the species does not matter)'),
            ('slowdown normal / switched off', r'$1/\Pi$ with the full mechanism and without the species'),
            ('rate', r'the rate at which a reaction makes (+) or consumes ($-$) the species, per site (s$^{-1}$)'),
        ],
        'figures': [
            {'key': 'binding_check', 'caption': 'Adsorption enthalpies: mechanism against literature.',
             'panels': [
                 ('the plot',
                  r'One row per adsorption reaction: the literature range (gray bar) and the mechanism\textquotesingle s adsorption enthalpy at 298 K (marker).',
                  r'Left is stronger binding (more negative $\Delta H$).',
                  [('a marker inside its bar', 'the mechanism agrees with the literature for that adsorption'),
                   (r'a red marker left of its bar ($\blacktriangleleft$)', r'binds too strongly: that adsorbate\textquotesingle s coverage, and the poisoning it causes, are over-predicted'),
                   (r'a blue marker right of its bar ($\blacktriangleright$)', 'binds too weakly: its coverage is under-predicted'),
                   (r'a marker $\delta$ kJ/mol outside the range', r'$\theta/\theta_\ast$ of that adsorbate is off by about $e^{\delta/RT}$: 10 kJ/mol at 650 K is a factor of 6')]),
             ]},
        ],
    },

    'convergence': {
        'title': 'Did the solver reach a true steady state?',
        'goal': r"""The reactor method integrates the surface in time from a clean start and accepts the
state when the residual is small and the rate and coverages stopped changing. This section
follows chosen conditions to the end of the integration, to check that the accepted state is
final and not a long-lived intermediate.""",
        'columns': [
            ('accepted at', r'surface time (s) at which the reactor method accepted the state (-- = never)'),
            ('residual at acceptance / at the end', r'normalized residual $\max_k |d\theta_k/dt|$ / allowed; below 1 passes'),
            ('largest after acceptance', r'the largest residual after the accepted time; above 1 means something was still changing'),
            ('TOF change', r'relative change of the dissociation TOF from acceptance to the end of the integration'),
        ],
        'figures': [
            {'key': 'convergence', 'caption': 'Reactor-method convergence from a clean surface, one row per traced condition.',
             'panels': [
                 ('(a) coverages',
                  r'The largest coverages against surface time (log scale) from a clean surface; dashed gray: empty sites. The dotted vertical line is where the reactor method accepted the state.',
                  'Read from left to right: the surface fills, then settles.',
                  [('every coverage flat after the dotted line', 'a true steady state'),
                   ('a coverage still moving late (around $10^6$--$10^8$ s)', 'a slow conversion: the accepted state was long-lived but not final')]),
                 ('(b) forward TOF', r'The dissociation rate against surface time.', 'Log scale.',
                  [('flat after the dotted line', 'the rate quoted in the tables is final'),
                   ('drifting after the dotted line', 'treat the rate at that condition as provisional')]),
                 ('(c) normalized residual',
                  r'$\max_k |d\theta_k/dt|$ divided by the allowed value; the dashed line is the limit (1).',
                  'Below the dashed line the state passes the residual test.',
                  [('falling below 1 and staying there', 'converged'),
                   ('climbing back above 1 after acceptance', 'something slow is still changing; extend reactor\\_t\\_max or report a slow drift'),
                   ("``never accepted'' written in the panel", 'no plateau within reactor\\_t\\_max: no steady state at that condition')]),
             ]},
        ],
    },

    'curves': {
        'title': 'Poisoning curves: how much impurity is tolerated?',
        'goal': r"""From what impurity level does each impurity slow \probe{} dissociation, and where is
the rate halved? At every point of both series the probe-gas pressure is kept and each
impurity alone is added over a range of levels $r = p_\mathrm{imp}/p(\probe{})$.""",
        'columns': [
            (r'$r$ at half rate, ppm', r'the impurity level that halves the dissociation rate ($\Pi = 0.5$), as a ratio and in ppm of the probe gas'),
            ('note', r"``already halved at the lowest level'': the true value is below it (an upper bound); ``never halved'': above the highest level"),
            (r'$\Pi$ at lowest / highest level', 'the poisoning factor at the two ends of the range of levels'),
            ('two-state levels', 'how many levels have two steady states (the start matters there)'),
        ],
        'figures': [
            {'key': 'poisoning_curves', 'caption': r'Poisoning curves: dissociation rate relative to the probe gas alone (top) and impurity-held sites (bottom) against the impurity level.',
             'panels': [
                 ('top row: poisoning factor',
                  r'$\Pi$ against the impurity level $r$, one line per impurity, one column per point.',
                  r'Log--log. The dashed line is $\Pi = 0.5$; a diamond marks where a line crosses it (the tolerated level of the table).',
                  [('a line staying near 1', 'that impurity does not poison in the range tested'),
                   (r'a line coming down as $r$ rises', 'more impurity, slower dissociation; the steeper the drop, the more threshold-like the poison'),
                   ('a drop by orders of magnitude within one decade', 'a threshold poison: harmless below that level, the surface is lost above it'),
                   ('the crossing moving to the right from one column to the next', 'more impurity is tolerated at that point (e.g. at higher temperature)'),
                   ('a ring, a dashed branch and a shaded band', 'two steady states: solid = a fresh (clean) catalyst, dashed = one that was covered first; inside the band the history decides')]),
                 ('bottom row: impurity-held sites',
                  r'$\theta_\mathrm{imp}$, the sites held by the impurity and its products, against $r$.',
                  r'Log--log; values below $10^{-12}$ are drawn at $10^{-12}$.',
                  [(r'a line rising with slope 1', 'dilute adsorption: the impurity coverage is proportional to its pressure'),
                   ('a line levelling off near 1', 'the surface is saturated with the impurity'),
                   ('the top line falling where the bottom line rises', 'the impurity slows the rate by taking sites')]),
             ]},
            {'key': 'site_blocking', 'caption': 'Site blocking: does the poisoning come from lost sites? One row per impurity at one point.',
             'panels': [
                 ('(a) sites', r'Empty sites (gray) and impurity-held sites (color) against the impurity level.', 'Log scale.',
                  [('empty sites falling as impurity-held sites rise', 'the impurity takes the sites hydrogen needs')]),
                 ('(b) site-count factor', r'$f = \sum_k \sigma_k C_k/\Gamma$: the number of sites Cantera implies per nominal site.', 'The dashed line is 1.',
                  [('$f = 1$', 'one-site adsorbates only'),
                   ('$f$ rising above 1', 'multi-site adsorbates inflate the site count; the fixed-site-count estimate differs')]),
                 ('(c) poisoning factor', r'$\Pi$ from the rates (dots), the full site-blocking prediction $(\theta_\ast/\theta_{\ast,\mathrm{ref}})^2(f/f_\mathrm{ref})^2$ (line) and the coverage-only one (dashed).', 'Log scale.',
                  [('dots on the line', 'the impurity acts only by taking sites'),
                   ('the dashed line below the solid one', 'the coverage-only estimate exaggerates the poisoning because of $f$'),
                   ('dots off the line', 'something besides site blocking: coverage-dependent kinetics, or a probe step that is not elementary')]),
             ]},
        ],
    },

    'sweep': {
        'title': 'Temperature sweep',
        'goal': r"""How do the dissociation rate and the poisoning change with temperature, on a finer
temperature grid than the main series, at the reference pressure? Where does the rate recover
to half of the clean value (the release temperature)?""",
        'columns': [
            ('release T', 'the lowest temperature at which an impurity case is back to half the clean rate ($\\Pi \\geq 0.5$)'),
            ('note', 'rate never below half (not poisoning), or still poisoned at the highest temperature'),
            (r'$E_\mathrm{app}$', r'local apparent activation energy $-R\,d\ln\mathrm{TOF}/d(1/T)$ (kJ/mol)'),
        ],
        'figures': [
            {'key': 'temperature_sweep', 'caption': 'Temperature sweep at the reference pressure.',
             'panels': [
                 ('(a) Arrhenius plot', r'Dissociation TOF (log) against $1000/T$, one line per case; hot is on the left.', 'The slope is $-E_\\mathrm{app}/R$.',
                  [('a straight line', 'a constant apparent activation energy'),
                   ('a bending line', 'the coverages change with temperature (the barrier the surface shows changes)'),
                   ('an impurity line approaching the reference on the left', 'poisoning lifts at high temperature'),
                   ('a gap in a line', 'no steady state at that temperature')]),
                 ('(b) apparent activation energy', r'Local $E_\mathrm{app}$ against $T$.', 'Positive: the rate rises with T.',
                  [('$E_\\mathrm{app}$ near 0', 'no barrier: the rate is set by how often the gas hits empty sites'),
                   ('an impurity case far above the reference', 'the impurity desorbs as T rises, freeing sites (expect a release temperature)'),
                   ('$E_\\mathrm{app}$ below 0', 'the rate falls as T rises (fewer empty sites, or a sticking coefficient falling with T)')]),
                 ('(c) poisoning factor', r'$\Pi$ against $T$ for each impurity case; the dashed line is $\Pi = 0.5$.', 'Log scale; a diamond marks the release temperature.',
                  [('a curve going up with T', 'poisoning lifts with temperature'),
                   ('a curve coming down with T', 'poisoning gets worse with temperature, e.g. a decomposition product binds more'),
                   ('a flat curve at 1', 'that impurity does not poison at any temperature')]),
             ]},
            {'key': 'sweep_coverages', 'caption': 'What covers the surface against temperature, one panel per case.',
             'panels': [
                 ('each panel (one per case)', 'Stacked bars of the coverages at each temperature of the sweep.', 'Follow one color across the bars.',
                  [('a color giving way to another as T rises', 'one adsorbate replaces another; this explains the bends in the Arrhenius plot'),
                   ('pale gray (empty sites) growing with T', 'the surface empties as adsorbates desorb'),
                   ('a missing bar', 'no steady state at that temperature')]),
             ]},
        ],
    },

    'arrhenius': {
        'title': 'Hydrogen rates in Arrhenius form',
        'goal': r"""The rates of hydrogen on the metal in the form $k = A\,\exp(-E_a/RT)$, so that metals
and impurities can be compared by two numbers. Two levels: the \emph{intrinsic} constants of
the \probe{} step itself (dissociation on an empty surface, recombination from a surface full
of \adsorbate{}, and their ratio, the equilibrium constant), which belong to the metal and do
not depend on impurities; and the \emph{observed} rates on the real, crowded surface of each
case, where the impurities act. An impurity\textquotesingle s effect is then a change of the
activation energy, $\Delta E_a$, and of the prefactor, $A/A_\mathrm{ref}$.""",
        'columns': [
            ('expression', r'$A\,\exp(-E_a/RT)$ with $E_a$ in kJ/mol and $R = 8.314\times10^{-3}$ kJ/(mol K); $A$ in s$^{-1}$ per site for rates, dimensionless for $K$ and for forward/reverse'),
            ('conditions', 'the gas and surface each fit belongs to; the temperature range is the range of the fit'),
            (r'$R^2$', r'quality of the straight-line fit of $\ln y$ against $1/T$; -- when $y$ is constant (e.g.\ forward/reverse = 1)'),
            (r'$\Delta E_a$, $A/A_\mathrm{ref}$', 'activation energy minus, and prefactor relative to, the reference case (probe gas alone) for the same quantity'),
            ('curved', r'$R^2$ below the threshold: one $E_a$ does not describe the whole range; use the local values'),
            (r'$k_\mathrm{diss}$', r'dissociation rate per site on an empty surface at the reference pressure (s$^{-1}$)'),
            (r'$k_\mathrm{rec}$', r'recombination rate per site on a surface full of \adsorbate{} (s$^{-1}$)'),
            ('$K$', r'equilibrium constant $\theta_\mathrm{H}^2/(\theta_\ast^2\,p/P^\circ)$, $P^\circ = 1$ bar; $\Delta H^\circ$, $\Delta S^\circ$ from the mechanism thermochemistry'),
        ],
        'figures': [
            {'key': 'arrhenius', 'caption': r'Hydrogen rates in Arrhenius form.',
             'panels': [
                 ('(a) intrinsic rates',
                  r'$k_\mathrm{diss}$ (empty surface) and $k_\mathrm{rec}$ (surface full of \adsorbate{}) per site against $1000/T$, with their fits; $E_a$ in the legend.',
                  'Hot on the left. A steeper line is a higher activation energy.',
                  [(r'$k_\mathrm{diss}$ nearly flat', r'no barrier to dissociation: the rate is the impingement flux ($\propto p/\sqrt{T}$) times the sticking coefficient'),
                   (r'$k_\mathrm{rec}$ falling steeply to the right', r'recombination is activated; its $E_a$ is close to the H$_2$ desorption energy'),
                   ('a higher line', 'a faster step at that temperature')]),
                 ('(b) equilibrium constant',
                  r'$K$ against $1000/T$ (van \textquotesingle t Hoff plot); $\Delta H^\circ$ and $\Delta S^\circ$ from the fit are written in the panel.',
                  r'Slope $= -\Delta H^\circ/R$.',
                  [(r'$K$ going up to the right (falling as T rises)', r'adsorption is exothermic ($\Delta H^\circ < 0$): hydrogen leaves the surface as T rises'),
                   (r'$K \gg 1$', r'the surface prefers \adsorbate{} (covered at that temperature)'),
                   (r'$K \ll 1$', r'\adsorbate{} is unstable; the surface stays mostly empty')]),
                 ('(c) observed forward TOF',
                  r'The dissociation TOF of each case on its steady-state surface against $1000/T$, with its fit line.',
                  'Each case at its own gas (see the conditions column).',
                  [('points on the line', 'the Arrhenius form holds over the range'),
                   ('points bending away from the line', "flagged ``curved'': use the local activation energies of the table"),
                   ('an impurity case below and steeper than the reference', 'the impurity blocks sites at low T and leaves at high T')]),
                 ('(d) activation energies',
                  r'$E_a$ of the dissociation (forward) and recombination (reverse) TOF per case; the dashed line is the intrinsic dissociation $E_a$.',
                  'Bars are fits over the whole range.',
                  [('observed $E_a$ above the intrinsic line', r'dissociation needs free sites, so the observed barrier includes the cost of freeing them (desorbing \adsorbate{} or the impurity)'),
                   (r'$\Delta E_a > 0$ for an impurity', 'the impurity adds an energy cost: it has to leave before hydrogen can dissociate'),
                   ('forward and reverse bars of equal height', r'the \probe{} step is at equilibrium: dissociation and recombination change together')]),
             ]},
        ],
    },

    'pressure_sweep': {
        'title': 'Pressure sweep',
        'goal': r"""How does the dissociation rate change with the probe-gas pressure at fixed
temperature, with every impurity keeping its ratio to the probe gas (the same purity at every
pressure)? Where does the poisoning factor cross one half (the half-rate pressure)?""",
        'columns': [
            ('half-rate pressure', r'the probe-gas pressure (bar) where $\Pi$ crosses 0.5'),
            ('direction', r'whether the rate recovers above half or falls below half as the pressure goes up'),
            ('slope', r'local $d\ln\mathrm{TOF}/d\ln p$: the order in the probe gas for the reference; the order in the total pressure at fixed purity for impurity cases'),
        ],
        'figures': [
            {'key': 'pressure_sweep', 'caption': 'Pressure sweep at fixed purity, one row per temperature.',
             'panels': [
                 ('(a) dissociation rate', r'TOF against $p(\probe{})$, log--log, one line per case.', 'A slope of 1 is first order.',
                  [('slope 1', r'few adsorbates: rate $\propto p$'),
                   ('a flat line', r'the surface is saturated with \adsorbate{}'),
                   ('an impurity line moving away from the reference as p rises', 'that impurity grows more harmful with pressure')]),
                 ('(b) local slope', r'$d\ln\mathrm{TOF}/d\ln p$ against $p$.', 'The order at each pressure.',
                  [('the reference near 1', 'few adsorbates on the surface'),
                   ('the reference near 0', r'the surface is saturated with \adsorbate{}'),
                   ('an impurity case near $-1$', 'a strong poison that scales with pressure: purity, not pressure, sets the poisoning')]),
                 ('(c) poisoning factor', r'$\Pi$ against $p$; dashed: $\Pi = 0.5$; a diamond marks the half-rate pressure.', 'Log scale.',
                  [('a curve going down with p', 'more pressure, more poisoning at the same purity'),
                   ('a curve going up with p', r'\probe{} wins the sites back at high pressure')]),
             ]},
            {'key': 'pressure_sweep_coverages', 'caption': 'What covers the surface against pressure, one panel per case.',
             'panels': [
                 ('each panel (one per case)', 'Stacked bars of the coverages at each pressure of the sweep.', 'Follow one color across the bars.',
                  [('a color growing with pressure', 'that adsorbate gains sites as everything is compressed'),
                   ('pale gray shrinking with pressure', 'the surface fills up')]),
             ]},
        ],
    },

    'orders': {
        'title': 'Reaction orders',
        'goal': r"""How sensitive is the dissociation rate to each partial pressure? The order
$n_j = d\ln\mathrm{TOF}/d\ln p_j$ is found by changing one partial pressure by $\pm$10\,\% with the
others fixed, at every point of both series.""",
        'columns': [
            (r'order in \probe{}', r'$d\ln\mathrm{TOF}/d\ln p(\probe{})$ with the impurity fixed'),
            ('order in the impurity', r'$d\ln\mathrm{TOF}/d\ln p_\mathrm{imp}$ with the probe gas fixed (-- for the reference)'),
            ('starts agree', 'all perturbed steady states had a single steady state'),
        ],
        'figures': [
            {'key': 'reaction_orders', 'caption': 'Reaction orders of the dissociation rate, one panel per point.',
             'panels': [
                 ('each panel (one per point)', r'Per case: circle = order in \probe{}, square = order in the impurity. Values are written where $|n| \geq 0.1$ and, for an impurity case, where they differ from the reference.', 'The gray line is order 0.',
                  [(r'order in \probe{} near 1', 'few adsorbates: more gas, proportionally faster'),
                   (r'order in \probe{} near 0', r'the surface is saturated with \adsorbate{}'),
                   (r'order in \probe{} rising with an impurity', r'the impurity removes \adsorbate{}, so the surface is less saturated'),
                   ('order in the impurity near 0', 'the impurity barely matters at this level'),
                   ('order in the impurity near $-2$', r'a strongly bound poison controls the empty sites ($\theta_\ast \propto 1/p_\mathrm{imp}$, rate $\propto \theta_\ast^2$)'),
                   ('a positive order in the impurity', 'the impurity speeds dissociation, e.g. by removing a site blocker')]),
             ]},
        ],
    },

    'rate_control': {
        'title': 'Rate control: which species and steps set the rate',
        'goal': r"""Which surface species hold the sites the reaction needs, and which elementary steps
limit the rate? The thermodynamic rate control $X_\mathrm{TRC,n} = \partial\ln r/\partial(-G_n/RT)$
measures how the dissociation rate responds when a species is made more stable; the degree of
rate control $X_\mathrm{DRC,i} = \partial\ln r/\partial\ln k_i$ how it responds to speeding up one
reaction.""",
        'columns': [
            ('top TRC species', r'the species with the most negative $X_\mathrm{TRC}$ (the main site holder) and its value'),
            ('TRC sum', r'sum of all $X_\mathrm{TRC}$; about $-2(1-\theta_\ast)$ for site blocking of a two-site step'),
            ('DRC of the probe step', r'$X_\mathrm{DRC}$ of the forward rate for the probe step itself (1: it alone limits the rate)'),
            ('largest other DRC step', r'the other reactions with the largest $|X_\mathrm{DRC}|$ of the forward rate (\#index: equation), and their values'),
        ],
        'figures': [
            {'key': 'rate_control', 'caption': r'Thermodynamic rate control of the dissociation rate (one panel per case).',
             'panels': [
                 ('each panel (one per case)', r'Bars of $X_\mathrm{TRC}$ per surface species with $|X| \geq 0.01$; the sum is under the axis.', 'Negative: making that species more stable slows dissociation.',
                  [('a large negative bar', 'that species holds sites the reaction needs'),
                   ('one bar near $-2$', 'that species holds nearly all sites: its binding energy controls the rate (check it in the binding check)'),
                   ('a positive bar', 'stabilizing that species speeds the rate'),
                   (r'bars summing to about $-2(1-\theta_\ast)$', 'the ordinary site-blocking picture of a two-site step')]),
             ]},
            {'key': 'rate_control_DRC', 'caption': r'Degree of rate control (cases where it is not trivial).',
             'panels': [
                 ('each panel (one per case)', r'Bars of $X_\mathrm{DRC}$ for the largest reactions (\#index: equation); blue: of the forward rate, orange: of the net \probe{} consumption.', 'Net values add up to 1.',
                  [('1 on the probe step and 0 elsewhere', 'the probe step is the slowest; the impurity acts only through the coverages'),
                   ('a large value on an impurity step', 'that step (e.g. the impurity\\textquotesingle s adsorption or removal) controls the rate'),
                   ('a negative value', 'speeding up that step slows dissociation (it produces a site blocker)')]),
             ]},
        ],
    },

    'paths': {
        'title': 'Reaction paths of hydrogen',
        'goal': r"""Where do the hydrogen atoms go at steady state? The diagrams show the flux of H atoms
between species, for each case at one point here, and at every other point in the appendix.""",
        'columns': [],
        'figures': [
            {'key': 'reaction_paths', 'caption': 'Hydrogen flux between species at steady state.',
             'panels': [
                 ('each diagram',
                  r'Arrows give the direction of flow; the label is the flux relative to the largest; the width is its logarithm. A node\textquotesingle s outline thickness is the logarithm of its amount (coverage or mole fraction; gray: below $10^{-6}$); nodes show the structures.',
                  r"``Scale'' is the largest flux (kmol m$^{-2}$ s$^{-1}$; divide by the site density for s$^{-1}$). Both directions are drawn (one-way flow); fluxes below $10^{-8}$ of the largest are hidden.",
                  [(r'a thick arrow from \probe{} to \adsorbate{}', 'dissociation is the main channel of hydrogen'),
                   ('two equal arrows in opposite directions', 'that pair is at equilibrium; its net flow is tiny'),
                   ('hydrogen flowing into an impurity-derived species', 'the impurity consumes hydrogen (e.g. hydrogenation of O or C)'),
                   ('a new arrow in an impurity case', 'the impurity opens a path that the probe gas alone does not have')]),
             ]},
        ],
    },

    'transients': {
        'title': 'Poisoning over time',
        'goal': r"""How fast does a working catalyst get poisoned? Each run starts from the steady state in
the probe gas alone; at $t = 0$ the impurity is added and the dissociation rate is followed as
the impurity builds up on the surface.""",
        'columns': [
            ('half time, tenth time', 'time (s) until the rate is half, and a tenth, of its value at $t = 0$ (-- = never)'),
            (r'deepest, final', 'the lowest and the final rate relative to $t = 0$'),
            ('settle time', 'time (s) after which the rate stays within 10\\,\\% of its final value'),
        ],
        'figures': [
            {'key': 'poisoning_over_time', 'caption': 'Poisoning over time after the impurity arrives.',
             'panels': [
                 ('top row: rate', r'TOF relative to $t = 0$ against the time after the impurity arrives, one column per point and impurity; light to dark: more impurity.', 'Log--log; dashed: half rate.',
                  [('a line crossing the dashed line', 'the time to half rate (the half time of the table)'),
                   ('darker lines falling earlier', 'more impurity, faster poisoning'),
                   ('a dip below the final level, then partial recovery', 'a transient worse than the steady state: the worst moment is the deepest value, not the final one'),
                   ('a line staying at 1', 'no poisoning at that level')]),
                 ('bottom row: impurity-held sites', r'$\theta_\mathrm{imp}$ against time.', 'Log--log.',
                  [('a rising line', 'the impurity builds up on the surface'),
                   ('the rise matching the drop above', 'the rate falls because the sites fill with impurity')]),
             ]},
        ],
    },

    'flux': {
        'title': 'Flux diagrams of the RMG reactors',
        'goal': r"""How do species flow in the RMG input file\textquotesingle s own batch reactors over time? These
are closed reactors whose gas changes as it reacts, unlike the fixed-gas steady states of the
other sections; one frame per reactor is shown, and the full video is in its folder.""",
        'columns': [],
        'figures': [
            {'key': 'flux', 'caption': 'RMG reactor flux diagram (one frame).',
             'panels': [
                 ('each frame', 'Species-to-species fluxes at one time in one reactor; edge width is the flux (log scale).', 'Compare node borders only within the gas or within the surface.',
                  [('a thick edge', 'a dominant conversion at that time'),
                   ('edges appearing in later frames', 'products accumulate and open new paths')]),
             ]},
        ],
    },
}


# --- Format 2: per area (per_area.py; docs/rates_units_equilibrium.md) -----------------------
CONVENTIONS = CONVENTIONS.replace(r"""\end{itemize}""", r"""\item \textbf{Two formats.} \emph{Format 1} gives rates per site (TOF, s$^{-1}$). \emph{Format 2}
  gives the same results per area of surface, the units of the hydrogen-surface literature:
  fluxes $J$ = TOF $\times\,\Gamma$ (mol \probe{} m$^{-2}$ s$^{-1}$), the effective adsorption
  coefficient $k_{d,\mathrm{eff}} = J_\mathrm{diss}/p$ (mol \probe{} m$^{-2}$ s$^{-1}$ Pa$^{-1}$), the
  surface recombination coefficient $k_r$ of $J_\mathrm{rec} = k_r(\theta_\mathrm{H}\Gamma)^2$
  (m$^2$ mol$^{-1}$ s$^{-1}$), the sticking probability $\gamma_\mathrm{eff} = J_\mathrm{diss}/Z$
  ($Z = p/\sqrt{2\pi MRT}$, the impingement flux) and the hydrogen on the surface
  $C_\mathrm{H} = \theta_\mathrm{H}\Gamma$ (mol H m$^{-2}$). Figures in Format 2 have the suffix
  \texttt{\_per\_area}; the site density $\Gamma$ they use is on their axis. Converting does not
  change any activation energy, only prefactors. $\Gamma$ belongs to an ideal facet, so ``per
  m$^2$'' is per atomic surface area: literature values per geometric area of a real, rough
  sample can be higher by a roughness factor (typically 1--10). Fluxes count \probe{}
  molecules, surface concentrations H atoms. Recombination is the surface step only; there is
  no subsurface or bulk hydrogen in the mechanism.
\end{itemize}""")

_F2_COLUMNS = [
    (r'$J_\mathrm{diss}$, $J_\mathrm{rec}$, $J_\mathrm{net}$', r'dissociation, recombination and net flux = TOF $\times\,\Gamma$ (mol \probe{} m$^{-2}$ s$^{-1}$)'),
    (r'$J_\mathrm{ex}$', r'exchange flux $\min(J_\mathrm{diss}, J_\mathrm{rec})$: the part that cycles between gas and surface; $J_0$ at equilibrium'),
    (r'$k_{d,\mathrm{eff}}$', r'effective adsorption coefficient $J_\mathrm{diss}/p$ (mol \probe{} m$^{-2}$ s$^{-1}$ Pa$^{-1}$) $= k_d\,\theta_v^2$'),
    (r'$\gamma_\mathrm{eff}$', r'sticking probability per collision, $J_\mathrm{diss}/Z$ (0 to 1)'),
    (r'$C_\mathrm{H}$', r'hydrogen on the surface, $\theta_\mathrm{H}\Gamma$ (mol H m$^{-2}$)'),
]
_SAME = ('(b), (c)', 'As in the per-site figure.', 'They have no units, so they are identical.', [])

SECTIONS['steady_state']['columns'] += _F2_COLUMNS
SECTIONS['steady_state']['figures'] += [
    {'key': 'rates_summary_vs_T_per_area', 'caption': 'Main steady states along the temperature series, per area (Format 2).',
     'panels': [
         ('(a) adsorption coefficient',
          r'The effective adsorption coefficient $k_{d,\mathrm{eff}} = J_\mathrm{diss}/p$ of each case, in the unit of literature adsorption coefficients; the axis gives the site density $\Gamma$ used.',
          r'It is the per-site TOF times $\Gamma/p$: the same pattern in per-area units. On a clean surface it would be $k_d = \gamma/\sqrt{2\pi MRT}$; on the steady-state surface $k_{d,\mathrm{eff}} = k_d\,\theta_v^2$.',
          [('a case below the reference', r'the impurity lowers the adsorption coefficient by taking empty sites ($\theta_v^2$)'),
           ('a value compared with a literature $k_d$', 'ours is per atomic area of an ideal facet; literature values per geometric area of a rough sample can be higher by a roughness factor of 1--10')]),
         _SAME]},
    {'key': 'rates_summary_vs_P_per_area', 'caption': 'Main steady states along the pressure series, per area (Format 2).',
     'panels': [
         ('(a) adsorption coefficient', r'$k_{d,\mathrm{eff}} = J_\mathrm{diss}/p$ of each case at each probe-gas pressure.',
          r'A coefficient: the pressure is divided out.',
          [('$k_{d,\\mathrm{eff}}$ coming down as the pressure rises', r'the surface fills ($\theta_v$ falls); first-order adsorption with a constant $k_d$, as the literature assumes, would give a flat line'),
           ('$k_{d,\\mathrm{eff}}$ flat with pressure', r'few adsorbates: $J = k_d\,p$ holds')]),
         _SAME]},
]

SECTIONS['sweep']['figures'] += [
    {'key': 'temperature_sweep_per_area', 'caption': 'Temperature sweep at the reference pressure, per area (Format 2).',
     'panels': [
         ('(a) Arrhenius plot per area', r'$k_{d,\mathrm{eff}}$ (log) against $1000/T$, one line per case.',
          r'The same slope as the per-site Arrhenius plot: converting does not change $E_\mathrm{app}$.',
          [('a straight line', 'one activation energy, the same as per site'),
           ('lines steeper than the clean-surface $k_d$ (Arrhenius section)', 'the step waits for free sites; the observed barrier includes freeing them')]),
         _SAME]},
]

SECTIONS['arrhenius']['columns'] += _F2_COLUMNS + [
    (r'$k_d$, $k_r$ (literature-style table)', r'the clean-surface adsorption coefficient (mol \probe{} m$^{-2}$ s$^{-1}$ Pa$^{-1}$; for a steady state: $k_{d,\mathrm{eff}}$) and the surface recombination coefficient (m$^2$ mol$^{-1}$ s$^{-1}$), each as $A\,\exp(-E_a/RT)$; $k_r$ is a property of the surface step, the same in every case'),
    (r'$\gamma_\mathrm{clean}$, $\theta_\mathrm{H,eq}$, $C_\mathrm{H,eq}$, $J_0$', r'on a surface holding only hydrogen at the reference pressure: the sticking probability (the mechanism\textquotesingle s sticking coefficient), the equilibrium coverage and amount of H, and the equilibrium exchange flux'),
    ('check $K$', r'$K$ from the rates divided by $K$ from the thermochemistry: 1 for a thermodynamically consistent mechanism'),
]
SECTIONS['arrhenius']['figures'] += [
    {'key': 'arrhenius_per_area', 'caption': 'Hydrogen rates in Arrhenius form, per area (Format 2).',
     'panels': [
         ('(a) adsorption coefficient, clean surface', r'$k_d$ (mol \probe{} m$^{-2}$ s$^{-1}$ Pa$^{-1}$) against $1000/T$ with its fit.',
          r'$k_d = \gamma/\sqrt{2\pi MRT}$: the sticking coefficient times a kinetic-theory factor.',
          [('nearly flat', r'no barrier to dissociation; a constant $\gamma$ gives $E_a \approx -RT/2$ (from $1/\sqrt{T}$)'),
           ('rising to the left (higher at high T)', 'activated sticking: the barrier in the sticking coefficient')]),
         ('(b) surface recombination coefficient', r'$k_r$ (m$^2$ mol$^{-1}$ s$^{-1}$) of $J_\mathrm{rec} = k_r(\theta_\mathrm{H}\Gamma)^2$ against $1000/T$.',
          r'A property of the surface step: the same in every case.',
          [('falling steeply to the right', 'recombination is activated; its $E_a$ is close to the H$_2$ desorption energy'),
           ('$E_a$(recombination) $-$ $E_a$(dissociation)', r'equals $-\Delta H^\circ$ of adsorption (thermodynamic consistency)')]),
         ('(c) effective adsorption coefficient', r'$k_{d,\mathrm{eff}}$ of each case on its steady-state surface, with fits.',
          r'These rows are the literature-style table: a per-area $k_d$ for each gas.',
          [('a case far below the clean-surface $k_d$', r'its surface has few empty sites ($k_{d,\mathrm{eff}} = k_d\,\theta_v^2$)'),
           ('an impurity case steeper than the reference', 'the impurity blocks sites at low T and leaves at high T')]),
         ('(d) activation energies', 'As in the per-site figure.', 'Converting does not change $E_a$.', [])]},
    {'key': 'surface_exchange_vs_T', 'caption': 'Sticking probability, exchange flux and hydrogen on the surface against temperature (temperature sweep).',
     'panels': [
         ('(a) sticking probability', r'$\gamma_\mathrm{eff} = J_\mathrm{diss}/Z$: the fraction of \probe{} molecules hitting the surface that dissociate. Dashed: the mechanism\textquotesingle s sticking coefficient (a clean surface).',
          'Log scale, between 0 and 1.',
          [('a case below the dashed line', r'fewer empty sites: $\gamma_\mathrm{eff} = \gamma\,\theta_v^2$'),
           ('going up with T', 'the surface empties as adsorbates desorb'),
           ('an impurity case below the reference', r'the impurity takes the sites \probe{} needs; the gap is the impurity effect as a probability')]),
         ('(b) exchange flux', r'$J_\mathrm{ex}$ against $1000/T$. Dashed: $J_0$ of a surface holding only hydrogen at equilibrium.',
          r'The rate at which hydrogen cycles between gas and surface, what isotope exchange (\probe{} + D$_2$) measures.',
          [('the probe-only line on the dashed line', 'that surface is at the equilibrium coverage'),
           ('an impurity case below the reference', 'fewer H atoms and empty sites to exchange: isotope exchange would be slower'),
           ('a steep line', 'activated exchange')]),
         ('(c) hydrogen on the surface', r'$C_\mathrm{H} = \theta_\mathrm{H}\Gamma$ (mol H m$^{-2}$). Dashed: $C_\mathrm{H,eq}$ of a surface holding only hydrogen.',
          'Linear scale.',
          [('coming down with T', 'hydrogen desorbs (adsorption is exothermic)'),
           ('an impurity case below the reference', 'the impurity displaces hydrogen from the surface')]),
     ]},
]

SECTIONS['pressure_sweep']['columns'] += _F2_COLUMNS
SECTIONS['pressure_sweep']['figures'] += [
    {'key': 'pressure_sweep_per_area', 'caption': 'Pressure sweep at fixed purity, per area (Format 2), one row per temperature.',
     'panels': [
         ('(a) dissociation flux', r'$J_\mathrm{diss}$ (mol \probe{} m$^{-2}$ s$^{-1}$) against $p$, log--log.',
          'A flux, so the slope is the reaction order, as per site.',
          [('slope 1', r'first order: $J = k_d\,p$, as the literature assumes'),
           ('the line flattening', r'the surface saturates with \adsorbate{}; $J = k_d\,p$ no longer holds')]),
         _SAME]},
    {'key': 'surface_exchange_vs_P', 'caption': 'Sticking probability, exchange flux and hydrogen on the surface against pressure, one row per temperature.',
     'panels': [
         ('(a) sticking probability', r'$\gamma_\mathrm{eff}$ against $p$; dashed: the clean-surface sticking coefficient.', 'Log--log.',
          [('coming down with pressure', r'the surface fills, so fewer collisions find two empty sites')]),
         ('(b) exchange flux', r'$J_\mathrm{ex}$ against $p$; dashed: $J_0$ at the equilibrium coverage.', 'Log--log.',
          [('rising more slowly than $p$', 'saturation: more pressure, but no more empty sites')]),
         ('(c) hydrogen on the surface', r'$C_\mathrm{H}$ against $p$; dashed: $C_\mathrm{H,eq}$.', 'Log--log.',
          [('rising as $\\sqrt{p}$ at low pressure', 'dissociative adsorption: one \\probe{} gives two H'),
           (r'levelling off near $\Gamma$', 'the surface is full of hydrogen')]),
     ]},
]

SECTIONS['equilibrium'] = {
    'title': 'Equilibrium coverage of hydrogen',
    'goal': r"""How much hydrogen sits on the clean metal at equilibrium with the gas, at each
temperature and pressure, and is the \probe{} step at equilibrium in every steady state? The
equilibrium constant $K$ of the step alone, from the mechanism\textquotesingle s
thermochemistry, fixes $\theta_\mathrm{H}/\theta_v = (K\,p/P^\circ)^{1/n}$. On a surface holding
only hydrogen this is the dissociative Langmuir isotherm. In any steady state the same ratio
holds whenever the step is at equilibrium, whatever else covers the surface: an impurity lowers
$\theta_\mathrm{H}$ and $\theta_v$ together but cannot change their ratio.""",
    'columns': [
        ('$K$', r'dimensionless equilibrium constant $\theta_\mathrm{H}^2/(\theta_v^2\,p/P^\circ)$, $P^\circ$ = 1 bar, from the thermochemistry'),
        ('check $K$', r'$K$ from the rates divided by $K$ from the thermochemistry (1 = consistent)'),
        (r'$\theta_\mathrm{H,eq}$, $\theta_{v,\mathrm{eq}}$, $C_\mathrm{H,eq}$', r'equilibrium coverage and amount of hydrogen (mol H m$^{-2}$) of a surface holding only hydrogen'),
        ('$J_0$', r'the exchange flux at that equilibrium (mol \probe{} m$^{-2}$ s$^{-1}$)'),
        ('$T_{1/2}$', r'the temperature at which that surface is half covered at the given pressure ($K\,p/P^\circ = 1$)'),
        ('$q$', r'$(\theta_\mathrm{H}/\theta_v)/(K\,p/P^\circ)^{1/n}$ of a steady state: 1 = the \probe{} step is at equilibrium there'),
    ],
    'figures': [
        {'key': 'equilibrium_coverage', 'caption': 'Equilibrium coverage of a surface holding only hydrogen, and the equilibrium check of the steady states.',
         'panels': [
             ('(a) isotherms', r'$\theta_\mathrm{H,eq}$ against $p$, a line per temperature; open circles: the steady states of the probe gas alone.',
              'Log pressure axis.',
              [('rising with $p$, as $\\sqrt{p}$ at low pressure', 'dissociative adsorption'),
               ('levelling off near 1', 'the surface is saturated with hydrogen'),
               ('the circles on the lines', 'the probe-gas-only steady states are at the equilibrium coverage')]),
             ('(b) against temperature', r'$\theta_\mathrm{H,eq}$ against $T$, a line per pressure; diamonds: $T_{1/2}$; circles: steady states.',
              r'S-shaped: covered at low $T$, empty at high $T$.',
              [('coming down with T', r'adsorption is exothermic ($\Delta H^\circ < 0$): hydrogen leaves as T rises'),
               ('a higher pressure shifting the curve to the right', r'more pressure keeps hydrogen on to a higher temperature ($T_{1/2}$ rises)')]),
             ('(c) equilibrium check', r'$q$ for each case at each main point; the dashed line is 1.', 'Linear scale around 1.',
              [('$q = 1$', r'the \probe{} step is at equilibrium in that steady state; an impurity acts only through $\theta_v$'),
               ('$q$ below 1', 'fewer H atoms than equilibrium: other reactions consume H faster than dissociation replaces it, so the step is not equilibrated and limits the rate'),
               ('$q$ above 1', r'more H than equilibrium: other reactions produce H (e.g.\ an impurity decomposing), which leaves as \probe{}')]),
         ]},
        {'key': 'equilibrium_coverage_per_area', 'caption': 'Equilibrium coverage per area (Format 2).',
         'panels': [
             ('(a), (b)', r'As in the per-site figure, as $C_\mathrm{H,eq}$ (mol H m$^{-2}$); the axis gives $\Gamma$.',
              r'The maximum, $\Gamma$, is a full surface.', [('a value', 'the hydrogen held on each m$^2$ of the metal at equilibrium')]),
             ('(c)', 'As in the per-site figure.', 'Identical: $q$ has no units.', []),
         ]},
    ],
}
