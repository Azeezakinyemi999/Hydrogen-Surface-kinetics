# Hydrogen-Surface-kinetics

Steady-state and transient analysis of surface reaction mechanisms with [Cantera](https://cantera.org): how fast a probe reaction runs (for example H₂ dissociation on a metal), what covers the surface, how impurity gases poison it, and which reactions and species control the rate.

It works with any surface mechanism in Cantera's YAML format, usually exported by RMG. A new mechanism, metal or impurity set needs only a new `settings.py`, which can sit next to the run's RMG files (RMG's own input is `input.py`).

Temperature and pressure effects are kept apart: a **temperature series** (`temperatures` at `reference_pressure`) and a **pressure series** (`pressures` at each of `pressure_temperatures`), never a grid. Every figure comes as `<name>_vs_T` and, with a pressure series, `<name>_vs_P`.

## What it computes

| Analysis | Answers |
|---|---|
| `steady_state` | forward/net probe TOF, coverages and poisoning factor for every case (probe gas only, each impurity, all together) and temperature, each from two starting surfaces |
| `composition` | which species cover the surface |
| `diagnostics` | binding enthalpies vs literature; what-if runs with a species switched off; main sources and sinks of a species |
| `convergence` | how the reactor method reached (or failed to reach) each steady state: coverages, TOF and residual against time |
| `curves` | rate vs impurity level and the level that halves the rate; both branches where two steady states exist; the site-blocking split |
| `sweep` | Arrhenius plot, apparent activation energy, the temperature at which poisoning lifts, surface coverage against T |
| `arrhenius` | the hydrogen rates as $A\,\exp(-E_a/RT)$: intrinsic dissociation, recombination and equilibrium constants of the metal, and the observed rates per case, with each impurity's $\Delta E_a$ and $A/A_{\mathrm{ref}}$ (no solving) |
| `pressure_sweep` | rate against H₂ pressure at fixed purity, its local slope, the half-rate pressure, surface coverage against p |
| `equilibrium` | equilibrium coverage of hydrogen from the thermochemistry (isotherms, the half-coverage temperature) and whether the H₂ step is at equilibrium in every steady state (no solving) |
| `orders` | reaction orders in the probe gas and the impurity |
| `rate_control` | thermodynamic rate control (TRC) and degree of rate control (DRC), one figure per (T, p) point |
| `paths` | element-flux diagrams at a steady state |
| `transients` | how fast a working catalyst is poisoned after an impurity arrives |
| `flux` | flux-diagram videos of the RMG input's batch reactors (optional, needs RMG-Py) |

## Install

Into an environment that has Cantera ≥ 3.0, NumPy, pandas and matplotlib (here the `rmg4.0` conda environment, which also has RMG-Py for `flux`):

```bash
cd Hydrogen-Surface-kinetics
~/anaconda3/envs/rmg4.0/bin/python -m pip install --no-deps -e .
```

`--no-deps` leaves the environment's existing packages untouched; `-e` means edits to the code take effect without reinstalling. Graphviz (`dot`) is needed for `paths` and `flux`, and ffmpeg for `flux`.

## Run

Each run is a folder with a `settings.py`. `runs/ni_first_try/` is a complete example (H₂ on Ni(111) with H₂O, NH₃, CH₄, CO and CO₂, temperature series only); `runs/try_*/<Cu|Fe|Ni>/` have settings with both series next to their RMG files.

```bash
bash run.sh runs/ni_first_try check                 # settings, mechanism, both series, solves needed (seconds)
bash run.sh runs/ni_first_try                       # every analysis listed in settings.py
bash run.sh runs/ni_first_try --only curves,sweep   # a subset
bash run.sh runs/ni_first_try status                # what is cached and saved
bash run.sh runs/ni_first_try report                # rebuild the PDF report from the saved results
bash run.sh runs/try_2/Cu runs/try_2/Fe runs/try_2/Ni   # several runs, one after another
```

Give each run as its folder or its `settings.py` (globs such as `runs/try_*/Ni` work), then `check` or `status` if wanted, then options; the options apply to every run. A failed run does not stop the others, and a summary at the end lists any that failed.

Tables (CSV) and figures (PNG, plus PDF in publication style) go to `runs/<name>/output/`, and the log to `runs/<name>/run.log`. Every rate is given per site (TOF, s⁻¹) and also per area in the units of the hydrogen-surface literature (mol H₂ m⁻² s⁻¹, mol H₂ m⁻² s⁻¹ Pa⁻¹, the sticking probability; `…_per_area` figures; see [docs/rates_units_equilibrium.md](docs/rates_units_equilibrium.md)). A run of every analysis in `settings.py` ends with a PDF report, `output/report/<run>_report.pdf`: one section per analysis with its goal, key numbers, bordered tables, figures and how to read each figure panel by panel (needs LaTeX; `--no-report` skips it). Every solve is cached, so a rerun only solves what changed. `run.sh` keeps the Mac awake during long runs; the full Ni analysis takes about 2 hours on 8 cores the first time.

Or from Python / Jupyter (see `notebooks/analyze.ipynb`):

```python
import hydrogen_surface_kinetics as hsk
from hydrogen_surface_kinetics import analyses

ctx = analyses.run('runs/ni_first_try/settings.py', only=['steady_state', 'curves'])
ctx.outputs['curves']['tolerance']
```

## A new run

Copy a `settings.py` (e.g. `runs/try_2/Ni/settings.py`) to `runs/<new name>/settings.py` and change what differs. The minimum is the mechanism file and the species roles:

```python
mechanism = 'chem_annotated.yaml'             # relative to settings.py
expected_site_density = (2.483e-9, 'mol/cm2')   # from the RMG input; catches unit errors

probe_gas = 'H2(6)'            # names as in the mechanism
empty_site = 'site(10)'
probe_adsorbate = 'HX(11)'     # the probe reaction  H2 + 2 site <=> 2 HX  is found automatically
impurities = {'CO': {'C': 1, 'O': 1}, 'NH3': {'N': 1, 'H': 3}}   # by composition or by name
impurity_ratio = 0.25                  # p(impurity)/p(H2); or one per impurity, e.g. {'CO': 1/190, 'NH3': 0.1/190}

temperatures = [650.0, 800.0, 1000.0]   # temperature series (K), at reference_pressure
reference_pressure = 1e5                # Pa
pressures = [1e3, 1e5, 1e7]             # pressure series (Pa) ...
pressure_temperatures = [650.0]         # ... at these temperatures
analyses = ['steady_state', 'composition', 'curves']
```

Everything else has defaults. Then `bash run.sh runs/<new name> check`: it lists both series and how many solves each analysis needs before anything runs.

## Documentation

[`docs/simulation_workflow.md`](docs/simulation_workflow.md) explains every stage, the equations, the design choices, every setting and every output column, and has a troubleshooting table.

## Tests

```bash
~/anaconda3/envs/rmg4.0/bin/python -m pytest
```

The tests use a small fictitious mechanism (`tests/data/toy_mechanism.yaml`) and take about 15 minutes, most of it one end-to-end run of every analysis and its PDF report.
