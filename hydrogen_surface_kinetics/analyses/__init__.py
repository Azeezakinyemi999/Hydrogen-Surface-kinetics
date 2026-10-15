"""The analyses, in the order they run. Each module has run(ctx)."""
import importlib
import time

from ..config import ANALYSIS_ORDER, load_config
from ..model import Model
from .common import Context


def module(name):
    return importlib.import_module(f'.{name}', __name__)


def run(cfg_or_path=None, only=None, ctx=None, log=print, report=None, **run_overrides):
    """Run the analyses chosen in the config (or `only`, a name or list) and return the
    Context, whose .outputs holds every table and result.
    Pass ctx= (a Context from an earlier call) to add analyses to the same loaded model
    and results, e.g. from a notebook. run_overrides: workers, job_timeout, recompute.
    A run of the full list (only=None) ends with the PDF report when report['enabled'] is
    set in the config; report=True/False overrides that."""
    if ctx is None:
        cfg = load_config(cfg_or_path) if isinstance(cfg_or_path, str) else cfg_or_path
        model = Model(cfg)
        log(model.summary())
        ctx = Context(model, run_overrides, log=log)
    else:
        ctx.run_overrides.update({k: v for k, v in run_overrides.items() if v is not None})
    cfg = ctx.cfg
    names = cfg.analyses if only is None else ([only] if isinstance(only, str) else list(only))
    unknown = [n for n in names if n not in ANALYSIS_ORDER]
    if unknown:
        raise ValueError(f'unknown analyses {unknown}; choose from {ANALYSIS_ORDER}')
    t0 = time.time()
    for name in ANALYSIS_ORDER:
        if name in names:
            module(name).run(ctx)
    if (cfg.report['enabled'] and only is None) if report is None else report:
        build_report(ctx)
    ctx.log(f'\ndone in {time.time() - t0:.0f} s; outputs in {cfg.output_dir}')
    return ctx


def build_report(ctx):
    """The run's PDF report. A failure is logged, never raised: the results are already saved
    and the report can be rebuilt from them (`bash run.sh runs/<name> report`)."""
    from .. import report
    ctx.log('\n== report ==')
    try:
        return report.build(ctx.cfg, ctx.model, log=ctx.log)
    except Exception as e:                               # noqa: BLE001
        ctx.log(f'  report not built: {type(e).__name__}: {str(e)[:2000]}\n'
                f'  the results are saved; rebuild it with: bash run.sh {ctx.cfg.run_dir} report')
        return None
