"""The analyses, in the order they run. Each module has run(ctx)."""
import importlib
import time

from ..config import ANALYSIS_ORDER, load_config
from ..model import Model
from .common import Context


def module(name):
    return importlib.import_module(f'.{name}', __name__)


def run(cfg_or_path, only=None, log=print, **run_overrides):
    """Run the analyses chosen in the config (or `only`, a name or list) and return the
    Context, whose .outputs holds every table and result.
    run_overrides: workers, job_timeout, recompute."""
    cfg = load_config(cfg_or_path) if isinstance(cfg_or_path, str) else cfg_or_path
    names = cfg.analyses if only is None else ([only] if isinstance(only, str) else list(only))
    unknown = [n for n in names if n not in ANALYSIS_ORDER]
    if unknown:
        raise ValueError(f'unknown analyses {unknown}; choose from {ANALYSIS_ORDER}')
    model = Model(cfg)
    log(model.summary())
    ctx = Context(model, run_overrides, log=log)
    t0 = time.time()
    for name in ANALYSIS_ORDER:
        if name in names:
            module(name).run(ctx)
    log(f'\nall done in {time.time() - t0:.0f} s; outputs in {cfg.output_dir}')
    return ctx
