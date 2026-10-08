"""Steady-state and transient analysis of surface reaction mechanisms with Cantera.

Typical use:

    import hydrogen_surface_kinetics as hsk
    cfg = hsk.load_config('runs/my_run/input.py')
    model = hsk.Model(cfg)
"""
import os as _os

# One thread per worker process. Numerical libraries read these when NumPy is first
# imported, so they must be set before that; existing values are respected.
for _var in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'):
    _os.environ.setdefault(_var, '1')

from .config import Config, ConfigError, load_config  # noqa: E402
from .model import Model, ModelError  # noqa: E402

__version__ = '0.1.0'
__all__ = ['Config', 'ConfigError', 'load_config', 'Model', 'ModelError']
