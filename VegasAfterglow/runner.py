"""Backwards-compatibility shim for ``from VegasAfterglow.runner import Fitter``.

The fitter implementation moved into :mod:`VegasAfterglow.fitting`; this module
re-exports the public surface so existing code keeps working unchanged.
"""

from .fitting import AfterglowLikelihood, Fitter  # noqa: F401

__all__ = ["Fitter", "AfterglowLikelihood"]
