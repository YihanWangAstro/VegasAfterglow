from ._version import __version__, __version_tuple__
from .native import NativeFunc, gil_free

# Core physics types (always available)
from .types import (
    ISM,
    Ejecta,
    FitResult,
    GaussianJet,
    Magnetar,
    Medium,
    Model,
    ModelParams,
    Observer,
    ParamDef,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    Scale,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
    logscale_screen,
)

# MCMC fitting (requires: pip install VegasAfterglow[mcmc])
try:
    from .runner import AfterglowLikelihood, Fitter
except ImportError:
    _MCMC_MISSING_MSG = (
        "MCMC fitting requires additional dependencies. "
        "Install them with: pip install VegasAfterglow[mcmc]"
    )

    def _mcmc_stub(name):
        class _Stub:
            def __init__(self, *args, **kwargs):
                raise ImportError(_MCMC_MISSING_MSG)

        _Stub.__name__ = _Stub.__qualname__ = name
        return _Stub

    Fitter = _mcmc_stub("Fitter")
    AfterglowLikelihood = _mcmc_stub("AfterglowLikelihood")

__all__ = [
    "__version__",
    "__version_tuple__",
    # Core physics
    "ISM",
    "Wind",
    "Medium",
    "TophatJet",
    "GaussianJet",
    "PowerLawJet",
    "PowerLawWing",
    "TwoComponentJet",
    "StepPowerLawJet",
    "Ejecta",
    "Model",
    "Radiation",
    "Observer",
    "Magnetar",
    # Utilities
    "NativeFunc",
    "gil_free",
    "logscale_screen",
    # Parameter types
    "ModelParams",
    "FitResult",
    "ParamDef",
    "Scale",
    # MCMC fitting (requires VegasAfterglow[mcmc])
    "Fitter",
    "AfterglowLikelihood",
]
