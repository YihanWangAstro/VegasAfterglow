from . import units
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


# MCMC fitting — lazy-loaded to avoid importing bilby/scipy on every startup.
# Actual import happens on first access to Fitter or AfterglowLikelihood.
def __getattr__(name):
    if name in ("Fitter", "AfterglowLikelihood"):
        try:
            from .runner import AfterglowLikelihood, Fitter
        except ImportError:
            raise ImportError(
                "MCMC fitting requires additional dependencies. "
                "Install them with: pip install VegasAfterglow[mcmc]"
            )
        # Cache in module namespace so __getattr__ is only called once
        globals()["Fitter"] = Fitter
        globals()["AfterglowLikelihood"] = AfterglowLikelihood
        return globals()[name]
    raise AttributeError(f"module 'VegasAfterglow' has no attribute {name!r}")


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
    # Units
    "units",
    # MCMC fitting (requires VegasAfterglow[mcmc])
    "Fitter",
    "AfterglowLikelihood",
]
