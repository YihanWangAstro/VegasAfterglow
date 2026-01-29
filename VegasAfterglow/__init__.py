from ._version import __version__, __version_tuple__

# Core physics types (always available)
from .types import (
    ISM,
    Ejecta,
    GaussianJet,
    Magnetar,
    Medium,
    Model,
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
)

# MCMC types (require: pip install VegasAfterglow[mcmc])
try:
    from .runner import AfterglowLikelihood, Fitter
    from .types import FitResult, ModelParams, ObsData, ParamDef, Scale, Setups, VegasMC
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
    ObsData = _mcmc_stub("ObsData")
    Setups = _mcmc_stub("Setups")
    VegasMC = _mcmc_stub("VegasMC")
    ModelParams = _mcmc_stub("ModelParams")
    FitResult = _mcmc_stub("FitResult")
    ParamDef = _mcmc_stub("ParamDef")
    Scale = _mcmc_stub("Scale")

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
    # MCMC (requires VegasAfterglow[mcmc])
    "ModelParams",
    "Setups",
    "ObsData",
    "VegasMC",
    "FitResult",
    "Fitter",
    "ParamDef",
    "Scale",
    "AfterglowLikelihood",
]
