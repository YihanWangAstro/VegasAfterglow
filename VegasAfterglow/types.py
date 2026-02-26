from dataclasses import dataclass
from enum import Enum
from typing import Sequence

import numpy as np

from .VegasAfterglowC import (  # noqa: F401
    ISM,
    Ejecta,
    GaussianJet,
    Magnetar,
    Medium,
    Model,
    ModelParams,
    ObsData,
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    Setups,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    VegasMC,
    Wind,
    powerLaw,
)


@dataclass
class FitResult:
    """
    The result of an MCMC fit.
    """

    samples: np.ndarray
    log_probs: np.ndarray
    labels: Sequence[str]
    top_k_params: np.ndarray
    top_k_log_probs: np.ndarray


class Scale(Enum):
    LINEAR = "linear"
    LOG = "log"
    FIXED = "fixed"


@dataclass
class ParamDef:
    """
    Single-parameter definition for MCMC.
    scale=LOG means we sample log10(x), then transform via 10**v.
    scale=FIXED means this param never appears in the sampler.
    """

    name: str
    lower: float
    upper: float
    scale: Scale = Scale.LINEAR
