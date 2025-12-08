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
)


@dataclass
class FitResult:
    """
    The result of a sampling fit.
    """

    samples: np.ndarray
    log_probs: np.ndarray
    labels: Sequence[str]
    latex_labels: Sequence[str] = None  # LaTeX labels for plotting
    top_k_params: np.ndarray = None
    top_k_log_probs: np.ndarray = None
    bilby_result: object = None  # Full bilby Result object (for diagnostics)


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
