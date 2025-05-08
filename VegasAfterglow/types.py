import numpy as np
from dataclasses import dataclass
from typing import Sequence, Tuple, Callable, Any, Optional
from enum import Enum
from .VegasAfterglowC import ModelParams, Setups, ObsData, VegasMC, Emission

@dataclass
class FitResult:
    """
    The result of an MCMC fit.
    """
    samples: np.ndarray
    log_probs: np.ndarray
    best_params: np.ndarray
    labels: Sequence[str]

class Scale(Enum):
    LINEAR = "linear"
    LOG    = "log"
    FIXED  = "fixed"

@dataclass
class ParamDef:
    """
    Single-parameter definition for MCMC.
    scale=LOG means we sample log10(x), then transform via 10**v.
    scale=FIXED means this param never appears in the sampler.
    """
    name:   str
    init:   float
    lower:  Optional[float] = None  # for FIXED these can be None
    upper:  Optional[float] = None
    scale:  Scale = Scale.LINEAR