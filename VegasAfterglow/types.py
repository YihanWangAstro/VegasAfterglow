import math
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
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
    logscale_screen,
)


class ModelParams:
    """Complete parameter set for GRB afterglow model configuration.

    All values are in user-facing units (CGS where applicable):
    angles in radians, energies in erg, densities in cm^-3, times in seconds.

    Custom attributes can be added freely for use with custom jet/medium factories
    (e.g., ``params.r_scale = 1.5``).
    """

    def __init__(self):
        self.theta_v = 0.0

        # External medium
        self.n_ism = 0.0
        self.n0 = math.inf
        self.A_star = 0.0
        self.k_m = 2.0

        # Jet core
        self.E_iso = 1e52
        self.Gamma0 = 300.0
        self.theta_c = 0.1
        self.k_e = 2.0
        self.k_g = 2.0
        self.tau = 1.0  # Central engine duration [seconds]

        # Jet wing
        self.E_iso_w = 1e52
        self.Gamma0_w = 300.0
        self.theta_w = math.pi / 2

        # Magnetar
        self.L0 = 0.0
        self.t0 = 1.0
        self.q = 2.0

        # Forward shock radiation
        self.p = 2.3
        self.eps_e = 0.1
        self.eps_B = 0.01
        self.xi_e = 1.0

        # Reverse shock radiation
        self.p_r = 2.3
        self.eps_e_r = 0.1
        self.eps_B_r = 0.01
        self.xi_e_r = 1.0


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
    initial: float = (
        None  # Initial value (in linear space, auto-converted for LOG scale)
    )
