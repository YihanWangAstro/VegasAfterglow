"""Configuration constants for afterglow model fitting."""

import math
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, FrozenSet, Tuple

import emcee
import numpy as np

from ..VegasAfterglowC import (
    ISM,
    GaussianJet,
    PowerLawJet,
    PowerLawWing,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
)

LATEX_LABELS = {
    "E_iso": r"$E_{\rm iso}$",
    "Gamma0": r"$\Gamma_0$",
    "theta_c": r"$\theta_c$",
    "theta_v": r"$\theta_v$",
    "theta_w": r"$\theta_w$",
    "k_e": r"$k_E$",
    "k_g": r"$k_\Gamma$",
    "E_iso_w": r"$E_{\rm iso,w}$",
    "Gamma0_w": r"$\Gamma_{0,w}$",
    "n_ism": r"$n_{\rm ISM}$",
    "A_star": r"$A_*$",
    "n0": r"$n_0$",
    "k_m": r"$k_m$",
    "p": r"$p$",
    "eps_e": r"$\epsilon_e$",
    "eps_B": r"$\epsilon_B$",
    "xi_e": r"$\xi_e$",
    "p_r": r"$p_r$",
    "eps_e_r": r"$\epsilon_{e,r}$",
    "eps_B_r": r"$\epsilon_{B,r}$",
    "xi_e_r": r"$\xi_{e,r}$",
    "tau": r"$\tau$",
    "L0": r"$L_0$",
    "t0": r"$t_0$",
    "q": r"$q$",
}


@dataclass(frozen=True)
class JetSpec:
    """Everything the fitting pipeline needs to know about one jet type.

    ``params`` is the ordered list of constructor arguments pulled from
    ``ModelParams`` -- it doubles as the set of sampler parameters this jet
    requires. ``fixed_kwargs`` are constructor arguments with a hard-coded
    value (not sampled). Parameters in ``_SWAPPABLE_JET_PARAMS`` that a jet
    does not require are forbidden for it (they would be silently unused).
    The opening angles theta_c / theta_w are never forbidden.
    """

    constructor: Callable
    params: Tuple[str, ...]
    fixed_kwargs: Dict[str, Any] = field(default_factory=dict)
    supports_magnetar: bool = True

    @property
    def required(self) -> FrozenSet[str]:
        return frozenset(self.params)

    @property
    def forbidden(self) -> FrozenSet[str]:
        return _SWAPPABLE_JET_PARAMS - self.required


@dataclass(frozen=True)
class MediumSpec:
    """Everything the fitting pipeline needs to know about one medium type.

    Unlike jets, ``required`` is not ``set(params)``: the wind medium always
    passes n_ism/n0/k_m through from ``ModelParams`` (which has physical
    defaults for them), so they are constructor arguments without being
    required sampler parameters.
    """

    constructor: Callable
    params: Tuple[str, ...]
    required: FrozenSet[str]
    forbidden: FrozenSet[str]


#: Structure parameters that switch on/off with the jet type; anything here
#: that a jet does not require is rejected by parameter validation.
_SWAPPABLE_JET_PARAMS = frozenset(
    {"E_iso", "Gamma0", "k_e", "k_g", "E_iso_w", "Gamma0_w"}
)

#: The jet registry: adding a jet type to the fitting pipeline is one entry here.
JETS = {
    "tophat": JetSpec(TophatJet, ("theta_c", "E_iso", "Gamma0")),
    "gaussian": JetSpec(GaussianJet, ("theta_c", "E_iso", "Gamma0")),
    "powerlaw": JetSpec(PowerLawJet, ("theta_c", "E_iso", "Gamma0", "k_e", "k_g")),
    "two_component": JetSpec(
        TwoComponentJet,
        ("theta_c", "E_iso", "Gamma0", "theta_w", "E_iso_w", "Gamma0_w"),
    ),
    "step_powerlaw": JetSpec(
        StepPowerLawJet,
        ("theta_c", "E_iso", "Gamma0", "E_iso_w", "Gamma0_w", "k_e", "k_g"),
    ),
    "powerlaw_wing": JetSpec(
        PowerLawWing,
        ("theta_c", "E_iso_w", "Gamma0_w", "k_e", "k_g"),
        supports_magnetar=False,
    ),
    "uniform": JetSpec(
        TophatJet,
        ("E_iso", "Gamma0"),
        fixed_kwargs={"theta_c": math.pi / 2},
    ),
}

#: The medium registry, same idea as JETS.
MEDIA = {
    "ism": MediumSpec(
        ISM,
        params=("n_ism",),
        required=frozenset({"n_ism"}),
        forbidden=frozenset({"A_star", "n0", "k_m"}),
    ),
    "wind": MediumSpec(
        Wind,
        params=("A_star", "n_ism", "n0", "k_m"),
        required=frozenset({"A_star"}),
        forbidden=frozenset(),
    ),
}

TOGGLE_RULES = {
    "forward_shock": ({"eps_e", "eps_B", "p"}, set()),
    "rvs_shock": (
        {"p_r", "eps_e_r", "eps_B_r", "tau"},
        {"p_r", "eps_e_r", "eps_B_r", "xi_e_r"},
    ),
    "magnetar": ({"L0", "t0", "q"}, {"L0", "t0", "q"}),
}

SAMPLER_DEFAULTS = {
    "dynesty": {
        "nlive": 500,
        "dlogz": 0.1,
        "sample": "rslice",
        "maxmcmc": 5000,
    },
    "emcee": {
        "nsteps": 5000,
        "nburn": 1000,
        "thin": 1,
        "moves": [
            (emcee.moves.DEMove(), 0.7),
            (emcee.moves.DESnookerMove(), 0.3),
        ],
    },
}


@dataclass
class _BandObs:
    """Band-integrated flux observation group."""

    nu_min: float
    nu_max: float
    num_points: int
    t: np.ndarray
    flux: np.ndarray
    err: np.ndarray
    weights: np.ndarray
    name: "str | None" = None
