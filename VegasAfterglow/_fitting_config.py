"""Configuration constants for afterglow model fitting."""

from dataclasses import dataclass

import emcee
import numpy as np

from .VegasAfterglowC import (
    GaussianJet,
    PowerLawJet,
    PowerLawWing,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
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

MEDIUM_RULES = {
    "ism": ({"n_ism"}, {"A_star", "n0", "k_m"}),
    "wind": ({"A_star"}, set()),
}

JET_RULES = {
    "tophat": ({"theta_c", "E_iso", "Gamma0"}, {"k_e", "k_g", "E_iso_w", "Gamma0_w"}),
    "gaussian": ({"theta_c", "E_iso", "Gamma0"}, {"k_e", "k_g", "E_iso_w", "Gamma0_w"}),
    "powerlaw": ({"theta_c", "E_iso", "Gamma0", "k_e", "k_g"}, {"E_iso_w", "Gamma0_w"}),
    "two_component": (
        {"theta_c", "E_iso", "Gamma0", "theta_w", "E_iso_w", "Gamma0_w"},
        {"k_e", "k_g"},
    ),
    "step_powerlaw": (
        {"theta_c", "E_iso", "Gamma0", "E_iso_w", "Gamma0_w", "k_e", "k_g"},
        set(),
    ),
    "powerlaw_wing": (
        {"theta_c", "E_iso_w", "Gamma0_w", "k_e", "k_g"},
        {"E_iso", "Gamma0"},
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

_JET_CONSTRUCTORS = {
    "tophat": (TophatJet, ["theta_c", "E_iso", "Gamma0"]),
    "gaussian": (GaussianJet, ["theta_c", "E_iso", "Gamma0"]),
    "powerlaw": (PowerLawJet, ["theta_c", "E_iso", "Gamma0", "k_e", "k_g"]),
    "two_component": (
        TwoComponentJet,
        ["theta_c", "E_iso", "Gamma0", "theta_w", "E_iso_w", "Gamma0_w"],
    ),
    "step_powerlaw": (
        StepPowerLawJet,
        ["theta_c", "E_iso", "Gamma0", "E_iso_w", "Gamma0_w", "k_e", "k_g"],
    ),
    "powerlaw_wing": (PowerLawWing, ["theta_c", "E_iso_w", "Gamma0_w", "k_e", "k_g"]),
    "uniform": (TophatJet, ["E_iso", "Gamma0"]),
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
