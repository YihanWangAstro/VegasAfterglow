"""
VegasAfterglow: Fast C++ afterglow physics library

A pure physics library for computing GRB afterglow light curves with:
- Multiple jet structures (tophat, gaussian, power law, etc.)
- ISM and Wind medium profiles
- Forward and reverse shock emission
- Inverse Compton scattering with Klein-Nishina corrections

For parameter estimation and MCMC fitting, use redback:
https://github.com/nikhil-sarin/redback
"""

from ._version import __version__, __version_tuple__
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

__all__ = [
    "__version__",
    "__version_tuple__",
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
]
