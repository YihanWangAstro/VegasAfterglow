"""
Type definitions and C++ bindings for VegasAfterglow.

This module imports the core C++ classes for afterglow physics calculations.
"""

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
)
