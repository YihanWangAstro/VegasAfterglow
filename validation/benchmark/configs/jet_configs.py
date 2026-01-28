"""
Jet configuration definitions for benchmarking.

This module defines all jet structure configurations to be tested
in the benchmark suite, covering all 6 jet types with various options.
"""

from dataclasses import dataclass
from typing import Any, Dict, List


@dataclass
class JetConfig:
    """Configuration for a jet structure."""
    name: str
    factory: str
    params: Dict[str, Any]
    description: str = ""
    supports_spreading: bool = True
    supports_magnetar: bool = True


# Standard parameters shared across jet types
STANDARD_JET_PARAMS = {
    "theta_c": 0.1,       # Core angle [rad] ~ 5.7 degrees
    "E_iso": 1e52,        # Isotropic equivalent energy [erg]
    "Gamma0": 300,        # Initial Lorentz factor
}

# Wing parameters for structured jets
WING_PARAMS = {
    "theta_w": 0.3,       # Wing angle [rad] ~ 17 degrees
    "E_iso_w": 1e50,      # Wing energy [erg]
    "Gamma0_w": 50,       # Wing Lorentz factor
}

# Power-law indices
POWERLAW_INDICES = {
    "k_e": 2.0,           # Energy angular falloff
    "k_g": 2.0,           # Lorentz factor angular falloff
}


# ============================================================================
# All 6 Jet Structure Types
# ============================================================================

JET_CONFIGS = [
    # 1. TophatJet - Uniform core, sharp edges
    JetConfig(
        name="tophat",
        factory="TophatJet",
        params={
            "theta_c": 0.1,
            "E_iso": 1e52,
            "Gamma0": 300,
        },
        description="Uniform top-hat jet with sharp edges",
    ),

    # 2. GaussianJet - Smooth Gaussian profile
    JetConfig(
        name="gaussian",
        factory="GaussianJet",
        params={
            "theta_c": 0.1,
            "E_iso": 1e52,
            "Gamma0": 300,
        },
        description="Gaussian angular energy/Lorentz factor profile",
    ),

    # 3. PowerLawJet - Power-law angular profile
    JetConfig(
        name="powerlaw",
        factory="PowerLawJet",
        params={
            "theta_c": 0.1,
            "E_iso": 1e52,
            "Gamma0": 300,
            "k_e": 2.0,
            "k_g": 2.0,
        },
        description="Power-law angular profile: E ∝ θ^-k_e, Γ ∝ θ^-k_g",
    ),

    # 4. TwoComponentJet - Core + uniform wing
    JetConfig(
        name="two_component",
        factory="TwoComponentJet",
        params={
            "theta_c": 0.05,
            "E_iso": 1e52,
            "Gamma0": 300,
            "theta_w": 0.3,
            "E_iso_w": 1e50,
            "Gamma0_w": 50,
        },
        description="Narrow core + wide uniform wing",
    ),
]


# ============================================================================
# Jet Variations for Comprehensive Testing
# ============================================================================

# Different core angles
CORE_ANGLE_VARIATIONS = [
    {"name": "narrow_core", "theta_c": 0.03},  # ~1.7 degrees
    {"name": "standard_core", "theta_c": 0.1},  # ~5.7 degrees
    {"name": "wide_core", "theta_c": 0.2},     # ~11.5 degrees
]

# Different initial Lorentz factors
LORENTZ_VARIATIONS = [
    {"name": "low_gamma", "Gamma0": 100},
    {"name": "medium_gamma", "Gamma0": 300},
    {"name": "high_gamma", "Gamma0": 1000},
]

# Different energy levels
ENERGY_VARIATIONS = [
    {"name": "low_energy", "E_iso": 1e51},
    {"name": "standard_energy", "E_iso": 1e52},
    {"name": "high_energy", "E_iso": 1e53},
]

# Power-law index variations
POWERLAW_INDEX_VARIATIONS = [
    {"name": "shallow_pl", "k_e": 1.0, "k_g": 1.0},
    {"name": "standard_pl", "k_e": 2.0, "k_g": 2.0},
    {"name": "steep_pl", "k_e": 4.0, "k_g": 4.0},
]


def create_jet(config: JetConfig, spreading: bool = False, magnetar=None):
    """
    Create a jet object from configuration.

    Parameters
    ----------
    config : JetConfig
        Jet configuration
    spreading : bool
        Enable lateral spreading
    magnetar : Magnetar, optional
        Magnetar energy injection

    Returns
    -------
    Ejecta
        Jet/ejecta object
    """
    import VegasAfterglow as va

    factory = getattr(va, config.factory)
    params = config.params.copy()

    # Add spreading if supported
    if config.supports_spreading:
        params["spreading"] = spreading

    # Add magnetar if supported and provided
    if config.supports_magnetar and magnetar is not None:
        params["magnetar"] = magnetar

    return factory(**params)


def create_jet_with_variations(base_config: JetConfig, **overrides):
    """
    Create a jet with parameter overrides.

    Parameters
    ----------
    base_config : JetConfig
        Base jet configuration
    **overrides : dict
        Parameters to override

    Returns
    -------
    Ejecta
        Jet object with overridden parameters
    """
    import VegasAfterglow as va

    factory = getattr(va, base_config.factory)
    params = base_config.params.copy()
    params.update(overrides)

    return factory(**params)


def get_all_jet_names() -> List[str]:
    """Return list of all jet configuration names."""
    return [c.name for c in JET_CONFIGS]


def get_jet_config(name: str) -> JetConfig:
    """Get jet configuration by name."""
    for config in JET_CONFIGS:
        if config.name == name:
            return config
    raise ValueError(f"Unknown jet type: {name}. Available: {get_all_jet_names()}")


def get_spreading_capable_jets() -> List[str]:
    """Return list of jet types that support spreading."""
    return [c.name for c in JET_CONFIGS if c.supports_spreading]


def get_magnetar_capable_jets() -> List[str]:
    """Return list of jet types that support magnetar injection."""
    return [c.name for c in JET_CONFIGS if c.supports_magnetar]
