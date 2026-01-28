"""
Physics configuration definitions for benchmarking.

This module defines medium types, radiation settings, resolution
configurations, and viewing angle setups for comprehensive testing.
"""

from dataclasses import dataclass
from typing import Any, Dict, List, Tuple


@dataclass
class MediumConfig:
    """Configuration for ambient medium."""
    name: str
    factory: str
    params: Dict[str, Any]
    description: str = ""


@dataclass
class RadiationConfig:
    """Configuration for radiation processes."""
    name: str
    params: Dict[str, Any]
    description: str = ""
    complexity: str = "low"  # low, medium, high (for timing expectations)


@dataclass
class ResolutionConfig:
    """Configuration for numerical resolution."""
    name: str
    values: Tuple[float, float, float]  # (phi_ppd, theta_ppd, t_ppd)
    description: str = ""


@dataclass
class ViewingAngleConfig:
    """Configuration for observer viewing angle."""
    name: str
    theta_obs: float
    description: str = ""


# ============================================================================
# Medium Configurations
# ============================================================================

MEDIUM_CONFIGS = [
    # Standard ISM
    MediumConfig(
        name="ISM",
        factory="ISM",
        params={"n_ism": 1.0},
        description="Uniform ISM with n=1 cm^-3",
    ),
    # Dense ISM (like molecular cloud)
    MediumConfig(
        name="dense_ISM",
        factory="ISM",
        params={"n_ism": 100.0},
        description="Dense ISM with n=100 cm^-3",
    ),
    # Diffuse ISM
    MediumConfig(
        name="diffuse_ISM",
        factory="ISM",
        params={"n_ism": 0.01},
        description="Diffuse ISM with n=0.01 cm^-3",
    ),
    # Standard stellar wind
    MediumConfig(
        name="wind",
        factory="Wind",
        params={"A_star": 0.1},
        description="Standard wind with A*=0.1",
    ),
    # Dense wind
    MediumConfig(
        name="dense_wind",
        factory="Wind",
        params={"A_star": 1.0},
        description="Dense wind with A*=1.0",
    ),
    # Light wind
    MediumConfig(
        name="light_wind",
        factory="Wind",
        params={"A_star": 0.01},
        description="Light wind with A*=0.01",
    ),
    # Wind with ISM floor
    MediumConfig(
        name="wind_ism_floor",
        factory="Wind",
        params={"A_star": 0.1, "n_ism": 0.1},
        description="Wind transitioning to ISM floor",
    ),
]


# ============================================================================
# Radiation Configurations
# ============================================================================

# Standard microphysics parameters
STANDARD_RADIATION = {
    "eps_e": 0.1,
    "eps_B": 0.01,
    "p": 2.3,
    "xi_e": 1.0,
}

RADIATION_CONFIGS = [
    # Basic synchrotron only
    RadiationConfig(
        name="synchrotron_only",
        params={
            **STANDARD_RADIATION,
            "ssc_cooling": False,
            "ssc": False,
            "kn": False,
        },
        description="Pure synchrotron emission",
        complexity="low",
    ),

    # Synchrotron with IC cooling
    RadiationConfig(
        name="with_ssc_cooling",
        params={
            **STANDARD_RADIATION,
            "ssc_cooling": True,
            "ssc": False,
            "kn": False,
        },
        description="Synchrotron with inverse Compton cooling",
        complexity="medium",
    ),

    # Full SSC emission (no KN)
    RadiationConfig(
        name="full_ssc",
        params={
            **STANDARD_RADIATION,
            "ssc_cooling": True,
            "ssc": True,
            "kn": False,
        },
        description="Synchrotron self-Compton emission",
        complexity="high",
    ),

    # Full SSC with Klein-Nishina
    RadiationConfig(
        name="ssc_kn",
        params={
            **STANDARD_RADIATION,
            "ssc_cooling": True,
            "ssc": True,
            "kn": True,
        },
        description="SSC with Klein-Nishina corrections",
        complexity="high",
    ),

    # High eps_B (fast cooling regime)
    RadiationConfig(
        name="fast_cooling",
        params={
            "eps_e": 0.1,
            "eps_B": 0.1,
            "p": 2.5,
            "xi_e": 1e-3,
            "ssc_cooling": False,
            "ssc": False,
            "kn": False,
        },
        description="High eps_B for fast cooling",
        complexity="low",
    ),

    # Steep electron spectrum
    RadiationConfig(
        name="steep_spectrum",
        params={
            "eps_e": 0.1,
            "eps_B": 0.01,
            "p": 2.8,
            "xi_e": 1.0,
            "ssc_cooling": False,
            "ssc": False,
            "kn": False,
        },
        description="Steep electron spectrum p=2.8",
        complexity="low",
    ),

    # Flat electron spectrum
    RadiationConfig(
        name="flat_spectrum",
        params={
            "eps_e": 0.1,
            "eps_B": 0.01,
            "p": 2.05,
            "xi_e": 1.0,
            "ssc_cooling": False,
            "ssc": False,
            "kn": False,
        },
        description="Flat electron spectrum p=2.05",
        complexity="low",
    ),

    # Partial electron acceleration
    RadiationConfig(
        name="partial_xi_e",
        params={
            "eps_e": 0.1,
            "eps_B": 0.01,
            "p": 2.3,
            "xi_e": 0.1,
            "ssc_cooling": False,
            "ssc": False,
            "kn": False,
        },
        description="Partial electron acceleration xi_e=0.1",
        complexity="low",
    ),
]


# ============================================================================
# Resolution Configurations (for convergence testing)
# ============================================================================

RESOLUTION_CONFIGS = [
    # Very coarse (fast, for quick checks)
    ResolutionConfig(
        name="coarse",
        values=(0.1, 0.3, 5),
        description="Very coarse resolution (fastest)",
    ),
    # Low resolution
    ResolutionConfig(
        name="low",
        values=(0.2, 0.3, 7),
        description="Low resolution",
    ),
    # Medium resolution (default)
    ResolutionConfig(
        name="medium",
        values=(0.3, 0.3, 10),
        description="Medium resolution (default)",
    ),
    # High resolution
    ResolutionConfig(
        name="high",
        values=(0.5, 0.3, 15),
        description="High resolution",
    ),
    # Very high resolution
    ResolutionConfig(
        name="very_high",
        values=(0.5, 1, 20),
        description="Very high resolution",
    ),
]


# ============================================================================
# Viewing Angle Configurations
# ============================================================================

VIEWING_ANGLES = [
    ViewingAngleConfig(
        name="on_axis",
        theta_obs=0.0,
        description="On-axis observation",
    ),
    ViewingAngleConfig(
        name="slight_off",
        theta_obs=0.05,
        description="Slightly off-axis (0.5 × θ_c)",
    ),
    ViewingAngleConfig(
        name="edge_core",
        theta_obs=0.1,
        description="Edge of core (1 × θ_c)",
    ),
    ViewingAngleConfig(
        name="off_axis_2x",
        theta_obs=0.2,
        description="Off-axis (2 × θ_c)",
    ),
    ViewingAngleConfig(
        name="off_axis_3x",
        theta_obs=0.3,
        description="Strongly off-axis (3 × θ_c)",
    ),
    ViewingAngleConfig(
        name="off_axis_5x",
        theta_obs=0.5,
        description="Very off-axis (5 × θ_c)",
    ),
]


# ============================================================================
# Observer Configuration
# ============================================================================

OBSERVER_CONFIGS = {
    "nearby": {"lumi_dist": 1e26, "z": 0.01},    # ~30 Mpc
    "standard": {"lumi_dist": 1e28, "z": 1.0},   # ~3 Gpc
    "distant": {"lumi_dist": 3e28, "z": 2.0},    # ~10 Gpc
}

# Default observer
OBSERVER_CONFIG = OBSERVER_CONFIGS["standard"]


# ============================================================================
# Frequency Bands for Testing
# ============================================================================

FREQUENCY_BANDS = {
    "radio_low": 1e9,      # 1 GHz
    "radio_high": 1e10,    # 10 GHz
    "mm": 1e11,            # 100 GHz (3mm)
    "submm": 3e11,         # 300 GHz (1mm)
    "ir": 3e13,            # Far IR
    "optical": 4.84e14,    # R-band
    "uv": 1e15,            # UV
    "xray_soft": 2.4e17,   # 1 keV
    "xray_hard": 2.4e18,   # 10 keV
    "gamma": 2.4e20,       # 1 MeV
}


# ============================================================================
# Factory Functions
# ============================================================================

def create_medium(config: MediumConfig):
    """Create a medium object from configuration."""
    import VegasAfterglow as va

    factory = getattr(va, config.factory)
    return factory(**config.params)


def create_radiation(config: RadiationConfig):
    """Create a radiation object from configuration."""
    import VegasAfterglow as va

    return va.Radiation(**config.params)


def create_observer(theta_obs: float = 0.0, config_name: str = "standard"):
    """Create an observer object."""
    import VegasAfterglow as va

    obs_config = OBSERVER_CONFIGS.get(config_name, OBSERVER_CONFIG)
    return va.Observer(
        lumi_dist=obs_config["lumi_dist"],
        z=obs_config["z"],
        theta_obs=theta_obs
    )


# ============================================================================
# Getter Functions
# ============================================================================

def get_medium_config(name: str) -> MediumConfig:
    """Get medium configuration by name."""
    for config in MEDIUM_CONFIGS:
        if config.name == name:
            return config
    available = [c.name for c in MEDIUM_CONFIGS]
    raise ValueError(f"Unknown medium type: {name}. Available: {available}")


def get_radiation_config(name: str) -> RadiationConfig:
    """Get radiation configuration by name."""
    for config in RADIATION_CONFIGS:
        if config.name == name:
            return config
    available = [c.name for c in RADIATION_CONFIGS]
    raise ValueError(f"Unknown radiation config: {name}. Available: {available}")


def get_resolution_config(name: str) -> ResolutionConfig:
    """Get resolution configuration by name."""
    for config in RESOLUTION_CONFIGS:
        if config.name == name:
            return config
    available = [c.name for c in RESOLUTION_CONFIGS]
    raise ValueError(f"Unknown resolution config: {name}. Available: {available}")


def get_viewing_angle_config(name: str) -> ViewingAngleConfig:
    """Get viewing angle configuration by name."""
    for config in VIEWING_ANGLES:
        if config.name == name:
            return config
    available = [c.name for c in VIEWING_ANGLES]
    raise ValueError(f"Unknown viewing angle: {name}. Available: {available}")


def get_all_medium_names() -> List[str]:
    """Return all medium configuration names."""
    return [c.name for c in MEDIUM_CONFIGS]


def get_all_radiation_names() -> List[str]:
    """Return all radiation configuration names."""
    return [c.name for c in RADIATION_CONFIGS]


def get_radiation_names_no_ssc() -> List[str]:
    """Return radiation configuration names excluding full SSC emission.

    Excludes configs where ssc=True (actual SSC emission).
    Keeps configs with ssc_cooling=True but ssc=False (cooling only).
    """
    return [c.name for c in RADIATION_CONFIGS if not c.params.get("ssc", False)]


def get_all_resolution_names() -> List[str]:
    """Return all resolution configuration names."""
    return [c.name for c in RESOLUTION_CONFIGS]


def get_all_viewing_angle_names() -> List[str]:
    """Return all viewing angle configuration names."""
    return [c.name for c in VIEWING_ANGLES]
