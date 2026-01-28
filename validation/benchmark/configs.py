"""Benchmark configuration: jet structures, media, radiation, and observer settings."""

from dataclasses import dataclass
from typing import Any, Dict, List

# ---------------------------------------------------------------------------
# Dataclasses
# ---------------------------------------------------------------------------

@dataclass
class JetConfig:
    name: str
    factory: str
    params: Dict[str, Any]
    supports_spreading: bool = True
    supports_magnetar: bool = True

@dataclass
class MediumConfig:
    name: str
    factory: str
    params: Dict[str, Any]

@dataclass
class RadiationConfig:
    name: str
    params: Dict[str, Any]

# ---------------------------------------------------------------------------
# Jet configurations
# ---------------------------------------------------------------------------

STANDARD_JET_PARAMS = {"theta_c": 0.1, "E_iso": 1e52, "Gamma0": 300}

JET_CONFIGS: Dict[str, JetConfig] = {c.name: c for c in [
    JetConfig("tophat",        "TophatJet",        {**STANDARD_JET_PARAMS}),
    JetConfig("gaussian",      "GaussianJet",      {**STANDARD_JET_PARAMS}),
    JetConfig("powerlaw",      "PowerLawJet",      {**STANDARD_JET_PARAMS, "k_e": 2.0, "k_g": 2.0}),
    JetConfig("two_component", "TwoComponentJet",  {"theta_c": 0.05, "E_iso": 1e52, "Gamma0": 300,
                                                     "theta_w": 0.3, "E_iso_w": 1e50, "Gamma0_w": 50}),
]}

# ---------------------------------------------------------------------------
# Medium configurations
# ---------------------------------------------------------------------------

_STD_RAD = {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3, "xi_e": 1.0}

MEDIUM_CONFIGS: Dict[str, MediumConfig] = {c.name: c for c in [
    MediumConfig("ISM",            "ISM",  {"n_ism": 1.0}),
    MediumConfig("dense_ISM",      "ISM",  {"n_ism": 100.0}),
    MediumConfig("diffuse_ISM",    "ISM",  {"n_ism": 0.01}),
    MediumConfig("wind",           "Wind", {"A_star": 0.1}),
    MediumConfig("dense_wind",     "Wind", {"A_star": 1.0}),
    MediumConfig("light_wind",     "Wind", {"A_star": 0.01}),
    MediumConfig("wind_ism_floor", "Wind", {"A_star": 0.1, "n_ism": 0.1}),
]}

# ---------------------------------------------------------------------------
# Radiation configurations
# ---------------------------------------------------------------------------

RADIATION_CONFIGS: Dict[str, RadiationConfig] = {c.name: c for c in [
    RadiationConfig("synchrotron_only", {**_STD_RAD, "ssc_cooling": False, "ssc": False, "kn": False}),
    RadiationConfig("with_ssc_cooling", {**_STD_RAD, "ssc_cooling": True,  "ssc": False, "kn": False}),
    RadiationConfig("full_ssc",         {**_STD_RAD, "ssc_cooling": True,  "ssc": True,  "kn": False}),
    RadiationConfig("ssc_kn",           {**_STD_RAD, "ssc_cooling": True,  "ssc": True,  "kn": True}),
    RadiationConfig("fast_cooling",     {"eps_e": 0.1, "eps_B": 0.1,  "p": 2.5, "xi_e": 1e-3, "ssc_cooling": False, "ssc": False, "kn": False}),
    RadiationConfig("steep_spectrum",   {"eps_e": 0.1, "eps_B": 0.01, "p": 2.8, "xi_e": 1.0,  "ssc_cooling": False, "ssc": False, "kn": False}),
    RadiationConfig("flat_spectrum",    {"eps_e": 0.1, "eps_B": 0.01, "p": 2.05,"xi_e": 1.0,  "ssc_cooling": False, "ssc": False, "kn": False}),
    RadiationConfig("partial_xi_e",     {"eps_e": 0.1, "eps_B": 0.01, "p": 2.3, "xi_e": 0.1,  "ssc_cooling": False, "ssc": False, "kn": False}),
]}

# ---------------------------------------------------------------------------
# Observer
# ---------------------------------------------------------------------------

OBSERVER_CONFIG = {"lumi_dist": 1e28, "z": 1.0}

# ---------------------------------------------------------------------------
# Lookup helpers
# ---------------------------------------------------------------------------

def _lookup(registry, name):
    if name not in registry:
        raise ValueError(f"Unknown key: {name}. Available: {list(registry.keys())}")
    return registry[name]

def get_jet_config(name: str) -> JetConfig:            return _lookup(JET_CONFIGS, name)
def get_medium_config(name: str) -> MediumConfig:      return _lookup(MEDIUM_CONFIGS, name)
def get_radiation_config(name: str) -> RadiationConfig: return _lookup(RADIATION_CONFIGS, name)

def get_all_jet_names() -> List[str]:       return list(JET_CONFIGS)
def get_all_medium_names() -> List[str]:    return list(MEDIUM_CONFIGS)
def get_all_radiation_names() -> List[str]: return list(RADIATION_CONFIGS)
def get_radiation_names_no_ssc() -> List[str]:
    return [n for n, c in RADIATION_CONFIGS.items() if not c.params.get("ssc", False)]

# ---------------------------------------------------------------------------
# Factory functions
# ---------------------------------------------------------------------------

def create_jet(config: JetConfig, spreading=False, magnetar=None):
    import VegasAfterglow as va
    params = config.params.copy()
    if config.supports_spreading:
        params["spreading"] = spreading
    if config.supports_magnetar and magnetar is not None:
        params["magnetar"] = magnetar
    return getattr(va, config.factory)(**params)

def create_medium(config: MediumConfig):
    import VegasAfterglow as va
    return getattr(va, config.factory)(**config.params)

def create_radiation(config: RadiationConfig):
    import VegasAfterglow as va
    return va.Radiation(**config.params)

def create_observer(theta_obs: float = 0.0):
    import VegasAfterglow as va
    return va.Observer(**OBSERVER_CONFIG, theta_obs=theta_obs)
