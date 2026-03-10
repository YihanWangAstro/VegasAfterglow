"""Shared constants, parsing utilities, and conversion functions."""

from __future__ import annotations

import re
from typing import Any

import numpy as np

from VegasAfterglow import units
from VegasAfterglow.cli import _parse_frequency
from VegasAfterglow.units import _NAMED_BANDS

from .schemas import SharedParams

# ---------------------------------------------------------------------------
# Constants (previously in constants.py)
# ---------------------------------------------------------------------------

UNIT_SUFFIXES = {
    "Hz": 1, "kHz": units.kHz, "MHz": units.MHz, "GHz": units.GHz,
    "eV": units.eV, "keV": units.keV, "MeV": units.MeV, "GeV": units.GeV,
}

OBS_FLUX_UNITS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm\u00b2/s"]

# ---------------------------------------------------------------------------
# Distance / redshift
# ---------------------------------------------------------------------------


def z_from_lumi_dist_mpc(d_L_mpc: float) -> float:
    """Invert the Hubble-law d_L(z) = c*z/H0*(1+z) to find z from d_L in Mpc."""
    c_km_s = 299792.458
    H0 = 67.4
    x = d_L_mpc * H0 / c_km_s
    return float((-1 + np.sqrt(1 + 4 * x)) / 2)


# ---------------------------------------------------------------------------
# Frequency parsing
# ---------------------------------------------------------------------------


def parse_entry(value: str):
    """Parse frequency entry: Hz value, unit suffix, filter name, named band, or [lo,hi]."""
    if value in _NAMED_BANDS:
        nu_min, nu_max = _NAMED_BANDS[value]
        return (nu_min, nu_max, value)
    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1]
        parts = inner.split(",")
        if len(parts) == 2:
            lo = parse_entry(parts[0].strip())
            hi = parse_entry(parts[1].strip())
            if isinstance(lo, (int, float)) and isinstance(hi, (int, float)):
                return (float(lo), float(hi), value)
        raise ValueError(f"Bad band format: {value}")
    match = re.match(r"^([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*([a-zA-Z]+)$", value)
    if match:
        num, suffix = float(match.group(1)), match.group(2)
        if suffix in UNIT_SUFFIXES:
            return num * UNIT_SUFFIXES[suffix]
    return _parse_frequency(value)


_FREQ_SPLIT = re.compile(r",(?![^\[]*\])")


def parse_frequency_input(frequencies_input: str) -> tuple[list[float], list[list[Any]], list[str]]:
    frequencies: list[float] = []
    bands: list[list[Any]] = []
    warnings: list[str] = []
    tokens = [tok.strip() for tok in _FREQ_SPLIT.split(frequencies_input) if tok.strip()]

    for token in tokens:
        try:
            entry = parse_entry(token)
            if isinstance(entry, tuple):
                bands.append([float(entry[0]), float(entry[1]), str(entry[2])])
            else:
                frequencies.append(float(entry))
        except Exception:
            warnings.append(f"Unknown frequency or filter: '{token}'")

    return sorted(frequencies), bands, warnings


# ---------------------------------------------------------------------------
# Shared params -> physics dict
# ---------------------------------------------------------------------------


def shared_to_physics(shared: SharedParams) -> dict[str, Any]:
    """Convert SharedParams to a physics-only dict (exclude UI fields, add derived values)."""
    physics = shared.model_dump(exclude={"d_L_mpc", "num_t"})
    physics["d_L_cm"] = shared.d_L_mpc * 3.0856775814913673e24
    physics["z"] = z_from_lumi_dist_mpc(shared.d_L_mpc)
    return physics
