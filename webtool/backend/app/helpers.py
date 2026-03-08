"""Pure utility functions for parsing and conversions."""

from __future__ import annotations

import re

import numpy as np

from VegasAfterglow.cli import _parse_frequency
from VegasAfterglow.units import _NAMED_BANDS

from .constants import UNIT_SUFFIXES


def z_from_lumi_dist_mpc(d_L_mpc: float) -> float:
    """Invert the Hubble-law d_L(z) = c*z/H0*(1+z) to find z from d_L in Mpc."""
    c_km_s = 299792.458
    H0 = 67.4
    x = d_L_mpc * H0 / c_km_s
    return float((-1 + np.sqrt(1 + 4 * x)) / 2)


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


def mag_to_cgs(mag: float, mag_err: float):
    """Convert AB magnitude (+ error) to F_ν in CGS (erg/s/cm²/Hz)."""
    f_nu = 10.0 ** (-(mag + 48.6) / 2.5)
    f_err = f_nu * (np.log(10.0) / 2.5) * abs(mag_err) if mag_err else 0.0
    return f_nu, f_err
