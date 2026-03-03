"""Pure utility functions: parsing, formatting, conversions."""

import re

import numpy as np
import streamlit as st

from VegasAfterglow.cli import (
    _broad_band,
    _format_band,
    _format_energy,
    _freq_colors,
    _parse_frequency,
)
from VegasAfterglow.units import _NAMED_BANDS

from .constants import TIME_PALETTE, UNIT_SUFFIXES


def z_from_lumi_dist_mpc(d_L_mpc):
    """Invert the Hubble-law d_L(z) = c*z/H0*(1+z) to find z from d_L in Mpc."""
    c_km_s = 299792.458
    H0 = 67.4  # km/s/Mpc  (must match _lumi_dist_from_z)
    x = d_L_mpc * H0 / c_km_s  # = z + z^2
    return (-1 + np.sqrt(1 + 4 * x)) / 2


def log_slider(label, min_exp, max_exp, default_exp, step=0.1, key=None):
    """Slider in exponent space. Label shows log10(name). Returns 10^x."""
    exp = st.slider(
        f"log10({label})",
        min_value=float(min_exp),
        max_value=float(max_exp),
        value=float(default_exp),
        step=float(step),
        key=key,
        format="%.1f",
    )
    return 10**exp


def smart_ylim(flux_arrays, scale=1.0):
    """Compute smart log y-axis limits from flux arrays."""
    all_pos = []
    for arr in flux_arrays:
        scaled = arr / scale
        pos = scaled[scaled > 0]
        if pos.size > 0:
            all_pos.append(pos)
    if not all_pos:
        return None, None
    combined = np.concatenate(all_pos)
    f_max, f_min = np.max(combined), np.min(combined)
    y_top = f_max * 10
    y_bot = f_min / 10
    if y_top / y_bot > 1e12:
        y_bot = y_top * 1e-12
    return y_bot, y_top


def parse_entry(value):
    """Parse a frequency entry: float Hz, unit suffix, filter name, named band, or [lo,hi]."""
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
    m = re.match(r"^([0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*([a-zA-Z]+)$", value)
    if m:
        num, suffix = float(m.group(1)), m.group(2)
        if suffix in UNIT_SUFFIXES:
            return num * UNIT_SUFFIXES[suffix]
    return _parse_frequency(value)


def freq_label(nu):
    """Human-readable frequency label for plot legend (plain text, no LaTeX)."""
    band = _broad_band(nu).replace(r"$\gamma$", "\u03b3")
    return f"{band} ({_format_energy(nu)})"


def band_label(nu_min, nu_max, name=None):
    """Human-readable band label for plot legend (plain text)."""
    return _format_band(nu_min, nu_max, name, latex=False)


def time_colors(times):
    """Map observer time snapshots to distinct colors."""
    n = len(times)
    return [TIME_PALETTE[i % len(TIME_PALETTE)] for i in range(n)]


def format_time_label(t_sec):
    """Format observer time for SED legend (plain text)."""
    if t_sec >= 365.25 * 86400:
        return f"{t_sec / (365.25 * 86400):.3g} yr"
    if t_sec >= 86400:
        return f"{t_sec / 86400:.3g} day"
    if t_sec >= 3600:
        return f"{t_sec / 3600:.3g} hr"
    if t_sec >= 60:
        return f"{t_sec / 60:.3g} min"
    return f"{t_sec:.3g} s"




def mag_to_cgs(mag, mag_err):
    """Convert AB magnitude (+ error) to F_ν in CGS (erg/s/cm²/Hz)."""
    f_nu = 10.0 ** (-(mag + 48.6) / 2.5)
    f_err = f_nu * (np.log(10.0) / 2.5) * abs(mag_err) if mag_err else 0.0
    return f_nu, f_err
