"""LaTeX / human-readable formatting for frequencies, bands, and parameters."""

import numpy as np

from ..units import _NAMED_BANDS, _c_A
from ..units import keV as _keV
from .colors import _broad_band


def _format_energy(nu, latex=False):
    """Format a frequency as a human-readable energy/wavelength/frequency string.

    Uses wavelength (nm) for optical/IR/UV, photon energy for X-ray and above,
    and frequency (GHz/MHz/Hz) for radio.
    """
    E_keV = nu / _keV
    lam_nm = _c_A / nu / 10

    def _v(val, unit):
        s = f"{val:.3g}"
        return rf"${s}$ {unit}" if latex else f"{s} {unit}"

    if E_keV >= 1e9:
        return _v(E_keV / 1e9, "TeV")
    if E_keV >= 1e6:
        return _v(E_keV / 1e6, "GeV")
    if E_keV >= 1e3:
        return _v(E_keV / 1e3, "MeV")
    if E_keV >= 0.1:
        return _v(E_keV, "keV")
    if 100 < lam_nm < 10000:
        return _v(lam_nm, "nm")
    if nu >= 1e12:
        return _v(nu / 1e12, "THz")
    if nu >= 1e9:
        s = f"{nu / 1e9:g}"
        return rf"${s}$ GHz" if latex else f"{s} GHz"
    if nu >= 1e6:
        s = f"{nu / 1e6:.0f}"
        return rf"${s}$ MHz" if latex else f"{s} MHz"
    s = f"{nu:.0f}"
    return rf"${s}$ Hz" if latex else f"{s} Hz"


def _format_nu_latex(nu, label=None):
    """Format a frequency as a LaTeX label with broad-band prefix for plots."""
    from .. import units

    band = _broad_band(nu)

    # Show filter name only if the user actually typed a filter name
    if label is not None:
        try:
            float(label)
        except ValueError:
            if label in units._VEGA_FILTERS:
                return rf"{band} (${label}$-band)"
            if label in units._ST_FILTERS:
                return rf"{band} ({label})"
            if label in units._SURVEY_FILTERS:
                return rf"{band} ({label})"

    return rf"{band} ({_format_energy(nu, latex=True)})"


def _format_band(nu_min, nu_max, name=None, latex=False):
    """Format a frequency band label (plain text or LaTeX)."""
    lo = _format_energy(nu_min, latex=latex)
    hi = _format_energy(nu_max, latex=latex)
    sep = "–" if latex else "-"
    range_str = f"{lo}{sep}{hi}"
    if name and name in _NAMED_BANDS:
        return f"{name} ({range_str})" if latex else f"{name}({range_str})"
    return range_str


def _sci_tex(val):
    """Format a number as LaTeX scientific notation: 1e52 -> 10^{52}."""
    if val == 0:
        return "0"
    exp = int(np.floor(np.log10(abs(val))))
    coeff = val / 10**exp
    if abs(coeff - 1) < 0.01:
        return rf"10^{{{exp}}}"
    return rf"{coeff:.1f}\times10^{{{exp}}}"
