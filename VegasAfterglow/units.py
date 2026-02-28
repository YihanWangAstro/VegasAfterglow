"""Unit conversion constants and magnitude converters for VegasAfterglow.

Each constant equals "1 of that unit in the API's base unit", so multiplying
raw data by a constant converts it to the base unit the API expects.

Base units: seconds, Hz, erg/cm^2/s/Hz, cm, radians.

Example::

    from VegasAfterglow.units import day, GHz, mJy, Mpc, deg

    # Input: multiply to convert TO base units
    obs = Observer(lumi_dist=100*Mpc, z=0.5, theta_obs=5*deg)
    flux = model.flux_density_grid(times * day, nus * GHz)

    # Output: divide to convert FROM base units
    flux_mJy = flux.total / mJy
    t_days = times / day

Magnitude converters::

    from VegasAfterglow.units import ABmag_to_cgs, Vegamag_to_cgs, filter_freq

    f_nu = ABmag_to_cgs(mag_array)          # AB magnitudes
    f_nu = Vegamag_to_cgs(mag_array, "R")   # Vega magnitudes in R band
    nu = filter("R")                         # effective frequency [Hz]
"""

import math

import numpy as np

# ---------------------------------------------------------------------------
# Time (base: seconds)
# ---------------------------------------------------------------------------
sec = 1.0
ms = 1e-3
minute = 60.0
hr = 3600.0
day = 86400.0
yr = 365.2425 * day

# ---------------------------------------------------------------------------
# Frequency (base: Hz)
# ---------------------------------------------------------------------------
Hz = 1.0
kHz = 1e3
MHz = 1e6
GHz = 1e9

# ---------------------------------------------------------------------------
# Photon energy -> frequency via E = h*nu (base: Hz)
# ---------------------------------------------------------------------------
_h_cgs = 6.62607015e-27  # Planck constant [erg·s]
_eV_cgs = 1.602176634e-12  # 1 eV [erg]

eV = _eV_cgs / _h_cgs
keV = 1e3 * eV
MeV = 1e6 * eV
GeV = 1e9 * eV

# ---------------------------------------------------------------------------
# Flux density (base: erg/cm^2/s/Hz)
# ---------------------------------------------------------------------------
Jy = 1e-23
mJy = 1e-26
uJy = 1e-29

# ---------------------------------------------------------------------------
# Distance (base: cm)
# ---------------------------------------------------------------------------
cm = 1.0
m = 1e2
km = 1e5
pc = 3.0856775814913673e18
kpc = 1e3 * pc
Mpc = 1e6 * pc
Gpc = 1e9 * pc

# ---------------------------------------------------------------------------
# Angle (base: radians)
# ---------------------------------------------------------------------------
rad = 1.0
deg = math.pi / 180
arcmin = deg / 60
arcsec = deg / 3600

# ===========================================================================
# Magnitude system converters
# ===========================================================================

_c_A = 2.99792458e18  # speed of light [Å/s]
_AB_F0 = 3.631e-20  # AB zero-point flux density [erg/cm²/s/Hz] (= 3631 Jy)
_ST_F0 = 3.631e-9  # ST zero-point flux density [erg/cm²/s/Å]

# ---------------------------------------------------------------------------
# Vega filter data: (zero-point flux [Jy], effective wavelength [Å])
#   Johnson-Cousins: Bessell & Murphy (2012)
#   2MASS: Cohen et al. (2003)
#   Swift UVOT: Breeveld et al. (2011), Poole et al. (2008)
#   SDSS: Fukugita et al. (1996); AB offsets from Blanton & Roweis (2007)
# ---------------------------------------------------------------------------
_VEGA_FILTERS = {
    # Johnson-Cousins
    "U": (1810.0, 3663.0),
    "B": (4130.0, 4361.0),
    "V": (3640.0, 5448.0),
    "R": (3080.0, 6407.0),
    "I": (2550.0, 7980.0),
    # 2MASS
    "J": (1594.0, 12350.0),
    "H": (1024.0, 16620.0),
    "Ks": (666.7, 21590.0),
    # Swift UVOT (lowercase convention)
    "v": (3665.0, 5468.0),
    "b": (4092.0, 4392.0),
    "u": (1420.0, 3465.0),
    "uvw1": (905.0, 2600.0),
    "uvm2": (766.0, 2246.0),
    "uvw2": (741.0, 1928.0),
    # SDSS (also used by WFST and Mephisto)
    "g": (3910.0, 4686.0),
    "r": (3130.0, 6166.0),
    "i": (2580.0, 7480.0),
    "z": (2210.0, 8932.0),
}

# ---------------------------------------------------------------------------
# ST (HST) filter data: pivot wavelength [Å]
#   WFC3/UVIS and WFC3/IR from instrument handbooks
# ---------------------------------------------------------------------------
_ST_FILTERS = {
    "F225W": 2372.0,
    "F275W": 2710.0,
    "F336W": 3355.0,
    "F438W": 4326.0,
    "F475W": 4773.0,
    "F555W": 5308.0,
    "F606W": 5889.0,
    "F625W": 6243.0,
    "F775W": 7651.0,
    "F814W": 8039.0,
    "F850LP": 9176.0,
    "F105W": 10552.0,
    "F110W": 11534.0,
    "F125W": 12486.0,
    "F140W": 13923.0,
    "F160W": 15369.0,
}

# ---------------------------------------------------------------------------
# Survey / telescope-specific filters: effective wavelength [Å]
#   SVOM VT: Götz et al. (2024)
#   WFST: Lei et al. (2023)
# ---------------------------------------------------------------------------
_SURVEY_FILTERS = {
    # SVOM Visible Telescope dichroic channels
    "VT_B": 5500.0,  # 450–650 nm blue channel
    "VT_R": 8200.0,  # 650–1000 nm red channel
    # WFST wide-band filter
    "w": 6200.0,  # ~390–820 nm wide filter
}


# ---------------------------------------------------------------------------
# AB magnitude (universal zero-point, no filter needed)
# ---------------------------------------------------------------------------


def ABmag_to_cgs(mag):
    """Convert AB magnitude(s) to flux density in erg/cm²/s/Hz.

    Parameters
    ----------
    mag : float or array_like
        AB magnitude(s).

    Returns
    -------
    float or ndarray
        Flux density in erg/cm²/s/Hz.
    """
    return _AB_F0 * np.power(10.0, -0.4 * np.asarray(mag, dtype=float))


def cgs_to_ABmag(f_nu):
    """Convert flux density in erg/cm²/s/Hz to AB magnitude(s).

    Parameters
    ----------
    f_nu : float or array_like
        Flux density in erg/cm²/s/Hz.

    Returns
    -------
    float or ndarray
        AB magnitude(s).
    """
    return -2.5 * np.log10(np.asarray(f_nu, dtype=float) / _AB_F0)


# ---------------------------------------------------------------------------
# Vega magnitude (per-filter zero-point)
# ---------------------------------------------------------------------------


def _get_filter(filt, table, system):
    if filt not in table:
        avail = ", ".join(sorted(table))
        raise ValueError(f"Unknown {system} filter '{filt}'. Available: {avail}")
    return table[filt]


def Vegamag_to_cgs(mag, filt):
    """Convert Vega magnitude(s) to flux density in erg/cm²/s/Hz.

    Parameters
    ----------
    mag : float or array_like
        Vega magnitude(s).
    filt : str
        Filter name (e.g. 'V', 'J', 'Ks', 'uvw1').

    Returns
    -------
    float or ndarray
        Flux density in erg/cm²/s/Hz.
    """
    zp_Jy, _ = _get_filter(filt, _VEGA_FILTERS, "Vega")
    return zp_Jy * 1e-23 * np.power(10.0, -0.4 * np.asarray(mag, dtype=float))


def cgs_to_Vegamag(f_nu, filt):
    """Convert flux density in erg/cm²/s/Hz to Vega magnitude(s).

    Parameters
    ----------
    f_nu : float or array_like
        Flux density in erg/cm²/s/Hz.
    filt : str
        Filter name (e.g. 'V', 'J', 'Ks', 'uvw1').

    Returns
    -------
    float or ndarray
        Vega magnitude(s).
    """
    zp_Jy, _ = _get_filter(filt, _VEGA_FILTERS, "Vega")
    return -2.5 * np.log10(np.asarray(f_nu, dtype=float) / (zp_Jy * 1e-23))


# ---------------------------------------------------------------------------
# ST magnitude (HST STMAG system, per-filter pivot wavelength)
# ---------------------------------------------------------------------------


def STmag_to_cgs(mag, filt):
    """Convert ST magnitude(s) to flux density f_nu in erg/cm²/s/Hz.

    Internally computes f_lambda from the ST zero-point, then converts to
    f_nu using the filter's pivot wavelength: f_nu = f_lambda * lambda^2 / c.

    Parameters
    ----------
    mag : float or array_like
        ST magnitude(s).
    filt : str
        HST filter name (e.g. 'F606W', 'F160W').

    Returns
    -------
    float or ndarray
        Flux density in erg/cm²/s/Hz.
    """
    lam_A = _get_filter(filt, _ST_FILTERS, "ST")
    f_lambda = _ST_F0 * np.power(10.0, -0.4 * np.asarray(mag, dtype=float))
    return f_lambda * lam_A**2 / _c_A


def cgs_to_STmag(f_nu, filt):
    """Convert flux density f_nu in erg/cm²/s/Hz to ST magnitude(s).

    Parameters
    ----------
    f_nu : float or array_like
        Flux density in erg/cm²/s/Hz.
    filt : str
        HST filter name (e.g. 'F606W', 'F160W').

    Returns
    -------
    float or ndarray
        ST magnitude(s).
    """
    lam_A = _get_filter(filt, _ST_FILTERS, "ST")
    f_lambda = np.asarray(f_nu, dtype=float) * _c_A / lam_A**2
    return -2.5 * np.log10(f_lambda / _ST_F0)


# ---------------------------------------------------------------------------
# Filter effective frequency lookup
# ---------------------------------------------------------------------------


_ALL_FILTER_WAVELENGTHS = {}
for _name, (_zp, _lam) in _VEGA_FILTERS.items():
    _ALL_FILTER_WAVELENGTHS[_name] = _lam
for _name, _lam in {**_ST_FILTERS, **_SURVEY_FILTERS}.items():
    _ALL_FILTER_WAVELENGTHS[_name] = _lam


def filter(filt):
    """Return the effective frequency [Hz] of a named filter.

    Works for Vega-system, ST-system, and survey telescope filters.

    Parameters
    ----------
    filt : str
        Filter name (e.g. 'V', 'J', 'F606W', 'g', 'VT_B').

    Returns
    -------
    float
        Effective frequency in Hz.
    """
    if filt in _ALL_FILTER_WAVELENGTHS:
        return _c_A / _ALL_FILTER_WAVELENGTHS[filt]
    avail = ", ".join(sorted(_ALL_FILTER_WAVELENGTHS))
    raise ValueError(f"Unknown filter '{filt}'. Available: {avail}")


# ---------------------------------------------------------------------------
# Named instrument bands: (nu_min_Hz, nu_max_Hz)
# ---------------------------------------------------------------------------

_NAMED_BANDS = {
    "XRT": (0.3 * keV, 10 * keV),  # Swift X-Ray Telescope
    "BAT": (15 * keV, 150 * keV),  # Swift Burst Alert Telescope
    "FXT": (0.3 * keV, 10 * keV),  # Einstein Probe FXT
    "WXT": (0.5 * keV, 4 * keV),  # Einstein Probe WXT
    "MXT": (0.2 * keV, 10 * keV),  # SVOM MXT
    "ECLAIRs": (4 * keV, 150 * keV),  # SVOM ECLAIRs
    "LAT": (100 * MeV, 300 * GeV),  # Fermi LAT
    "GBM": (8 * keV, 40 * MeV),  # Fermi GBM
}


def band(name):
    """Return (nu_min, nu_max) in Hz for a named instrument band.

    Parameters
    ----------
    name : str
        Band name (e.g. 'XRT', 'LAT', 'BAT').

    Returns
    -------
    tuple of float
        (nu_min, nu_max) in Hz.
    """
    if name not in _NAMED_BANDS:
        avail = ", ".join(sorted(_NAMED_BANDS))
        raise ValueError(f"Unknown band '{name}'. Available: {avail}")
    return _NAMED_BANDS[name]
