"""Shared constants for the VegasAfterglow webapp."""

from VegasAfterglow import units
from VegasAfterglow.cli import _FLUX_SCALES, _TIME_SCALES

UNIT_SUFFIXES = {
    "Hz": 1, "kHz": units.kHz, "MHz": units.MHz, "GHz": units.GHz,
    "eV": units.eV, "keV": units.keV, "MeV": units.MeV, "GeV": units.GeV,
}

FREQ_SCALES = {"Hz": 1, "GHz": units.GHz, "keV": units.keV, "MeV": units.MeV}

# Instrument sensitivities: (nu_min_Hz, nu_max_Hz, sensitivity, kind)
# kind="Fnu": flux density (erg/cm²/s/Hz) — dashed line on primary y-axis
# kind="Fband": band-integrated flux (erg/cm²/s) — solid line on secondary y-axis
INSTRUMENTS = {
    # -- Radio (Fnu, ~1h integration) --
    "VLA":           (1e9,     5e10,    5e-29,   "Fnu"),
    "ALMA":          (8.4e10,  9.5e11,  3e-28,   "Fnu"),
    "MeerKAT":       (9e8,     1.67e9,  5e-29,   "Fnu"),
    "ngVLA":         (1e9,     1.16e11, 3e-30,   "Fnu"),
    # -- Optical/IR (Fnu) --
    "Rubin/LSST":    (3.3e14,  1e15,    3.6e-31, "Fnu"),
    "JWST":          (6e13,    2.5e14,  9e-33,   "Fnu"),
    "WFST":          (3e14,    9.4e14,  3.6e-30, "Fnu"),
    "SVOM/VT":       (4.6e14,  7.5e14,  3.6e-30, "Fnu"),
    # -- X-ray (Fband, erg/cm²/s) --
    "Swift/XRT":     (7.25e16, 2.42e18, 1e-13,   "Fband"),
    "Chandra":       (1.21e17, 1.93e18, 1e-15,   "Fband"),
    "EP/WXT":        (1.21e17, 9.67e17, 1e-10,   "Fband"),
    "EP/FXT":        (7.25e16, 2.42e18, 2e-14,   "Fband"),
    "SVOM/MXT":      (4.84e16, 2.42e18, 1e-12,   "Fband"),
    "SVOM/ECLAIRs":  (9.67e17, 3.63e19, 5e-10,   "Fband"),
    # Swift BAT 5-sigma ~2e-8 erg/cm2/s (15-150 keV, 1 s trigger sensitivity).
    "Swift/BAT":     (3.63e18, 3.63e19, 2e-8,    "Fband"),
    # -- Gamma-ray (Fband, erg/cm²/s) --
    # Fermi GBM onboard threshold ~0.7 ph/cm2/s (50-300 keV, 1 s), converted
    # to an approximate energy flux scale ~2e-7 erg/cm2/s.
    "Fermi/GBM":     (1.93e18, 9.67e21, 2e-7,    "Fband"),
    "Fermi/LAT":     (4.84e21, 7.25e25, 1e-12,   "Fband"),
    "CTA":           (4.84e24, 7.25e28, 1e-13,   "Fband"),
}

OBS_FLUX_UNITS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm\u00b2/s"]

# Re-export cli constants for convenience
FLUX_SCALES = _FLUX_SCALES
TIME_SCALES = _TIME_SCALES
