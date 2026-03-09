"""Shared constants for the VegasAfterglow webapp."""

from VegasAfterglow import units

UNIT_SUFFIXES = {
    "Hz": 1, "kHz": units.kHz, "MHz": units.MHz, "GHz": units.GHz,
    "eV": units.eV, "keV": units.keV, "MeV": units.MeV, "GeV": units.GeV,
}

OBS_FLUX_UNITS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm\u00b2/s"]
