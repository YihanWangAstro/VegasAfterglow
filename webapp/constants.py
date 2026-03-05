"""Shared constants for the VegasAfterglow webapp."""

from VegasAfterglow import units
from VegasAfterglow.cli import (
    _COMP_LABELS,
    _FLUX_SCALES,
    _TIME_LABELS,
    _TIME_SCALES,
)

# Plotly dash equivalents for matplotlib line styles
PLOTLY_DASH = {
    "total": "solid",
    "fwd_sync": "dash",
    "fwd_ssc": "dot",
    "rvs_sync": "dashdot",
    "rvs_ssc": "longdashdot",
}

UNIT_SUFFIXES = {
    "Hz": 1, "kHz": units.kHz, "MHz": units.MHz, "GHz": units.GHz,
    "eV": units.eV, "keV": units.keV, "MeV": units.MeV, "GeV": units.GeV,
}

FREQ_SCALES = {"Hz": 1, "GHz": units.GHz, "keV": units.keV, "MeV": units.MeV}

# Plotly-compatible flux labels (HTML tags for Plotly axis titles)
PLOTLY_FLUX_LABELS = {
    "mJy": "mJy",
    "Jy": "Jy",
    "uJy": "\u03bcJy",
    "cgs": "erg cm<sup>\u22122</sup> s<sup>\u22121</sup> Hz<sup>\u22121</sup>",
    "AB mag": "AB mag",
}

# Instrument sensitivities: (nu_min_Hz, nu_max_Hz, sensitivity, kind, color)
# kind="Fnu": flux density (erg/cm²/s/Hz) — dashed line on primary y-axis
# kind="Fband": band-integrated flux (erg/cm²/s) — solid line on secondary y-axis
INSTRUMENTS = {
    # -- Radio (Fnu, ~1h integration) --
    "VLA":           (1e9,     5e10,    5e-29,   "Fnu",   "#E63946"),
    "ALMA":          (8.4e10,  9.5e11,  3e-28,   "Fnu",   "#457B9D"),
    "MeerKAT":       (9e8,     1.67e9,  5e-29,   "Fnu",   "#2A9D8F"),
    "ngVLA":         (1e9,     1.16e11, 3e-30,   "Fnu",   "#E9C46A"),
    # -- Optical/IR (Fnu) --
    "Rubin/LSST":    (3.3e14,  1e15,    3.6e-31, "Fnu",   "#F4A261"),
    "JWST":          (6e13,    2.5e14,  9e-33,   "Fnu",   "#264653"),
    "WFST":          (3e14,    9.4e14,  3.6e-30, "Fnu",   "#8B5E3C"),
    "SVOM/VT":       (4.6e14,  7.5e14,  3.6e-30, "Fnu",   "#D4A017"),
    # -- X-ray (Fband, erg/cm²/s) --
    "Swift/XRT":     (7.25e16, 2.42e18, 1e-13,   "Fband", "#A8DADC"),
    "Chandra":       (1.21e17, 1.93e18, 1e-15,   "Fband", "#1D3557"),
    "EP/WXT":        (1.21e17, 9.67e17, 1e-10,   "Fband", "#FF6B6B"),
    "EP/FXT":        (7.25e16, 2.42e18, 2e-14,   "Fband", "#4ECDC4"),
    "SVOM/MXT":      (4.84e16, 2.42e18, 1e-12,   "Fband", "#45B7D1"),
    "SVOM/ECLAIRs":  (9.67e17, 3.63e19, 5e-10,   "Fband", "#96CEB4"),
    # -- Gamma-ray (Fband, erg/cm²/s) --
    "Fermi/LAT":     (4.84e21, 7.25e25, 1e-12,   "Fband", "#6A0572"),
    "CTA":           (4.84e24, 7.25e28, 1e-13,   "Fband", "#C44536"),
}

FBAND_TITLE = "F (erg cm<sup>\u22122</sup> s<sup>\u22121</sup>)"

TIME_PALETTE = ["#E03530", "#E8872E", "#2AB07E", "#2878B5", "#7B3FA0", "#D4A017"]

OBS_COLORS = [
    "#000000", "#E03530", "#2878B5", "#2AB07E", "#7B3FA0",
    "#E8872E", "#D4A017", "#555555", "#C44536", "#1D3557",
]

OBS_FLUX_UNITS = ["mJy", "Jy", "uJy", "cgs", "AB mag", "erg/cm\u00b2/s"]

# Publication-quality axis style (matches CLI _setup_plot_style)
AXIS_COMMON = dict(
    showline=True, linewidth=0.8, linecolor="#000", mirror=True,
    ticks="inside", ticklen=5, tickwidth=0.8, tickcolor="#000",
    minor=dict(ticks="inside", ticklen=2.5, showgrid=True,
               gridcolor="rgba(0,0,0,0.06)", griddash="dot", gridwidth=0.3),
    showgrid=True, gridcolor="rgba(0,0,0,0.10)", griddash="dot", gridwidth=0.3,
    title_font=dict(size=13, color="#000"), tickfont=dict(size=11, color="#000"),
    exponentformat="power",
)

LEGEND_COMMON = dict(
    orientation="v",
    yanchor="bottom",
    y=0.02,
    xanchor="left",
    x=0.02,
    font=dict(size=12, color="#333"),
    tracegroupgap=1,
    borderwidth=0,
    bgcolor="rgba(255,255,255,0.75)",
)

LAYOUT_COMMON = dict(
    template="plotly_white",
    hovermode="closest",
    height=465,
    width=620,
    margin=dict(l=65, r=20, t=15, b=55),
    plot_bgcolor="#ffffff",
    paper_bgcolor="#ffffff",
)

# Re-export cli constants for convenience
COMP_LABELS = _COMP_LABELS
FLUX_SCALES = _FLUX_SCALES
TIME_SCALES = _TIME_SCALES
TIME_LABELS = _TIME_LABELS
