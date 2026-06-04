"""Color registry and frequency/filter color helpers.

Used by both the ``vegasgen`` CLI plots and ``Fitter.draw_best_fit``. Three
layers stacked from most-specific to least-specific:

1. ``_FILTER_COLOR_MAP`` — canonical colors for known filter / instrument
   names (Johnson, SDSS, 2MASS, HST WFC3, Swift, EP, Fermi, …).
2. ``_emission_spectrum_color(nu)`` — frequency-based EM-spectrum gradient.
3. ``FREQ_PALETTE`` / ``_freq_colors(nus)`` — discrete qualitative palette
   for the CLI when no filter names are involved.

``_broad_band(nu)`` lives here too because the legend / label code treats
"Radio / Optical / X-ray / …" as a broadband classification on top of the
same frequency table.
"""

import numpy as np

# Discrete qualitative palette ordered warm → cool (radio → X-ray/gamma).
FREQ_PALETTE = [
    "#E03530",  # red
    "#E8872E",  # orange
    "#D4A017",  # gold
    "#8C564B",  # brown
    "#BCBD22",  # olive
    "#2AB07E",  # green
    "#17BECF",  # cyan
    "#2878B5",  # blue
    "#1D3557",  # navy
    "#7B3FA0",  # purple
    "#E377C2",  # pink
    "#555555",  # gray
]


def _freq_colors(nus):
    """Map frequencies to distinct colors ordered warm (low ν) → cool (high ν).

    Assigns colors from a discrete qualitative palette by frequency rank,
    spread evenly across the full palette so even a few frequencies use
    the full warm→cool range.
    """
    n = len(nus)
    if n == 0:
        return []
    P = len(FREQ_PALETTE)
    order = np.argsort(nus)
    colors = [""] * n
    step = (P - 1) / max(1, n - 1) if n <= P else 1
    for rank, idx in enumerate(order):
        ci = round(rank * step) if n <= P else rank % P
        colors[idx] = FREQ_PALETTE[ci]
    return colors


# Canonical colors for well-known photometric filters and X-ray / gamma
# instruments. Visible filters use the rough hue of their bandpass; non-visible
# bands use a stylised EM-spectrum gradient (UV → violet, X-ray → navy,
# γ-ray → very dark purple, IR → dark red / brown, radio → maroon).
_FILTER_COLOR_MAP = {
    # SDSS / Pan-STARRS / ZTF / DECam ugrizy
    "u": "#7B1FA2",  # near-UV
    "g": "#3F8F3F",  # green
    "r": "#D32F2F",  # red
    "i": "#8B1A1A",  # dark red
    "z": "#BFA28F",  # NIR light tan
    "y": "#9A7B5E",  # NIR medium tan
    # Johnson-Cousins UBVRI
    "U": "#7B1FA2",
    "B": "#1976D2",  # blue
    "V": "#3F8F3F",  # green
    "R": "#D32F2F",
    "I": "#8B1A1A",
    # 2MASS / UKIDSS near-IR: warm-brown ramp tuned for ΔE > 13 between any
    # adjacent pair, so H (mid-NIR) reads as clearly brown rather than the
    # red-brown that previously collided with the visible-red filters (i, R).
    "J": "#704435",  # medium warm brown
    "H": "#4A2C24",  # dark warm brown
    "K": "#2A1F1B",  # near-black brown
    "Ks": "#2A1F1B",
    # HST common WFC3/ACS broad filters
    "F275W": "#9C27B0",  # near-UV
    "F336W": "#7B1FA2",
    "F438W": "#1976D2",
    "F555W": "#3F8F3F",
    "F606W": "#7CB342",  # broad-V (lighter green to differ from g/V)
    "F814W": "#8B1A1A",
    "F125W": "#A0521C",
    "F160W": "#7C3A14",
    # X-ray and γ-ray instruments — spread across blue → cyan → indigo →
    # purple so several instruments can coexist on one plot without colour
    # clashes. Lightness still trends darker with photon energy within a hue
    # family, but the family changes when ΔE would otherwise drop below ~10.
    "WXT": "#4FC3F7",  # Einstein Probe WXT (soft, 0.5–4 keV) — bright cyan-blue
    "NICER": "#2196F3",  # NICER (0.2–12 keV) — medium blue
    "XRT": "#0D47A1",  # Swift XRT (0.3–10 keV) — deep blue
    "FXT": "#00838F",  # Einstein Probe FXT — teal
    "Chandra": "#1A237E",  # Chandra (0.5–7 keV) — dark indigo
    "BAT": "#5E35B1",  # Swift BAT (15–150 keV) — purple
    "GBM": "#311B92",  # Fermi GBM (8 keV – 40 MeV) — deep indigo-purple
    "LAT": "#0D0029",  # Fermi LAT (100 MeV+) — near-black purple
}


def _broad_band(nu):
    """Return the broad-band name (Radio / IR / Optical / …) for a frequency."""
    from ..units import eV as _eV_Hz

    E_eV = nu / _eV_Hz

    if nu < 1e12:  # < 1 THz
        return "Radio"
    if nu < 4e14:  # 1 THz - 400 THz  (λ > 750 nm)
        return "IR"
    if nu < 7.5e14:  # 400 - 750 THz  (400-750 nm)
        return "Optical"
    if E_eV < 100:  # 750 THz - 0.1 keV  (λ > 12 nm)
        return "UV"
    if E_eV < 1e5:  # 0.1 - 100 keV
        return "X-ray"
    if E_eV < 1e9:  # 100 keV - 1 GeV
        return r"$\gamma$-ray"
    if E_eV < 1e12:  # 1 GeV - 1 TeV
        return "GeV"
    return "TeV"


def _emission_spectrum_color(nu):
    """Frequency-based fallback color for filters not in ``_FILTER_COLOR_MAP``.

    Coarse mapping of nu (Hz) to a color along the EM spectrum:
    radio → maroon, IR → red-brown, visible → matching hue (red→violet),
    UV → violet, X-ray → navy, γ-ray → near-black.
    """
    if not np.isfinite(nu) or nu <= 0:
        return "#555555"
    if nu < 1e10:
        return "#3E1F0A"  # very long radio
    if nu < 1e12:
        return "#5C2A12"  # radio / sub-mm
    if nu < 1e14:
        return "#8B4513"  # far / mid IR
    if nu < 3e14:
        return "#A04A1B"  # near IR
    if nu < 4e14:
        return "#C04A14"  # red end of visible
    if nu < 5e14:
        return "#D32F2F"  # red
    if nu < 5.7e14:
        return "#E68A1A"  # orange
    if nu < 6.0e14:
        return "#D4AC0D"  # yellow
    if nu < 6.5e14:
        return "#3F8F3F"  # green
    if nu < 7.2e14:
        return "#1976D2"  # blue
    if nu < 3e15:
        return "#7B1FA2"  # UV
    if nu < 3e16:
        return "#5A1A8C"  # FUV / EUV
    if nu < 3e18:
        return "#1B3E72"  # soft X-ray
    if nu < 3e19:
        return "#0D2E5C"  # hard X-ray
    return "#1A0F4A"  # gamma-ray


def _auto_filter_name(nu, tol=0.01):
    """Closest canonical filter name within ``tol`` relative frequency, or None.

    Matches ``nu`` against the photometric filter table in ``units.py``
    (Johnson-Cousins, SDSS, 2MASS, Swift UVOT, HST WFC3, SVOM VT, WFST).
    Used by ``Fitter.draw_best_fit`` as a fallback when the user adds data via
    ``nu=filter('r')`` without an explicit ``label=`` — a 1% tolerance is
    tight enough to distinguish neighbouring filters (e.g. SDSS r vs Johnson
    R, ~4% apart) but loose enough to forgive small floating-point drift.
    """
    if not np.isfinite(nu) or nu <= 0:
        return None
    from ..units import _ALL_FILTER_WAVELENGTHS, _c_A

    log_nu = np.log10(nu)
    best_name, best_dist = None, float("inf")
    for name, lam in _ALL_FILTER_WAVELENGTHS.items():
        dist = abs(np.log10(_c_A / lam) - log_nu)
        if dist < best_dist:
            best_name, best_dist = name, dist
    return best_name if best_dist < np.log10(1.0 + tol) else None


def _filter_color(name, nu):
    """Return a stable color for a filter, preferring the name registry.

    Lookup order:
      1. Exact match in ``_FILTER_COLOR_MAP``.
      2. Suffix after the last ``_`` or ``-`` (handles e.g. ``VT_R`` → ``R``,
         ``EP_WXT`` → ``WXT``, ``ZTF-r`` → ``r``).
      3. ``_emission_spectrum_color(nu)`` (frequency-based fallback).
    """
    if name:
        key = str(name).strip()
        if key in _FILTER_COLOR_MAP:
            return _FILTER_COLOR_MAP[key]
        for sep in ("_", "-"):
            if sep in key:
                tail = key.rsplit(sep, 1)[1]
                if tail in _FILTER_COLOR_MAP:
                    return _FILTER_COLOR_MAP[tail]
    return _emission_spectrum_color(nu)


def _disambiguate_filter_colors(items):
    """Differentiate same-color filters at different frequencies via lightness.

    ``items`` is a sequence of ``(key, base_color, freq)`` tuples. Returns a
    dict ``{key: hex_color}``. Entries sharing the same ``base_color`` are
    ordered by ``freq`` and spread in HLS lightness so the lowest-frequency
    one stays darkest and the highest-frequency one becomes lighter — keeps
    each filter family recognisable while breaking ties.
    """
    import colorsys

    from matplotlib.colors import to_hex, to_rgb

    by_color = {}
    for key, color, freq in items:
        by_color.setdefault(color, []).append((freq, key))

    out = {}
    for color, group in by_color.items():
        if len(group) == 1:
            out[group[0][1]] = color
            continue
        group.sort(key=lambda x: x[0])  # ascending frequency
        r, g, b = to_rgb(color)
        h, l, s = colorsys.rgb_to_hls(r, g, b)
        n = len(group)
        for i, (_freq, key) in enumerate(group):
            f = i / (n - 1) if n > 1 else 0.5  # 0..1
            # ± 0.2 lightness spread around the base, clamped sensibly
            new_l = max(0.15, min(0.85, l + (f - 0.5) * 0.4))
            r2, g2, b2 = colorsys.hls_to_rgb(h, new_l, s)
            out[key] = to_hex((r2, g2, b2))
    return out
