"""Host-galaxy dust extinction laws for the afterglow fitter.

Pei (1992) Table 4 fits for SMC, LMC, and MW. Returns ``k(lambda) = A(lambda)/A_V``,
renormalized so ``k(V) = k(5500 A) == 1`` exactly. Truncated to zero above the
Lyman limit (lambda < 912 A): the Pei92 fits are calibrated UV-IR only, the
Drude functional form has no K-edge structure, and X-ray opacity is dominated
by gas (N_H) not dust -- so any nonzero value the formula returns at EUV / X-ray
wavelengths is a numerical artifact of the BKG component's tail, not physics.
On the long-wavelength side Pei92 already self-truncates (``k ~ 1e-10`` at radio).
"""

import numpy as np

# Pei (1992) Table 4 -- six-component fits to xi(lambda)/A_B + R_V values.
# Components: background (0.04 um), FUV (0.08 um), 2175 A bump (0.22 um),
# and three FIR features (9.7, 18, 25 um).
_PEI92_COEFFS = {
    "smc": {
        "a": np.array([185.0, 27.0, 0.005, 0.010, 0.012, 0.030]),
        "lam_um": np.array([0.042, 0.08, 0.22, 9.7, 18.0, 25.0]),
        "b": np.array([90.0, 5.50, -1.95, -1.95, -1.80, 0.0]),
        "n": np.array([2.0, 4.0, 2.0, 2.0, 2.0, 2.0]),
        "R_V": 2.93,
    },
    "lmc": {
        "a": np.array([175.0, 19.0, 0.023, 0.005, 0.006, 0.020]),
        "lam_um": np.array([0.046, 0.08, 0.22, 9.7, 18.0, 25.0]),
        "b": np.array([90.0, 5.50, -1.95, -1.95, -1.80, 0.0]),
        "n": np.array([2.0, 4.5, 2.0, 2.0, 2.0, 2.0]),
        "R_V": 3.16,
    },
    "mw": {
        "a": np.array([165.0, 14.0, 0.045, 0.002, 0.002, 0.012]),
        "lam_um": np.array([0.047, 0.08, 0.22, 9.7, 18.0, 25.0]),
        "b": np.array([90.0, 4.0, -1.95, -1.95, -1.80, 0.0]),
        "n": np.array([2.0, 6.5, 2.0, 2.0, 2.0, 2.0]),
        "R_V": 3.08,
    },
}

_V_BAND_CM = 5.5e-5  # 5500 A -- fiducial V-band wavelength for renormalization
_LYMAN_LIMIT_CM = 9.12e-6  # 912 A -- ionization edge; Pei92 not valid below


def _pei92_raw(lam_cm, profile):
    """Raw Pei (1992) eq. 20 evaluation, before V-band renormalization.

    Zero above the Lyman limit (``lam_cm < 9.12e-6``).
    """
    c = _PEI92_COEFFS[profile]
    lam_arr = np.asarray(lam_cm, dtype=np.float64)
    lam_um = lam_arr * 1e4
    with np.errstate(divide="ignore", invalid="ignore"):
        x = lam_um[..., None] / c["lam_um"]
        denom = x ** c["n"] + (1.0 / x) ** c["n"] + c["b"]
        xi_over_AB = (c["a"] / denom).sum(axis=-1)
    raw = ((c["R_V"] + 1.0) / c["R_V"]) * xi_over_AB
    return np.where(lam_arr < _LYMAN_LIMIT_CM, 0.0, raw)


# Precomputed V-band values, so k(V) is identically 1 (Pei92 fit residual at V
# is ~3 percent for SMC; this absorbs it into the normalization).
_V_NORM = {p: float(_pei92_raw(_V_BAND_CM, p)) for p in _PEI92_COEFFS}


def pei92(lam_cm, profile="smc") -> np.ndarray:
    """k(lambda) = A(lambda)/A_V from Pei (1992) Table 4.

    Vectorized over ``lam_cm``; output has the same shape as the input.

    Parameters
    ----------
    lam_cm : float or array_like
        Rest-frame wavelength in cm.
    profile : {'smc', 'lmc', 'mw'}
        Pei (1992) extinction-curve profile.

    Returns
    -------
    ndarray
        k(lambda) = A(lambda)/A_V. ``k(V) == 1`` exactly. Returns 0 below the
        Lyman limit (912 A) -- see module docstring for why.
    """
    return _pei92_raw(lam_cm, profile) / _V_NORM[profile]


def smc(lam_cm) -> np.ndarray:
    """Pei (1992) SMC extinction law. See :func:`pei92`."""
    return pei92(lam_cm, "smc")


def lmc(lam_cm) -> np.ndarray:
    """Pei (1992) LMC extinction law. See :func:`pei92`."""
    return pei92(lam_cm, "lmc")


def mw(lam_cm) -> np.ndarray:
    """Pei (1992) Milky Way extinction law. See :func:`pei92`."""
    return pei92(lam_cm, "mw")


BUILTIN_LAWS = {"smc": smc, "lmc": lmc, "mw": mw}
