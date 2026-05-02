"""Host-galaxy dust extinction laws for GRB afterglow MCMC fitting.

Implements the Pei (1992) six-component extinction curves for the SMC, LMC,
and Milky Way environments.  All functions return the dimensionless attenuation
``k(λ) = A(λ)/A_V`` normalised so that ``k(V) = 1`` exactly (V-band reference
wavelength 0.55 μm).  Wavelengths below the Lyman limit (912 Å) are set to
zero, so that X-ray data passed to the fitter is not spuriously attenuated.

The ``BUILTIN_LAWS`` mapping is used by :class:`~VegasAfterglow.Fitter` when
``extinction="smc"`` / ``"lmc"`` / ``"mw"`` is requested.

Usage
-----
Direct evaluation::

    import numpy as np
    from VegasAfterglow.extinction import smc, lmc, mw, pei92

    lam_cm = np.logspace(-8, -3, 200)   # 1 Å … 10 μm
    k_smc  = smc(lam_cm)                # A(λ)/A_V for SMC

Integrated with the fitter (preferred interface)::

    from VegasAfterglow import Fitter, ParamDef, Scale

    fitter = Fitter(z=1.58, lumi_dist=3.364e28,
                    jet="tophat", medium="ism",
                    extinction="smc")

    params = [
        ...
        ParamDef("A_V", 0, 3, Scale.linear),
        ...
    ]

References
----------
Pei (1992), ApJ 395, 130.
"""

import numpy as np

# ---------------------------------------------------------------------------
# Physical constants / reference wavelengths
# ---------------------------------------------------------------------------

_CM_PER_UM: float = 1e-4          # conversion factor: 1 μm = 1e-4 cm
_LYMAN_LIMIT_CM: float = 912e-8   # Lyman limit (912 Å) in cm
_V_BAND_UM: float = 0.55          # V-band reference wavelength [μm]

# ---------------------------------------------------------------------------
# Pei (1992) six-component parameter tables
# ---------------------------------------------------------------------------
# Each row: [amplitude, peak wavelength (μm), exponent n, offset b]
# Formula per component:
#   f_i(λ) = amp_i / ( (λ/λ_i)^n_i + (λ_i/λ)^n_i + b_i )
#
# Parameters from Table 1 of Pei (1992, ApJ 395, 130), following the
# component assignments in Gordon et al. (dust_extinction):
#   BKG  – far-UV background continuum (n=2)
#   FUV  – far-UV bump at 0.08 μm       (n=4)
#   NUV  – near-UV / 2175 Å feature    (n=2)
#   SIL1 – silicate feature at 9.7 μm  (n=2)
#   SIL2 – silicate feature at 18  μm  (n=2)
#   FIR  – far-IR component at 25  μm  (n=2)
# ---------------------------------------------------------------------------

_PEI92_TABLES: dict = {
    "smc": np.array([
        # amp     λ_i (μm)   n      b
        [185.0,   0.042,   2.0,  90.0 ],  # BKG
        [ 27.0,   0.080,   4.0,   5.5 ],  # FUV
        [  0.005, 0.220,   2.0,  -1.95],  # NUV
        [  0.010, 9.700,   2.0,   4.0 ],  # SIL1
        [  0.012, 18.00,   2.0,   4.0 ],  # SIL2
        [  0.030, 25.00,   2.0,   2.0 ],  # FIR
    ]),
    "lmc": np.array([
        [175.0,   0.042,   2.0,  90.0 ],  # BKG
        [ 19.0,   0.080,   4.0,   5.5 ],  # FUV
        [  0.023, 0.220,   2.0,  -1.95],  # NUV
        [  0.005, 9.700,   2.0,   4.0 ],  # SIL1
        [  0.006, 18.00,   2.0,   4.0 ],  # SIL2
        [  0.020, 25.00,   2.0,   2.0 ],  # FIR
    ]),
    "mw": np.array([
        [165.0,   0.042,   2.0,  90.0 ],  # BKG
        [ 14.0,   0.080,   4.0,   5.5 ],  # FUV
        [  0.045, 0.220,   2.0,  -1.95],  # NUV
        [  0.002, 9.700,   2.0,   4.0 ],  # SIL1
        [  0.002, 18.00,   2.0,   4.0 ],  # SIL2
        [  0.012, 25.00,   2.0,   2.0 ],  # FIR
    ]),
}


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _pei92_raw(lam_um: np.ndarray, table: np.ndarray) -> np.ndarray:
    """Evaluate the raw (unnormalised) Pei (1992) sum at *lam_um* [μm]."""
    result = np.zeros_like(lam_um, dtype=np.float64)
    for amp, lam_i, n, b in table:
        x = lam_um / lam_i
        result += amp / (x**n + x**(-n) + b)
    return result


# Pre-compute V-band normalization constants (once at import time)
_K_V: dict = {
    law: float(_pei92_raw(np.array([_V_BAND_UM]), table)[0])
    for law, table in _PEI92_TABLES.items()
}


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def pei92(lam_cm, law: str) -> np.ndarray:
    """Pei (1992) dust extinction law: k(λ) = A(λ)/A_V.

    Evaluates the six-component Pei (1992) extinction curve for the chosen
    environment and normalises so that ``k(V) = 1`` exactly (V-band
    reference wavelength 0.55 μm).  Wavelengths below the Lyman limit
    (912 Å = 9.12 × 10⁻⁸ cm) are set to zero so that X-ray data is not
    attenuated.

    Parameters
    ----------
    lam_cm : array_like
        Rest-frame wavelengths in **centimetres**.
    law : {"smc", "lmc", "mw"}
        Dust environment.  ``"smc"`` (steep UV rise, negligible 2175 Å
        bump) is the most common choice for high-redshift GRB hosts.

    Returns
    -------
    k : ndarray
        Dimensionless extinction curve A(λ)/A_V, same shape as *lam_cm*.
        Zero below 912 Å.

    References
    ----------
    Pei (1992), ApJ 395, 130.

    Examples
    --------
    >>> import numpy as np
    >>> from VegasAfterglow.extinction import pei92
    >>> lam = np.array([5.5e-5])   # V band (0.55 μm) in cm
    >>> pei92(lam, "smc")
    array([1.])
    """
    lam_cm = np.asarray(lam_cm, dtype=np.float64)
    law = law.lower()
    if law not in _PEI92_TABLES:
        raise ValueError(
            f"Unknown Pei (1992) law: {law!r}. "
            f"Expected one of {sorted(_PEI92_TABLES)}."
        )
    table = _PEI92_TABLES[law]
    k_v = _K_V[law]

    lam_um = lam_cm / _CM_PER_UM          # cm → μm
    k = _pei92_raw(lam_um, table) / k_v

    # Zero below the Lyman limit (gas-phase X-ray opacity is handled separately)
    k[lam_cm < _LYMAN_LIMIT_CM] = 0.0
    return k


def smc(lam_cm) -> np.ndarray:
    """Pei (1992) SMC dust law: k(λ) = A(λ)/A_V.

    Steep far-UV rise, negligible 2175 Å bump.  The most commonly used
    extinction law for high-redshift GRB host galaxies.

    Parameters
    ----------
    lam_cm : array_like
        Rest-frame wavelengths in centimetres.

    Returns
    -------
    k : ndarray
        A(λ)/A_V, same shape as *lam_cm*.  Zero below 912 Å.

    References
    ----------
    Pei (1992), ApJ 395, 130.
    """
    return pei92(lam_cm, "smc")


def lmc(lam_cm) -> np.ndarray:
    """Pei (1992) LMC dust law: k(λ) = A(λ)/A_V.

    Moderate far-UV rise and a weak 2175 Å bump.

    Parameters
    ----------
    lam_cm : array_like
        Rest-frame wavelengths in centimetres.

    Returns
    -------
    k : ndarray
        A(λ)/A_V, same shape as *lam_cm*.  Zero below 912 Å.

    References
    ----------
    Pei (1992), ApJ 395, 130.
    """
    return pei92(lam_cm, "lmc")


def mw(lam_cm) -> np.ndarray:
    """Pei (1992) Milky Way dust law: k(λ) = A(λ)/A_V.

    Strong 2175 Å feature and moderate far-UV rise.  Suitable for
    MW-like dust in GRB host galaxies.

    Parameters
    ----------
    lam_cm : array_like
        Rest-frame wavelengths in centimetres.

    Returns
    -------
    k : ndarray
        A(λ)/A_V, same shape as *lam_cm*.  Zero below 912 Å.

    References
    ----------
    Pei (1992), ApJ 395, 130.
    """
    return pei92(lam_cm, "mw")


# ---------------------------------------------------------------------------
# Registry used by Fitter._resolve_extinction
# ---------------------------------------------------------------------------

BUILTIN_LAWS: dict = {
    "smc": smc,
    "lmc": lmc,
    "mw":  mw,
}
