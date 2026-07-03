"""Pei (1992) extinction laws in VegasAfterglow.extinction."""
import numpy as np
import pytest

from VegasAfterglow import extinction

V_BAND_CM = 5477e-8  # matches extinction._V_BAND_CM normalization wavelength
LYMAN_LIMIT_CM = 912e-8


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_k_at_V_is_unity(law):
    """Each Pei92 law (SMC/LMC/MW) normalizes to k(λ)=1 at the V-band wavelength to 1e-12 relative precision."""
    assert extinction.pei92(extinction._V_BAND_CM, law) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_zero_below_lyman_limit(law):
    """Each Pei92 law returns exactly zero extinction for wavelengths below the Lyman limit (912 Å)."""
    lam = np.array([0.5, 0.9, 0.99]) * LYMAN_LIMIT_CM
    np.testing.assert_array_equal(extinction.pei92(lam, law), 0.0)


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_uv_extinction_exceeds_optical(law):
    """Each Pei92 law rises toward the UV: k(2000 Å) exceeds k(7000 Å), and both are positive."""
    # Dust extinction rises toward the UV for all Pei92 profiles
    k_uv = extinction.pei92(2000e-8, law)
    k_opt = extinction.pei92(7000e-8, law)
    assert k_uv > k_opt > 0


def test_vectorized_shape_preserved():
    """pei92 on a 2-D wavelength array preserves the input shape and returns only finite, non-negative values."""
    lam = np.geomspace(1000e-8, 3e-4, 25).reshape(5, 5)
    out = extinction.pei92(lam, "smc")
    assert out.shape == lam.shape
    assert np.all(np.isfinite(out))
    assert np.all(out >= 0)


def test_wrapper_functions_match_pei92():
    """The smc/lmc/mw convenience wrappers return values identical to pei92 called with the matching law name."""
    lam = np.geomspace(1500e-8, 1e-4, 10)
    np.testing.assert_array_equal(extinction.smc(lam), extinction.pei92(lam, "smc"))
    np.testing.assert_array_equal(extinction.lmc(lam), extinction.pei92(lam, "lmc"))
    np.testing.assert_array_equal(extinction.mw(lam), extinction.pei92(lam, "mw"))


def test_builtin_laws_registry():
    """BUILTIN_LAWS registry contains exactly smc/lmc/mw, each mapping to a callable positive at the V band."""
    assert set(extinction.BUILTIN_LAWS) == {"smc", "lmc", "mw"}
    for fn in extinction.BUILTIN_LAWS.values():
        assert fn(V_BAND_CM) > 0
