"""Pei (1992) extinction laws in VegasAfterglow.extinction."""
import numpy as np
import pytest

from VegasAfterglow import extinction

V_BAND_CM = 5477e-8  # matches extinction._V_BAND_CM normalization wavelength
LYMAN_LIMIT_CM = 912e-8


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_k_at_V_is_unity(law):
    assert extinction.pei92(extinction._V_BAND_CM, law) == pytest.approx(1.0, rel=1e-12)


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_zero_below_lyman_limit(law):
    lam = np.array([0.5, 0.9, 0.99]) * LYMAN_LIMIT_CM
    np.testing.assert_array_equal(extinction.pei92(lam, law), 0.0)


@pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
def test_uv_extinction_exceeds_optical(law):
    # Dust extinction rises toward the UV for all Pei92 profiles
    k_uv = extinction.pei92(2000e-8, law)
    k_opt = extinction.pei92(7000e-8, law)
    assert k_uv > k_opt > 0


def test_vectorized_shape_preserved():
    lam = np.geomspace(1000e-8, 3e-4, 25).reshape(5, 5)
    out = extinction.pei92(lam, "smc")
    assert out.shape == lam.shape
    assert np.all(np.isfinite(out))
    assert np.all(out >= 0)


def test_wrapper_functions_match_pei92():
    lam = np.geomspace(1500e-8, 1e-4, 10)
    np.testing.assert_array_equal(extinction.smc(lam), extinction.pei92(lam, "smc"))
    np.testing.assert_array_equal(extinction.lmc(lam), extinction.pei92(lam, "lmc"))
    np.testing.assert_array_equal(extinction.mw(lam), extinction.pei92(lam, "mw"))


def test_builtin_laws_registry():
    assert set(extinction.BUILTIN_LAWS) == {"smc", "lmc", "mw"}
    for fn in extinction.BUILTIN_LAWS.values():
        assert fn(V_BAND_CM) > 0
