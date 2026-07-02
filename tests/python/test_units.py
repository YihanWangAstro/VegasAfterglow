"""Magnitude/flux conversions and named bands in VegasAfterglow.units."""
import numpy as np
import pytest

from VegasAfterglow import units


# ----------------- AB magnitudes -----------------

def test_ABmag_zero_is_AB_zeropoint():
    # AB mag 0 corresponds to 3631 Jy = 3.631e-20 erg/cm^2/s/Hz
    assert units.ABmag_to_cgs(0.0) == pytest.approx(3.631e-20, rel=1e-6)


def test_ABmag_roundtrip_scalar():
    for mag in (-5.0, 0.0, 17.3, 30.0):
        assert units.cgs_to_ABmag(units.ABmag_to_cgs(mag)) == pytest.approx(mag, abs=1e-12)


def test_ABmag_roundtrip_array():
    mags = np.linspace(10, 25, 7)
    out = units.cgs_to_ABmag(units.ABmag_to_cgs(mags))
    np.testing.assert_allclose(out, mags, atol=1e-12)


def test_ABmag_five_mags_is_factor_100():
    assert units.ABmag_to_cgs(10.0) / units.ABmag_to_cgs(15.0) == pytest.approx(100.0, rel=1e-12)


# ----------------- Vega / ST magnitudes -----------------

def test_Vegamag_roundtrip():
    filt = next(iter(units._VEGA_FILTERS))
    for mag in (0.0, 12.5, 20.0):
        assert units.cgs_to_Vegamag(units.Vegamag_to_cgs(mag, filt), filt) == pytest.approx(mag, abs=1e-10)


def test_STmag_roundtrip():
    filt = next(iter(units._ST_FILTERS))
    for mag in (0.0, 12.5, 20.0):
        assert units.cgs_to_STmag(units.STmag_to_cgs(mag, filt), filt) == pytest.approx(mag, abs=1e-10)


def test_unknown_filter_raises():
    with pytest.raises((KeyError, ValueError)):
        units.Vegamag_to_cgs(15.0, "definitely-not-a-filter")


# ----------------- named bands / filters -----------------

def test_band_returns_positive_range():
    nu_min, nu_max = units.band("XRT")
    assert 0 < nu_min < nu_max


def test_band_unknown_name_raises():
    with pytest.raises(ValueError, match="Unknown band"):
        units.band("definitely-not-a-band")


def test_filter_returns_positive_frequency():
    filt = next(iter(units._VEGA_FILTERS))
    nu = units.filter(filt)
    assert np.all(np.asarray(nu) > 0)
