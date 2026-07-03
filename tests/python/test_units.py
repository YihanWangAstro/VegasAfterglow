"""Magnitude/flux conversions and named bands in VegasAfterglow.units."""
import numpy as np
import pytest

from VegasAfterglow import units


# ----------------- AB magnitudes -----------------

def test_ABmag_zero_is_AB_zeropoint():
    """AB magnitude 0 converts to the AB zero-point flux density of 3631 Jy (3.631e-20 erg/cm^2/s/Hz)."""
    # AB mag 0 corresponds to 3631 Jy = 3.631e-20 erg/cm^2/s/Hz
    assert units.ABmag_to_cgs(0.0) == pytest.approx(3.631e-20, rel=1e-6)


def test_ABmag_roundtrip_scalar():
    """cgs_to_ABmag inverts ABmag_to_cgs to within 1e-12 mag for scalar magnitudes from -5 to 30."""
    for mag in (-5.0, 0.0, 17.3, 30.0):
        assert units.cgs_to_ABmag(units.ABmag_to_cgs(mag)) == pytest.approx(mag, abs=1e-12)


def test_ABmag_roundtrip_array():
    """AB magnitude to flux density round-trip recovers a numpy array of magnitudes elementwise to within relative tolerance 1e-7 (atol 1e-12 mag)."""
    mags = np.linspace(10, 25, 7)
    out = units.cgs_to_ABmag(units.ABmag_to_cgs(mags))
    np.testing.assert_allclose(out, mags, atol=1e-12)


def test_ABmag_five_mags_is_factor_100():
    """A 5-magnitude difference in AB magnitude corresponds to exactly a factor of 100 in flux density."""
    assert units.ABmag_to_cgs(10.0) / units.ABmag_to_cgs(15.0) == pytest.approx(100.0, rel=1e-12)


# ----------------- Vega / ST magnitudes -----------------

def test_Vegamag_roundtrip():
    """cgs_to_Vegamag inverts Vegamag_to_cgs to within 1e-10 mag for a registered Vega-system filter."""
    filt = next(iter(units._VEGA_FILTERS))
    for mag in (0.0, 12.5, 20.0):
        assert units.cgs_to_Vegamag(units.Vegamag_to_cgs(mag, filt), filt) == pytest.approx(mag, abs=1e-10)


def test_STmag_roundtrip():
    """cgs_to_STmag inverts STmag_to_cgs to within 1e-10 mag for a registered ST-system filter."""
    filt = next(iter(units._ST_FILTERS))
    for mag in (0.0, 12.5, 20.0):
        assert units.cgs_to_STmag(units.STmag_to_cgs(mag, filt), filt) == pytest.approx(mag, abs=1e-10)


def test_unknown_filter_raises():
    """Vegamag_to_cgs raises KeyError or ValueError when given an unrecognized filter name."""
    with pytest.raises((KeyError, ValueError)):
        units.Vegamag_to_cgs(15.0, "definitely-not-a-filter")


# ----------------- named bands / filters -----------------

def test_band_returns_positive_range():
    """band('XRT') returns a strictly ordered positive frequency range with 0 < nu_min < nu_max."""
    nu_min, nu_max = units.band("XRT")
    assert 0 < nu_min < nu_max


def test_band_unknown_name_raises():
    """band raises ValueError with an 'Unknown band' message for an unrecognized band name."""
    with pytest.raises(ValueError, match="Unknown band"):
        units.band("definitely-not-a-band")


def test_filter_returns_positive_frequency():
    """filter returns a strictly positive effective frequency in Hz for a registered Vega filter."""
    filt = next(iter(units._VEGA_FILTERS))
    nu = units.filter(filt)
    assert np.all(np.asarray(nu) > 0)
