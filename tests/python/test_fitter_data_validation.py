"""Data-add boundary validation: Fitter.add_flux_density / add_spectrum / add_flux."""
import numpy as np
import pytest

from VegasAfterglow import Fitter
from VegasAfterglow.units import band as named_band


NAN = float("nan")
INF = float("inf")


def _fitter():
    return Fitter(z=0.5, lumi_dist=1e28)


# ----------------- happy paths -----------------

def test_add_flux_density_happy():
    """add_flux_density accepts a positive scalar frequency with matching finite t, f_nu, and err arrays without raising."""
    f = _fitter()
    f.add_flux_density(
        nu=5e14,
        t=np.array([1e3, 1e4, 1e5]),
        f_nu=np.array([1e-26, 1e-26, 1e-26]),
        err=np.array([1e-28, 1e-28, 1e-28]),
    )


def test_add_flux_density_with_weights():
    """add_flux_density accepts an optional positive weights array of the same shape as the data without raising."""
    f = _fitter()
    f.add_flux_density(
        nu=5e14,
        t=np.array([1e3, 1e4]),
        f_nu=np.array([1e-26, 1e-26]),
        err=np.array([1e-28, 1e-28]),
        weights=np.array([1.0, 2.0]),
    )


def test_add_spectrum_happy():
    """add_spectrum accepts a positive scalar time with matching finite nu, f_nu, and err arrays without raising."""
    f = _fitter()
    f.add_spectrum(
        t=1e4,
        nu=np.array([1e14, 5e14, 1e15]),
        f_nu=np.array([2e-26, 1e-26, 5e-27]),
        err=np.array([2e-28, 1e-28, 5e-29]),
    )


def test_add_flux_happy():
    """add_flux accepts a named instrument band (XRT) resolved to a frequency range with matching t, flux, and err arrays without raising."""
    f = _fitter()
    f.add_flux(
        band=named_band("XRT"),
        t=np.array([1e3, 1e4, 1e5]),
        flux=np.array([1e-12, 1e-12, 1e-12]),
        err=np.array([1e-14, 1e-14, 1e-14]),
    )


# ----------------- add_flux_density failure modes -----------------

@pytest.mark.parametrize("bad_nu", [0, -1e14, NAN, INF])
def test_add_flux_density_bad_nu(bad_nu):
    """add_flux_density raises ValueError naming nu for each non-positive or non-finite frequency (zero, negative, NaN, Inf)."""
    f = _fitter()
    with pytest.raises(ValueError, match="nu"):
        f.add_flux_density(nu=bad_nu, t=np.array([1e3]), f_nu=np.array([1e-26]), err=np.array([1e-28]))


def test_add_flux_density_shape_mismatch():
    """add_flux_density raises ValueError mentioning 'same shape' when f_nu is shorter than the t and err arrays."""
    f = _fitter()
    with pytest.raises(ValueError, match="same shape"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4, 1e5]),
            f_nu=np.array([1e-26, 1e-26]),  # one short
            err=np.array([1e-28, 1e-28, 1e-28]),
        )


def test_add_flux_density_empty():
    """add_flux_density raises ValueError mentioning 'empty' when the t, f_nu, and err arrays contain no points."""
    f = _fitter()
    with pytest.raises(ValueError, match="empty"):
        f.add_flux_density(nu=5e14, t=np.array([]), f_nu=np.array([]), err=np.array([]))


def test_add_flux_density_nan_flux():
    """add_flux_density raises ValueError mentioning 'non-finite' when the f_nu array contains a NaN entry."""
    f = _fitter()
    with pytest.raises(ValueError, match="non-finite"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4]),
            f_nu=np.array([1e-26, NAN]),
            err=np.array([1e-28, 1e-28]),
        )


@pytest.mark.parametrize("bad_err", [
    np.array([1e-28, 0.0]),
    np.array([1e-28, -1e-28]),
    np.array([1e-28, NAN]),
    np.array([1e-28, INF]),
])
def test_add_flux_density_bad_err(bad_err):
    """add_flux_density raises ValueError naming err for each measurement uncertainty that is zero, negative, NaN, or Inf."""
    f = _fitter()
    with pytest.raises(ValueError, match="err"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4]),
            f_nu=np.array([1e-26, 1e-26]),
            err=bad_err,
        )


def test_add_flux_density_bad_weights_shape():
    """add_flux_density raises ValueError naming weights when the weights array is longer than the data arrays."""
    f = _fitter()
    with pytest.raises(ValueError, match="weights"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4]),
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
            weights=np.array([1.0, 2.0, 3.0]),  # one too many
        )


def test_add_flux_density_negative_weights():
    """add_flux_density raises ValueError naming weights when the weights array contains a negative entry."""
    f = _fitter()
    with pytest.raises(ValueError, match="weights"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4]),
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
            weights=np.array([1.0, -0.5]),
        )


# ----------------- add_spectrum failure modes -----------------

@pytest.mark.parametrize("bad_t", [0, -1, NAN, INF])
def test_add_spectrum_bad_t(bad_t):
    """add_spectrum raises a ValueError whose message starts with 'add_spectrum: t' for each non-positive or non-finite observation time (zero, negative, NaN, Inf)."""
    f = _fitter()
    with pytest.raises(ValueError, match="^add_spectrum: t"):
        f.add_spectrum(
            t=bad_t,
            nu=np.array([1e14, 5e14]),
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
        )


def test_add_spectrum_bad_nu():
    """add_spectrum raises ValueError naming nu when the frequency array contains a negative value."""
    f = _fitter()
    with pytest.raises(ValueError, match="nu"):
        f.add_spectrum(
            t=1e4,
            nu=np.array([1e14, -5e14]),  # negative
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
        )


def test_add_spectrum_shape_mismatch():
    """add_spectrum raises ValueError mentioning 'same shape' when f_nu is shorter than the nu and err arrays."""
    f = _fitter()
    with pytest.raises(ValueError, match="same shape"):
        f.add_spectrum(
            t=1e4,
            nu=np.array([1e14, 5e14, 1e15]),
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28, 1e-28]),
        )


# ----------------- add_flux failure modes -----------------

@pytest.mark.parametrize("bad_band", [None, 5e14, (5e14,), (1, 2, 3)])
def test_add_flux_bad_band_format(bad_band):
    """add_flux raises ValueError naming band for each argument that is not a two-element (nu_min, nu_max) pair (None, bare scalar, 1-tuple, 3-tuple)."""
    f = _fitter()
    with pytest.raises(ValueError, match="band"):
        f.add_flux(band=bad_band, t=np.array([1e3]), flux=np.array([1e-12]), err=np.array([1e-14]))


@pytest.mark.parametrize("band", [
    (5e14, 1e14),     # reversed
    (-1e14, 1e14),    # negative
    (0, 1e14),        # zero
    (1e14, 1e14),     # equal
    (NAN, 1e15),
    (1e14, INF),
])
def test_add_flux_bad_band_values(band):
    """add_flux raises ValueError naming band for each frequency pair that is not a strictly increasing positive finite range (reversed, negative, zero, equal edges, NaN, Inf)."""
    f = _fitter()
    with pytest.raises(ValueError, match="band"):
        f.add_flux(band=band, t=np.array([1e3]), flux=np.array([1e-12]), err=np.array([1e-14]))


def test_add_flux_bad_num_points():
    """add_flux raises ValueError naming num_points when only one in-band integration frequency is requested."""
    f = _fitter()
    with pytest.raises(ValueError, match="num_points"):
        f.add_flux(
            band=(1e14, 1e15),
            t=np.array([1e3]),
            flux=np.array([1e-12]),
            err=np.array([1e-14]),
            num_points=1,
        )


def test_add_flux_negative_err():
    """add_flux raises ValueError naming err when the error array contains a negative entry."""
    f = _fitter()
    with pytest.raises(ValueError, match="err"):
        f.add_flux(
            band=(1e14, 1e15),
            t=np.array([1e3, 1e4]),
            flux=np.array([1e-12, 1e-12]),
            err=np.array([1e-14, -1e-14]),
        )


def test_add_flux_shape_mismatch():
    """add_flux raises ValueError mentioning 'same shape' when flux is shorter than the t and err arrays."""
    f = _fitter()
    with pytest.raises(ValueError, match="same shape"):
        f.add_flux(
            band=(1e14, 1e15),
            t=np.array([1e3, 1e4, 1e5]),
            flux=np.array([1e-12, 1e-12]),
            err=np.array([1e-14, 1e-14, 1e-14]),
        )
