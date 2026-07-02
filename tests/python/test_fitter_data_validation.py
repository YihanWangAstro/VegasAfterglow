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
    f = _fitter()
    f.add_flux_density(
        nu=5e14,
        t=np.array([1e3, 1e4, 1e5]),
        f_nu=np.array([1e-26, 1e-26, 1e-26]),
        err=np.array([1e-28, 1e-28, 1e-28]),
    )


def test_add_flux_density_with_weights():
    f = _fitter()
    f.add_flux_density(
        nu=5e14,
        t=np.array([1e3, 1e4]),
        f_nu=np.array([1e-26, 1e-26]),
        err=np.array([1e-28, 1e-28]),
        weights=np.array([1.0, 2.0]),
    )


def test_add_spectrum_happy():
    f = _fitter()
    f.add_spectrum(
        t=1e4,
        nu=np.array([1e14, 5e14, 1e15]),
        f_nu=np.array([2e-26, 1e-26, 5e-27]),
        err=np.array([2e-28, 1e-28, 5e-29]),
    )


def test_add_flux_happy():
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
    f = _fitter()
    with pytest.raises(ValueError, match="nu"):
        f.add_flux_density(nu=bad_nu, t=np.array([1e3]), f_nu=np.array([1e-26]), err=np.array([1e-28]))


def test_add_flux_density_shape_mismatch():
    f = _fitter()
    with pytest.raises(ValueError, match="same shape"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4, 1e5]),
            f_nu=np.array([1e-26, 1e-26]),  # one short
            err=np.array([1e-28, 1e-28, 1e-28]),
        )


def test_add_flux_density_empty():
    f = _fitter()
    with pytest.raises(ValueError, match="empty"):
        f.add_flux_density(nu=5e14, t=np.array([]), f_nu=np.array([]), err=np.array([]))


def test_add_flux_density_nan_flux():
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
    f = _fitter()
    with pytest.raises(ValueError, match="err"):
        f.add_flux_density(
            nu=5e14,
            t=np.array([1e3, 1e4]),
            f_nu=np.array([1e-26, 1e-26]),
            err=bad_err,
        )


def test_add_flux_density_bad_weights_shape():
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
    f = _fitter()
    with pytest.raises(ValueError, match="^add_spectrum: t"):
        f.add_spectrum(
            t=bad_t,
            nu=np.array([1e14, 5e14]),
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
        )


def test_add_spectrum_bad_nu():
    f = _fitter()
    with pytest.raises(ValueError, match="nu"):
        f.add_spectrum(
            t=1e4,
            nu=np.array([1e14, -5e14]),  # negative
            f_nu=np.array([1e-26, 1e-26]),
            err=np.array([1e-28, 1e-28]),
        )


def test_add_spectrum_shape_mismatch():
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
    f = _fitter()
    with pytest.raises(ValueError, match="band"):
        f.add_flux(band=band, t=np.array([1e3]), flux=np.array([1e-12]), err=np.array([1e-14]))


def test_add_flux_bad_num_points():
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
    f = _fitter()
    with pytest.raises(ValueError, match="err"):
        f.add_flux(
            band=(1e14, 1e15),
            t=np.array([1e3, 1e4]),
            flux=np.array([1e-12, 1e-12]),
            err=np.array([1e-14, -1e-14]),
        )


def test_add_flux_shape_mismatch():
    f = _fitter()
    with pytest.raises(ValueError, match="same shape"):
        f.add_flux(
            band=(1e14, 1e15),
            t=np.array([1e3, 1e4, 1e5]),
            flux=np.array([1e-12, 1e-12]),
            err=np.array([1e-14, 1e-14, 1e-14]),
        )
