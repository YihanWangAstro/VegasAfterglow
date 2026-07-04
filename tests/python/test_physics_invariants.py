"""Exact physical invariants of the flux pipeline.

These hold to numerical precision by construction (verified empirically to
<= 1e-9 relative); any violation indicates a real bug, not physics tolerance.
"""
import numpy as np
import pytest

import VegasAfterglow as va
from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet

# Exact analytic invariants verify to float precision on the exact-libm build;
# under AFTERGLOW_FAST_MATH the two compared paths evaluate the polynomial
# kernels (rel err ~1e-7) with different arguments, so the achievable agreement
# is bounded by the kernel accuracy, not by the invariant.
EXACT_RTOL = 1e-6 if getattr(va.VegasAfterglowC, "fast_math_enabled", False) else 1e-9

P = 2.5
T = np.logspace(3.5, 5.5, 16)
NU = np.full_like(T, 1e15)


def _model(lumi_dist=3e28, z=0.5, axisymmetric=True, ssc=False, rvs=False):
    kw = {}
    if rvs:
        kw["rvs_rad"] = Radiation(eps_e=0.1, eps_B=0.01, p=P)
    return Model(
        jet=TophatJet(theta_c=0.3, E_iso=1e53, Gamma0=300),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=lumi_dist, z=z, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=P, ssc=ssc),
        axisymmetric=axisymmetric,
        **kw,
    )


def test_flux_scales_exactly_with_inverse_distance_squared():
    """Doubling the luminosity distance reduces the total flux density by exactly a factor of 4 (F proportional to 1/D_L^2) to 1e-9 relative."""
    f1 = np.asarray(_model(lumi_dist=3e28).flux_density(T, NU).total)
    f2 = np.asarray(_model(lumi_dist=6e28).flux_density(T, NU).total)
    np.testing.assert_allclose(f1 / f2, 4.0, rtol=EXACT_RTOL)


def test_total_equals_sum_of_components():
    """With SSC and reverse shock enabled, the total flux equals the sum of forward and reverse synchrotron plus SSC components to 1e-12 relative."""
    f = _model(ssc=True, rvs=True).flux_density(T, NU)
    total = np.asarray(f.total)
    parts = (np.asarray(f.fwd.sync) + np.asarray(f.fwd.ssc)
             + np.asarray(f.rvs.sync) + np.asarray(f.rvs.ssc))
    np.testing.assert_allclose(total, parts, rtol=1e-12)


def test_axisymmetric_flag_consistent_on_axis():
    """For an on-axis observer (theta_obs=0) the axisymmetric fast path and the full 3D integration yield the same flux to 1e-9 relative."""
    fa = np.asarray(_model(axisymmetric=True).flux_density(T, NU).total)
    fb = np.asarray(_model(axisymmetric=False).flux_density(T, NU).total)
    np.testing.assert_allclose(fa, fb, rtol=EXACT_RTOL)


def test_redshift_transformation_invariance():
    """F(nu, t; z2) == F(nu (1+z2)/(1+z1), t (1+z1)/(1+z2); z1) * (1+z2)/(1+z1)
    at fixed luminosity distance -- the exact cosmological transformation."""
    z1, z2 = 0.2, 1.4
    m1, m2 = _model(z=z1), _model(z=z2)
    scale = (1 + z2) / (1 + z1)
    f_transformed = np.asarray(m1.flux_density(T / scale, NU * scale).total) * scale
    f_direct = np.asarray(m2.flux_density(T, NU).total)
    np.testing.assert_allclose(f_transformed, f_direct, rtol=EXACT_RTOL)


def test_disabled_components_are_zero():
    """With SSC and the reverse shock disabled, the fwd.ssc, rvs.sync, and rvs.ssc components are exactly zero while fwd.sync remains strictly positive."""
    f = _model(ssc=False, rvs=False).flux_density(T, NU)
    assert np.all(np.asarray(f.fwd.ssc) == 0)
    assert np.all(np.asarray(f.rvs.sync) == 0)
    assert np.all(np.asarray(f.rvs.ssc) == 0)
    assert np.all(np.asarray(f.fwd.sync) > 0)


def test_series_and_grid_evaluations_agree():
    """flux_density evaluated along a (t, nu) series matches the corresponding row of flux_density_grid at the same frequency to 1e-12 relative."""
    m = _model()
    a = np.asarray(m.flux_density(T, NU).total)
    b = np.asarray(m.flux_density_grid(T, np.array([1e15])).total)[0]
    np.testing.assert_allclose(a, b, rtol=1e-12)


def test_band_flux_matches_integrated_flux_density():
    """Band-integrated flux from flux() agrees with the trapezoidal frequency integral of flux_density_grid over 1e14-1e15 Hz within 1% (measured 1.7e-4), pinning the quadrature contract."""
    m = _model()
    tb = np.logspace(3, 5, 8)
    band = np.asarray(m.flux(tb, 1e14, 1e15, 16).total)
    nu_fine = np.logspace(14, 15, 60)
    grid = np.asarray(m.flux_density_grid(tb, nu_fine).total)
    integrated = np.trapezoid(grid, nu_fine, axis=0)
    # measured agreement 1.7e-4; 1% guards the quadrature contract
    np.testing.assert_allclose(band, integrated, rtol=1e-2)


def test_flux_positive_and_finite_everywhere():
    """With SSC and reverse shock enabled, the total flux density is finite and strictly positive at every sampled time."""
    f = _model(ssc=True, rvs=True).flux_density(T, NU)
    total = np.asarray(f.total)
    assert np.all(np.isfinite(total))
    assert np.all(total > 0)
