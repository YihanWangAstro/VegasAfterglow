"""Corner sweep of the model parameter space.

Every jet variant, medium, radiation switch, observer geometry, and output
method the Python API exposes is exercised once with small grids. Assertions
are structural: outputs are finite, shapes match the inputs, and the total
flux is positive where physics guarantees it. Quantitative regression checks
live in test_golden.py.
"""

import math

import numpy as np
import pytest

import VegasAfterglow as va

pytestmark = pytest.mark.corners

T = np.logspace(2, 7, 20)
NU = np.full_like(T, 1e14)
NUS = np.array([1e9, 1e14, 1e17])


def _observer():
    return va.Observer(lumi_dist=3e28, z=1.0, theta_obs=0.0)


def _observer_off_axis():
    return va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.4)


JETS = {
    "tophat": lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300),
    "tophat_spread": lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300, spreading=True),
    "tophat_thick": lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300, duration=1000),
    "tophat_magnetar": lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300,
                                            magnetar=va.Magnetar(L0=1e47, t0=1e3, q=2)),
    "gaussian": lambda: va.GaussianJet(theta_c=0.05, E_iso=1e53, Gamma0=300),
    "powerlaw": lambda: va.PowerLawJet(theta_c=0.05, E_iso=1e53, Gamma0=300, k_e=2, k_g=1),
    "two_component": lambda: va.TwoComponentJet(theta_c=0.05, E_iso=1e53, Gamma0=300,
                                                theta_w=0.3, E_iso_w=1e51, Gamma0_w=30),
    "step_powerlaw": lambda: va.StepPowerLawJet(theta_c=0.05, E_iso=1e53, Gamma0=300,
                                                E_iso_w=1e51, Gamma0_w=30, k_e=2, k_g=1),
    "powerlaw_wing": lambda: va.PowerLawWing(theta_c=0.05, E_iso_w=1e52, Gamma0_w=100, k_e=2, k_g=1),
    "custom_ejecta": lambda: va.Ejecta(
        E_iso=lambda phi, theta: 1e53 * math.exp(-((theta / 0.1) ** 2)),
        Gamma0=lambda phi, theta: 1 + 299 * math.exp(-((theta / 0.1) ** 2)),
        sigma0=lambda phi, theta: 0.1,
        E_dot=lambda phi, theta, t: 1e47 if t < 100 else 0.0,
        M_dot=lambda phi, theta, t: 1e20 if t < 100 else 0.0,
        duration=100,
    ),
}

MEDIA = {
    "ism": lambda: va.ISM(n_ism=1.0),
    "ism_thin": lambda: va.ISM(n_ism=1e-4),
    "wind": lambda: va.Wind(A_star=0.1),
    "wind_full": lambda: va.Wind(A_star=0.5, n_ism=1.0, n0=1e6, k_m=1.5),
    "custom": lambda: va.Medium(rho=lambda phi, theta, r: 1.67e-24 / (1 + (r / 1e17) ** 2)),
}

RADS = {
    "plain": lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    "ssc": lambda: va.Radiation(eps_e=0.1, eps_B=1e-4, p=2.3, ssc=True),
    "ssc_kn": lambda: va.Radiation(eps_e=0.1, eps_B=1e-4, p=2.3, ssc=True, kn=True),
    "p_near2": lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.02),
    "p_steep": lambda: va.Radiation(eps_e=0.3, eps_B=0.3, p=2.9),
    "xi_e": lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=0.1),
}


def _build_model(jet="tophat", medium="ism", rad="plain", observer=_observer, rvs=None, **kwargs):
    rvs_rad = RADS[rvs]() if rvs is not None else None
    return va.Model(jet=JETS[jet](), medium=MEDIA[medium](), observer=observer(),
                    fwd_rad=RADS[rad](), rvs_rad=rvs_rad, **kwargs)


def _assert_light_curve(flux):
    total = np.asarray(flux.total)
    assert total.shape == T.shape
    assert np.all(np.isfinite(total))
    assert np.all(total > 0)
    return total


@pytest.mark.parametrize("jet_name", sorted(JETS))
def test_jet_corner(jet_name):
    """Each jet variant (tophat, spreading, thick-shell, magnetar, Gaussian, power-law, two-component, step power-law, wing, custom Ejecta) produces a light curve matching the input time shape that is finite and strictly positive."""
    model = _build_model(jet=jet_name)
    _assert_light_curve(model.flux_density(T, NU))


@pytest.mark.parametrize("medium_name", sorted(MEDIA))
def test_medium_corner(medium_name):
    """Each ambient-medium variant (ISM, thin ISM, stellar wind, hybrid wind, custom density profile) produces a finite, strictly positive light curve."""
    model = _build_model(medium=medium_name)
    _assert_light_curve(model.flux_density(T, NU))


@pytest.mark.parametrize("rad_name", sorted(RADS))
def test_radiation_corner(rad_name):
    """Each forward-shock radiation configuration (synchrotron only, SSC, SSC with Klein-Nishina, p near 2, steep p, reduced xi_e) produces a finite, strictly positive light curve."""
    model = _build_model(rad=rad_name)
    _assert_light_curve(model.flux_density(T, NU))


@pytest.mark.parametrize("rad_name", sorted(RADS))
def test_reverse_shock_corner(rad_name):
    """A thick-shell jet with each reverse-shock radiation configuration keeps the total flux finite and positive while the reverse-shock synchrotron component is finite, non-negative, and shaped like the input times."""
    model = _build_model(jet="tophat_thick", rvs=rad_name)
    flux = model.flux_density(T, NU)
    _assert_light_curve(flux)
    rvs_sync = np.asarray(flux.rvs.sync)
    assert rvs_sync.shape == T.shape
    assert np.all(np.isfinite(rvs_sync))
    assert np.all(rvs_sync >= 0)


@pytest.mark.parametrize("jet_name", ["tophat", "gaussian", "two_component"])
def test_off_axis_corner(jet_name):
    """For an observer at theta_obs=0.4 outside each jet's core (tophat, Gaussian, two-component), both the light curve and the (n_nu, n_t) flux grid are finite and strictly positive with the expected shapes."""
    model = _build_model(jet=jet_name, observer=_observer_off_axis)
    _assert_light_curve(model.flux_density(T, NU))
    grid = np.asarray(model.flux_density_grid(T, NUS).total)
    assert grid.shape == (NUS.size, T.size)
    assert np.all(np.isfinite(grid))
    assert np.all(grid > 0)


def test_non_axisymmetric_custom_jet():
    """A custom Ejecta jet solved with axisymmetric=False and an off-axis observer produces a finite, strictly positive light curve."""
    model = _build_model(jet="custom_ejecta", observer=_observer_off_axis, axisymmetric=False)
    _assert_light_curve(model.flux_density(T, NU))


def test_high_resolution_tight_rtol():
    """Raising the grid resolutions above their defaults and loosening the solver rtol to 1e-4 still yields a finite, strictly positive light curve."""
    model = _build_model(resolutions=(0.3, 1.0, 20), rtol=1e-4)
    _assert_light_curve(model.flux_density(T, NU))


def test_method_flux_density():
    """flux_density on the baseline model returns a total-flux array matching the input time shape that is finite and strictly positive."""
    model = _build_model()
    _assert_light_curve(model.flux_density(T, NU))


def test_method_flux_density_grid():
    """flux_density_grid returns a (n_nu, n_t) total-flux grid that is finite and strictly positive."""
    model = _build_model()
    grid = np.asarray(model.flux_density_grid(T, NUS).total)
    assert grid.shape == (NUS.size, T.size)
    assert np.all(np.isfinite(grid))
    assert np.all(grid > 0)


def test_method_flux_band():
    """Band-integrated flux over 1e17-1e19 Hz matches the input time shape and is finite and strictly positive."""
    model = _build_model()
    band = np.asarray(model.flux(T, 1e17, 1e19, 8).total)
    assert band.shape == T.shape
    assert np.all(np.isfinite(band))
    assert np.all(band > 0)


def test_method_details():
    """details returns forward-shock dynamics with finite Lorentz factor Gamma >= 1, strictly positive radius, and finite observer times."""
    model = _build_model()
    details = model.details(t_min=1e2, t_max=1e6)
    gamma = np.asarray(details.fwd.Gamma)
    radius = np.asarray(details.fwd.r)
    t_obs = np.asarray(details.fwd.t_obs)
    assert np.all(np.isfinite(gamma))
    assert np.all(gamma >= 1)
    assert np.all(radius > 0)
    assert np.all(np.isfinite(t_obs))


def test_method_flux_density_exposures():
    """Exposure-averaged flux over finite 10 s exposures returns one finite, strictly positive value per epoch."""
    model = _build_model()
    n_expo = 5
    flux = model.flux_density_exposures(T[:n_expo], NU[:n_expo], np.full(n_expo, 10.0), 4)
    total = np.asarray(flux.total)
    assert total.shape == (n_expo,)
    assert np.all(np.isfinite(total))
    assert np.all(total > 0)


def test_jet_property_accessors():
    """jet_E_iso and jet_Gamma0 on a Gaussian jet return arrays matching the theta grid with E_iso > 0 and Gamma0 > 1 everywhere."""
    model = _build_model(jet="gaussian")
    theta = np.array([0.0, 0.02, 0.05])
    e_iso = model.jet_E_iso(0.0, theta)
    gamma0 = model.jet_Gamma0(0.0, theta)
    assert np.asarray(e_iso).shape == theta.shape
    assert np.asarray(gamma0).shape == theta.shape
    assert np.all(np.asarray(e_iso) > 0)
    assert np.all(np.asarray(gamma0) > 1)


def test_sky_image():
    """sky_image for an off-axis observer returns a finite intensity map with at least one lit pixel, the requested pixel dimension, a positive pixel solid angle, and at least four extent entries."""
    import numpy as np

    import VegasAfterglow as va

    m = va.Model(
        jet=va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300),
        medium=va.ISM(n_ism=1.0),
        observer=va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.2),
        fwd_rad=va.Radiation(eps_e=0.1, eps_B=1e-3, p=2.3),
    )
    img = m.sky_image(np.array([1e5]), 1e15, fov=3e18, npixel=32)
    arr = np.asarray(img.image)
    assert np.all(np.isfinite(arr))
    assert np.any(arr > 0)
    assert arr.shape[-1] == 32
    assert img.pixel_solid_angle > 0
    assert len(np.asarray(img.extent)) >= 4


# ---------------------------------------------------------------------------
# Extreme parameter values: validity-boundary corners. Every entry was probed
# to produce finite, positive flux over t = 1 s .. 3e2 yr before being added.
# ---------------------------------------------------------------------------

T_WIDE = np.logspace(0, 10, 30)  # coasting through deep-Newtonian
NU_MID = np.full_like(T_WIDE, 1e15)

EXTREMES = {
    "Gamma0_transrel": dict(jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=2)),
    "Gamma0_5000": dict(jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=5000)),
    "theta_c_needle": dict(jet=lambda: va.TophatJet(theta_c=0.01, E_iso=1e53, Gamma0=300)),
    "theta_c_spherical": dict(jet=lambda: va.TophatJet(theta_c=np.pi / 2, E_iso=1e53, Gamma0=300)),
    "E_iso_1e48": dict(jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e48, Gamma0=300)),
    "E_iso_1e55": dict(jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e55, Gamma0=300)),
    "obs_on_jet_edge": dict(observer=lambda: va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.1)),
    "obs_equatorial": dict(observer=lambda: va.Observer(lumi_dist=3e28, z=0.5, theta_obs=1.5)),
    "z_zero": dict(observer=lambda: va.Observer(lumi_dist=3e28, z=0.0, theta_obs=0.0)),
    "z_ten": dict(observer=lambda: va.Observer(lumi_dist=3e29, z=10.0, theta_obs=0.0)),
    "n_ism_igm": dict(medium=lambda: va.ISM(n_ism=1e-6)),
    "n_ism_dense": dict(medium=lambda: va.ISM(n_ism=1e6)),
    "A_star_heavy": dict(medium=lambda: va.Wind(A_star=10)),
    "A_star_light": dict(medium=lambda: va.Wind(A_star=1e-3)),
    "eps_near_one": dict(fwd_rad=lambda: va.Radiation(eps_e=0.99, eps_B=0.99, p=2.3)),
    "eps_tiny": dict(fwd_rad=lambda: va.Radiation(eps_e=1e-6, eps_B=1e-6, p=2.3)),
    "p_2.01": dict(fwd_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.01)),
    "p_3.5": dict(fwd_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=3.5)),
    "p_hard_1.5": dict(fwd_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=1.5)),
    "xi_e_1e-3": dict(fwd_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=1e-3)),
    "ultralong_rvs": dict(
        jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300, duration=1e5),
        rvs_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    ),
    "sigma_fwd_only": dict(
        jet=lambda: va.Ejecta(E_iso=lambda p_, th: 1e53, Gamma0=lambda p_, th: 300.0,
                              sigma0=lambda p_, th: 10.0),
    ),
}


def _build_extreme(overrides):
    parts = dict(
        jet=lambda: va.TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300),
        medium=lambda: va.ISM(n_ism=1.0),
        observer=lambda: va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=lambda: va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    )
    parts.update({k: v for k, v in overrides.items() if k != "rvs_rad"})
    kw = {k: f() for k, f in parts.items()}
    if "rvs_rad" in overrides:
        kw["rvs_rad"] = overrides["rvs_rad"]()
    return va.Model(**kw)


@pytest.mark.parametrize("name", sorted(EXTREMES))
def test_extreme_corner(name):
    """Each validity-boundary extreme (trans-relativistic to Gamma0=5000, needle to spherical jets, IGM to dense media, hard to steep p, magnetized shell, ultralong reverse shock) yields finite, strictly positive flux from 1 s through the deep-Newtonian phase."""
    model = _build_extreme(EXTREMES[name])
    f = np.asarray(model.flux_density(T_WIDE, NU_MID).total)
    assert np.all(np.isfinite(f)), name
    assert np.all(f > 0), name


@pytest.mark.parametrize("sigma", [0.1, 10.0])
def test_magnetized_ejecta_with_reverse_shock_decays(sigma):
    """Regression for the sigma>0 + rvs_rad runaway (fixed 2026-07-02): at early
    times Gamma ~= Gamma4 makes the magnetized jump ratio -> 1, and the
    reverse-shock crossing rate divided by (Gamma*comp_ratio/Gamma4 - 1) -> 0/0;
    the NaN then silently froze dGamma/dt at zero (eternal coasting). Fixed by
    the shock-penetration gate in FRShockEqn::compute_dx3_dt."""
    m = va.Model(
        jet=va.Ejecta(E_iso=lambda p_, th: 1e53, Gamma0=lambda p_, th: 300.0,
                      sigma0=lambda p_, th, s=sigma: s),
        medium=va.ISM(n_ism=1.0),
        observer=va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
        rvs_rad=va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    )
    f = np.asarray(m.flux_density(T_WIDE, NU_MID).total)
    assert np.all(np.isfinite(f))
    # a physical afterglow must decay after its peak
    assert f[-1] < np.max(f) * 1e-2


@pytest.mark.parametrize("sigma", [0.3, 1.0, 2.0])
@pytest.mark.parametrize("npts", [15, 30, 60])
def test_magnetized_knife_edge_sigmas(sigma, npts):
    """Regression for the penetration->0+ knife edge (fixed 2026-07-02): for
    magnetized shells the Zhang-Kobayashi crossing rate diverges as the
    penetration factor (Gamma*comp_ratio/Gamma4 - 1) -> 0+, and whether a run
    hit the singularity depended on the time grid's t0 (sigma=0.999 passed
    while sigma=1.0 failed). Fixed by capping the comoving consumption rate at
    the fast magnetosonic speed of the magnetized upstream in
    FRShockEqn::compute_dx3_dt (unreachable for sigma=0, where the penetration
    factor is algebraically >= 1). Parametrized over time grids because the
    original failure was grid-t0 sensitive."""
    t = np.logspace(0, 10, npts)
    nu = np.full_like(t, 1e15)
    m = va.Model(
        jet=va.Ejecta(E_iso=lambda p_, th: 1e53, Gamma0=lambda p_, th: 300.0,
                      sigma0=lambda p_, th, s=sigma: s),
        medium=va.ISM(n_ism=1.0),
        observer=va.Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
        rvs_rad=va.Radiation(eps_e=0.1, eps_B=0.01, p=2.3),
    )
    f = np.asarray(m.flux_density(t, nu).total)
    assert np.all(np.isfinite(f))
    assert np.max(f) > 0
    assert f[-1] < np.max(f) * 1e-2  # decays after peak
    assert t[int(np.argmax(f))] < 1e5  # peaks early, no late runaway
