"""Pybind constructor input validation: NaN / inf / non-positive must raise early."""
import math

import numpy as np
import pytest

from VegasAfterglow import (
    GaussianJet,
    ISM,
    Magnetar,
    Model,
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
)


NAN = float("nan")
INF = float("inf")


# ----------------- TophatJet -----------------

def test_tophat_happy_path():
    """TophatJet constructs with positive theta_c, positive E_iso, and Gamma0 > 1 without raising."""
    TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)


@pytest.mark.parametrize("bad", [-0.1, 0.0, NAN, INF, math.pi])
def test_tophat_bad_theta_c(bad):
    """TophatJet raises ValueError naming theta_c for each half-opening angle that is negative, zero, NaN, Inf, or as large as pi."""
    with pytest.raises(ValueError, match="theta_c"):
        TophatJet(theta_c=bad, E_iso=1e52, Gamma0=300)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_tophat_bad_E_iso(bad):
    """TophatJet raises ValueError naming E_iso for each non-positive or non-finite isotropic-equivalent energy."""
    with pytest.raises(ValueError, match="E_iso"):
        TophatJet(theta_c=0.1, E_iso=bad, Gamma0=300)


@pytest.mark.parametrize("bad", [0.0, 0.5, 1.0, NAN, INF, -10.0])
def test_tophat_bad_Gamma0(bad):
    """TophatJet raises ValueError naming Gamma0 for each initial bulk Lorentz factor that is <= 1, negative, or non-finite."""
    with pytest.raises(ValueError, match="Gamma0"):
        TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=bad)


# ----------------- GaussianJet -----------------

def test_gaussian_happy_path():
    """GaussianJet constructs with positive theta_c, positive E_iso, and Gamma0 > 1 without raising."""
    GaussianJet(theta_c=0.1, E_iso=1e52, Gamma0=300)


def test_gaussian_rejects_nan_E_iso():
    """GaussianJet raises ValueError naming E_iso when the isotropic-equivalent energy is NaN."""
    with pytest.raises(ValueError, match="E_iso"):
        GaussianJet(theta_c=0.1, E_iso=NAN, Gamma0=300)


# ----------------- PowerLawJet -----------------

def test_powerlaw_happy_path():
    """PowerLawJet constructs with valid core parameters and positive power-law indices k_e and k_g without raising."""
    PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=2.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_powerlaw_bad_k_e(bad):
    """PowerLawJet raises ValueError naming k_e for each non-positive or non-finite energy power-law index."""
    with pytest.raises(ValueError, match="k_e"):
        PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=bad, k_g=2.0)


def test_powerlaw_bad_k_g():
    """PowerLawJet raises ValueError naming k_g when the Lorentz-factor power-law index is NaN."""
    with pytest.raises(ValueError, match="k_g"):
        PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=NAN)


# ----------------- TwoComponentJet -----------------

def test_two_component_happy_path():
    """TwoComponentJet constructs when the wing angle theta_w exceeds the core angle theta_c and all core/wing parameters are valid, without raising."""
    TwoComponentJet(
        theta_c=0.05, E_iso=1e52, Gamma0=300,
        theta_w=0.3, E_iso_w=1e51, Gamma0_w=10,
    )


def test_two_component_rejects_theta_w_le_theta_c():
    """TwoComponentJet raises ValueError naming theta_w when the wing angle does not exceed the core angle theta_c."""
    # Should reject when theta_w <= theta_c
    with pytest.raises(ValueError, match="theta_w"):
        TwoComponentJet(
            theta_c=0.1, E_iso=1e52, Gamma0=300,
            theta_w=0.05, E_iso_w=1e51, Gamma0_w=10,
        )


def test_two_component_rejects_nan_E_iso_w():
    """TwoComponentJet raises ValueError naming E_iso_w when the wing isotropic-equivalent energy is NaN."""
    with pytest.raises(ValueError, match="E_iso_w"):
        TwoComponentJet(
            theta_c=0.05, E_iso=1e52, Gamma0=300,
            theta_w=0.3, E_iso_w=NAN, Gamma0_w=10,
        )


# ----------------- StepPowerLawJet -----------------

def test_step_powerlaw_happy_path():
    """StepPowerLawJet constructs with valid core and wing parameters and positive indices k_e and k_g without raising."""
    StepPowerLawJet(
        theta_c=0.05, E_iso=1e52, Gamma0=300,
        E_iso_w=1e51, Gamma0_w=10, k_e=2.0, k_g=2.0,
    )


def test_step_powerlaw_rejects_Gamma0_w_le_one():
    """StepPowerLawJet raises ValueError naming Gamma0_w when the wing bulk Lorentz factor equals 1."""
    with pytest.raises(ValueError, match="Gamma0_w"):
        StepPowerLawJet(
            theta_c=0.05, E_iso=1e52, Gamma0=300,
            E_iso_w=1e51, Gamma0_w=1.0, k_e=2.0, k_g=2.0,
        )


# ----------------- PowerLawWing -----------------

def test_powerlaw_wing_happy_path():
    """PowerLawWing constructs with positive theta_c, E_iso_w, Gamma0_w > 1, and indices k_e and k_g without raising."""
    PowerLawWing(theta_c=0.05, E_iso_w=1e51, Gamma0_w=10, k_e=2.0, k_g=2.0)


# ----------------- ISM / Wind -----------------

def test_ism_happy_path():
    """ISM constructs with a positive ambient number density n_ism without raising."""
    ISM(n_ism=1.0)


def test_ism_zero_density_allowed():
    """ISM accepts n_ism = 0 (no ambient density floor) without raising."""
    # n_ism = 0 is valid (no ambient floor); just no exception
    ISM(n_ism=0.0)


@pytest.mark.parametrize("bad", [-1.0, NAN, INF])
def test_ism_bad_n_ism(bad):
    """ISM raises ValueError naming n_ism for each negative, NaN, or Inf number density."""
    with pytest.raises(ValueError, match="n_ism"):
        ISM(n_ism=bad)


def test_wind_happy_path():
    """Wind constructs with a positive stellar-wind parameter A_star without raising."""
    Wind(A_star=1.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_wind_bad_A_star(bad):
    """Wind raises ValueError naming A_star for each non-positive or non-finite wind parameter."""
    with pytest.raises(ValueError, match="A_star"):
        Wind(A_star=bad)


def test_wind_bad_k_m():
    """Wind raises ValueError naming k_m when the density power-law slope is negative."""
    with pytest.raises(ValueError, match="k_m"):
        Wind(A_star=1.0, k_m=-1.0)


def test_wind_n0_inf_allowed():
    """Wind accepts n0 = +inf as the no-density-floor sentinel without raising."""
    # +inf is the "no floor" sentinel; should not raise
    Wind(A_star=1.0, n0=INF)


def test_wind_n0_zero_rejected():
    """Wind raises ValueError naming n0 when the density floor is zero."""
    with pytest.raises(ValueError, match="n0"):
        Wind(A_star=1.0, n0=0.0)


def test_wind_n0_negative_rejected():
    """Wind raises ValueError naming n0 when the density floor is negative."""
    with pytest.raises(ValueError, match="n0"):
        Wind(A_star=1.0, n0=-1.0)


# ----------------- Magnetar -----------------

def test_magnetar_happy_path():
    """Magnetar constructs with positive spin-down luminosity L0, timescale t0, and decay index q without raising."""
    Magnetar(L0=1e46, t0=100, q=2.0)


@pytest.mark.parametrize("bad_field,kwargs", [
    ("L0", {"L0": NAN, "t0": 100, "q": 2}),
    ("L0", {"L0": -1.0, "t0": 100, "q": 2}),
    ("L0", {"L0": 0.0, "t0": 100, "q": 2}),
    ("t0", {"L0": 1e46, "t0": NAN, "q": 2}),
    ("t0", {"L0": 1e46, "t0": -1.0, "q": 2}),
    ("q",  {"L0": 1e46, "t0": 100, "q": NAN}),
    ("q",  {"L0": 1e46, "t0": 100, "q": -1.0}),
])
def test_magnetar_rejects_bad_inputs(bad_field, kwargs):
    """Magnetar raises ValueError naming the offending parameter for each case that passes a NaN, negative, or zero value for L0, t0, or q."""
    with pytest.raises(ValueError, match=bad_field):
        Magnetar(**kwargs)


# ----------------- Observer -----------------

def test_observer_happy_path():
    """Observer constructs with positive lumi_dist, non-negative z, and theta_obs inside [0, pi] without raising."""
    Observer(lumi_dist=1e28, z=0.5, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_observer_bad_lumi_dist(bad):
    """Observer raises ValueError naming lumi_dist for each non-positive or non-finite luminosity distance."""
    with pytest.raises(ValueError, match="lumi_dist"):
        Observer(lumi_dist=bad, z=0.5, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-0.1, NAN, INF])
def test_observer_bad_z(bad):
    """Observer raises ValueError naming z for each negative, NaN, or Inf redshift."""
    with pytest.raises(ValueError, match="z"):
        Observer(lumi_dist=1e28, z=bad, theta_obs=0.1)


def test_observer_z_zero_allowed():
    """Observer accepts redshift z = 0 without raising."""
    Observer(lumi_dist=1e28, z=0.0, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-0.1, math.pi + 0.1, NAN, INF])
def test_observer_bad_theta_obs(bad):
    """Observer raises ValueError naming theta_obs for each viewing angle that is negative, greater than pi, or non-finite."""
    with pytest.raises(ValueError, match="theta_obs"):
        Observer(lumi_dist=1e28, z=0.5, theta_obs=bad)


def test_observer_theta_obs_bounds_inclusive():
    """Observer accepts the inclusive viewing-angle bounds theta_obs = 0 and theta_obs = pi without raising."""
    Observer(lumi_dist=1e28, z=0.5, theta_obs=0.0)
    Observer(lumi_dist=1e28, z=0.5, theta_obs=math.pi)


# ----------------- Radiation -----------------

def test_radiation_happy_path():
    """Radiation constructs with eps_e and eps_B in (0, 1] and electron spectral index p > 1 without raising."""
    Radiation(eps_e=0.1, eps_B=0.01, p=2.3)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5, NAN, INF])
def test_radiation_bad_eps_e(bad):
    """Radiation raises ValueError naming eps_e for each electron energy fraction that is zero, negative, above 1, or non-finite."""
    with pytest.raises(ValueError, match="eps_e"):
        Radiation(eps_e=bad, eps_B=0.01, p=2.3)


@pytest.mark.parametrize("bad", [0.0, -1e-5, 1.5, NAN, INF])
def test_radiation_bad_eps_B(bad):
    """Radiation raises ValueError naming eps_B for each magnetic energy fraction that is zero, negative, above 1, or non-finite."""
    with pytest.raises(ValueError, match="eps_B"):
        Radiation(eps_e=0.1, eps_B=bad, p=2.3)


@pytest.mark.parametrize("bad", [0.5, 1.0, -1.0, NAN, INF])
def test_radiation_bad_p(bad):
    """Radiation raises ValueError naming p for each electron spectral index that is non-finite or violates the strict p > 1 bound, including p = 1 exactly."""
    # p > 1 is the strict bound; p = 1.0 is rejected; p in (1, 2] (fast cooling) is fine.
    with pytest.raises(ValueError, match="p"):
        Radiation(eps_e=0.1, eps_B=0.01, p=bad)


def test_radiation_p_slow_cooling_allowed():
    """Radiation accepts p = 1.5, confirming hard electron spectral indices in 1 < p < 2 are not rejected."""
    # Fast-cooling regime 1 < p < 2 should not be rejected
    Radiation(eps_e=0.1, eps_B=0.01, p=1.5)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5, NAN, INF])
def test_radiation_bad_xi_e(bad):
    """Radiation raises ValueError naming xi_e for each accelerated-electron fraction that is zero, negative, above 1, or non-finite."""
    with pytest.raises(ValueError, match="xi_e"):
        Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=bad)


def test_radiation_xi_e_equals_one_allowed():
    """Radiation accepts the inclusive upper bound xi_e = 1 without raising."""
    Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=1.0)


# ----------------- Model constructor -----------------

def _good_components():
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
    medium = ISM(n_ism=1.0)
    obs = Observer(lumi_dist=1e28, z=0.5, theta_obs=0.0)
    rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.3)
    return jet, medium, obs, rad


def test_model_happy_path():
    """Model constructs from valid jet, medium, observer, and radiation components without raising."""
    Model(*_good_components())


@pytest.mark.parametrize("bad", [0.0, -1e-5, 1.0, 1.1, NAN, INF])
def test_model_bad_rtol(bad):
    """Model raises ValueError naming rtol for each relative tolerance that is zero, negative, >= 1, or non-finite."""
    jet, medium, obs, rad = _good_components()
    with pytest.raises(ValueError, match="rtol"):
        Model(jet, medium, obs, rad, rtol=bad)


@pytest.mark.parametrize("res", [(-0.1, 0.25, 10), (0.1, 0.0, 10), (0.1, 0.25, NAN), (0.1, 0.25, INF)])
def test_model_bad_resolutions(res):
    """Model raises ValueError matching 'resol' for each resolutions triple containing a negative, zero, NaN, or Inf entry."""
    jet, medium, obs, rad = _good_components()
    with pytest.raises(ValueError, match="resol"):
        Model(jet, medium, obs, rad, resolutions=res)


# ----------------- flux_density_exposures per-element -----------------

def test_flux_density_exposures_rejects_bad_expo_time():
    """flux_density_exposures raises ValueError naming the offending element expo_time[1] when one exposure time is negative."""
    model = Model(*_good_components())
    t = np.array([1e3, 1e4, 1e5])
    nu = np.array([5e14, 5e14, 5e14])
    bad_expo = np.array([100.0, -10.0, 100.0])  # negative in middle
    with pytest.raises(ValueError, match=r"expo_time\[1\]"):
        model.flux_density_exposures(t, nu, bad_expo, num_points=5)


def test_flux_density_exposures_rejects_nan_expo_time():
    """flux_density_exposures raises ValueError naming the offending element expo_time[1] when one exposure time is NaN."""
    model = Model(*_good_components())
    t = np.array([1e3, 1e4])
    nu = np.array([5e14, 5e14])
    bad_expo = np.array([100.0, NAN])
    with pytest.raises(ValueError, match=r"expo_time\[1\]"):
        model.flux_density_exposures(t, nu, bad_expo, num_points=5)


# ----------------- generic Ejecta / Medium -----------------

def test_ejecta_happy_path():
    """Ejecta constructs from callable E_iso and Gamma0 angular profiles without raising."""
    from VegasAfterglow import Ejecta
    Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0)


def test_ejecta_rejects_non_callable():
    """Ejecta raises ValueError saying 'must be callable' when E_iso or Gamma0 is passed as a scalar instead of a profile function."""
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="E_iso must be callable"):
        Ejecta(E_iso=1e52, Gamma0=lambda phi, theta: 300.0)
    with pytest.raises(ValueError, match="Gamma0 must be callable"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=300.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_ejecta_bad_duration(bad):
    """Ejecta raises ValueError naming duration for each non-positive or non-finite ejection duration."""
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="duration"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0, duration=bad)


@pytest.mark.parametrize("bad", [NAN, INF, -1e52])
def test_ejecta_rejects_broken_E_iso_profile(bad):
    """Ejecta raises ValueError naming E_iso when the energy profile callable returns NaN, Inf, or a negative value."""
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="E_iso"):
        Ejecta(E_iso=lambda phi, theta: bad, Gamma0=lambda phi, theta: 300.0)


@pytest.mark.parametrize("bad", [NAN, INF, 0.5, -2.0])
def test_ejecta_rejects_broken_Gamma0_profile(bad):
    """Ejecta raises ValueError naming Gamma0 when the Lorentz-factor profile callable returns a non-finite value or one that does not exceed 1."""
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="Gamma0"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: bad)


def test_ejecta_rejects_broken_sigma0_profile():
    """Ejecta raises ValueError naming sigma0 when the magnetization profile callable returns NaN."""
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="sigma0"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0,
               sigma0=lambda phi, theta: NAN)


def test_medium_happy_path():
    """Medium constructs from a callable mass-density profile rho(phi, theta, r) without raising."""
    from VegasAfterglow import Medium
    Medium(rho=lambda phi, theta, r: 1.67e-24)


def test_medium_rejects_non_callable():
    """Medium raises ValueError saying 'rho must be callable' when the density profile is passed as a scalar."""
    from VegasAfterglow import Medium
    with pytest.raises(ValueError, match="rho must be callable"):
        Medium(rho=1.0)


@pytest.mark.parametrize("bad", [NAN, INF, -1e-24])
def test_medium_rejects_broken_rho_profile(bad):
    """Medium raises ValueError naming rho when the density profile callable returns NaN, Inf, or a negative value."""
    from VegasAfterglow import Medium
    with pytest.raises(ValueError, match="rho"):
        Medium(rho=lambda phi, theta, r: bad)
