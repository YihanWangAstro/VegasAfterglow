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
    TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)


@pytest.mark.parametrize("bad", [-0.1, 0.0, NAN, INF, math.pi])
def test_tophat_bad_theta_c(bad):
    with pytest.raises(ValueError, match="theta_c"):
        TophatJet(theta_c=bad, E_iso=1e52, Gamma0=300)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_tophat_bad_E_iso(bad):
    with pytest.raises(ValueError, match="E_iso"):
        TophatJet(theta_c=0.1, E_iso=bad, Gamma0=300)


@pytest.mark.parametrize("bad", [0.0, 0.5, 1.0, NAN, INF, -10.0])
def test_tophat_bad_Gamma0(bad):
    with pytest.raises(ValueError, match="Gamma0"):
        TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=bad)


# ----------------- GaussianJet -----------------

def test_gaussian_happy_path():
    GaussianJet(theta_c=0.1, E_iso=1e52, Gamma0=300)


def test_gaussian_rejects_nan_E_iso():
    with pytest.raises(ValueError, match="E_iso"):
        GaussianJet(theta_c=0.1, E_iso=NAN, Gamma0=300)


# ----------------- PowerLawJet -----------------

def test_powerlaw_happy_path():
    PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=2.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_powerlaw_bad_k_e(bad):
    with pytest.raises(ValueError, match="k_e"):
        PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=bad, k_g=2.0)


def test_powerlaw_bad_k_g():
    with pytest.raises(ValueError, match="k_g"):
        PowerLawJet(theta_c=0.1, E_iso=1e52, Gamma0=300, k_e=2.0, k_g=NAN)


# ----------------- TwoComponentJet -----------------

def test_two_component_happy_path():
    TwoComponentJet(
        theta_c=0.05, E_iso=1e52, Gamma0=300,
        theta_w=0.3, E_iso_w=1e51, Gamma0_w=10,
    )


def test_two_component_rejects_theta_w_le_theta_c():
    # Should reject when theta_w <= theta_c
    with pytest.raises(ValueError, match="theta_w"):
        TwoComponentJet(
            theta_c=0.1, E_iso=1e52, Gamma0=300,
            theta_w=0.05, E_iso_w=1e51, Gamma0_w=10,
        )


def test_two_component_rejects_nan_E_iso_w():
    with pytest.raises(ValueError, match="E_iso_w"):
        TwoComponentJet(
            theta_c=0.05, E_iso=1e52, Gamma0=300,
            theta_w=0.3, E_iso_w=NAN, Gamma0_w=10,
        )


# ----------------- StepPowerLawJet -----------------

def test_step_powerlaw_happy_path():
    StepPowerLawJet(
        theta_c=0.05, E_iso=1e52, Gamma0=300,
        E_iso_w=1e51, Gamma0_w=10, k_e=2.0, k_g=2.0,
    )


def test_step_powerlaw_rejects_Gamma0_w_le_one():
    with pytest.raises(ValueError, match="Gamma0_w"):
        StepPowerLawJet(
            theta_c=0.05, E_iso=1e52, Gamma0=300,
            E_iso_w=1e51, Gamma0_w=1.0, k_e=2.0, k_g=2.0,
        )


# ----------------- PowerLawWing -----------------

def test_powerlaw_wing_happy_path():
    PowerLawWing(theta_c=0.05, E_iso_w=1e51, Gamma0_w=10, k_e=2.0, k_g=2.0)


# ----------------- ISM / Wind -----------------

def test_ism_happy_path():
    ISM(n_ism=1.0)


def test_ism_zero_density_allowed():
    # n_ism = 0 is valid (no ambient floor); just no exception
    ISM(n_ism=0.0)


@pytest.mark.parametrize("bad", [-1.0, NAN, INF])
def test_ism_bad_n_ism(bad):
    with pytest.raises(ValueError, match="n_ism"):
        ISM(n_ism=bad)


def test_wind_happy_path():
    Wind(A_star=1.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_wind_bad_A_star(bad):
    with pytest.raises(ValueError, match="A_star"):
        Wind(A_star=bad)


def test_wind_bad_k_m():
    with pytest.raises(ValueError, match="k_m"):
        Wind(A_star=1.0, k_m=-1.0)


def test_wind_n0_inf_allowed():
    # +inf is the "no floor" sentinel; should not raise
    Wind(A_star=1.0, n0=INF)


def test_wind_n0_zero_rejected():
    with pytest.raises(ValueError, match="n0"):
        Wind(A_star=1.0, n0=0.0)


def test_wind_n0_negative_rejected():
    with pytest.raises(ValueError, match="n0"):
        Wind(A_star=1.0, n0=-1.0)


# ----------------- Magnetar -----------------

def test_magnetar_happy_path():
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
    with pytest.raises(ValueError, match=bad_field):
        Magnetar(**kwargs)


# ----------------- Observer -----------------

def test_observer_happy_path():
    Observer(lumi_dist=1e28, z=0.5, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_observer_bad_lumi_dist(bad):
    with pytest.raises(ValueError, match="lumi_dist"):
        Observer(lumi_dist=bad, z=0.5, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-0.1, NAN, INF])
def test_observer_bad_z(bad):
    with pytest.raises(ValueError, match="z"):
        Observer(lumi_dist=1e28, z=bad, theta_obs=0.1)


def test_observer_z_zero_allowed():
    Observer(lumi_dist=1e28, z=0.0, theta_obs=0.1)


@pytest.mark.parametrize("bad", [-0.1, math.pi + 0.1, NAN, INF])
def test_observer_bad_theta_obs(bad):
    with pytest.raises(ValueError, match="theta_obs"):
        Observer(lumi_dist=1e28, z=0.5, theta_obs=bad)


def test_observer_theta_obs_bounds_inclusive():
    Observer(lumi_dist=1e28, z=0.5, theta_obs=0.0)
    Observer(lumi_dist=1e28, z=0.5, theta_obs=math.pi)


# ----------------- Radiation -----------------

def test_radiation_happy_path():
    Radiation(eps_e=0.1, eps_B=0.01, p=2.3)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5, NAN, INF])
def test_radiation_bad_eps_e(bad):
    with pytest.raises(ValueError, match="eps_e"):
        Radiation(eps_e=bad, eps_B=0.01, p=2.3)


@pytest.mark.parametrize("bad", [0.0, -1e-5, 1.5, NAN, INF])
def test_radiation_bad_eps_B(bad):
    with pytest.raises(ValueError, match="eps_B"):
        Radiation(eps_e=0.1, eps_B=bad, p=2.3)


@pytest.mark.parametrize("bad", [0.5, 1.0, -1.0, NAN, INF])
def test_radiation_bad_p(bad):
    # p > 1 is the strict bound; p = 1.0 is rejected; p in (1, 2] (fast cooling) is fine.
    with pytest.raises(ValueError, match="p"):
        Radiation(eps_e=0.1, eps_B=0.01, p=bad)


def test_radiation_p_slow_cooling_allowed():
    # Fast-cooling regime 1 < p < 2 should not be rejected
    Radiation(eps_e=0.1, eps_B=0.01, p=1.5)


@pytest.mark.parametrize("bad", [0.0, -0.1, 1.5, NAN, INF])
def test_radiation_bad_xi_e(bad):
    with pytest.raises(ValueError, match="xi_e"):
        Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=bad)


def test_radiation_xi_e_equals_one_allowed():
    Radiation(eps_e=0.1, eps_B=0.01, p=2.3, xi_e=1.0)


# ----------------- Model constructor -----------------

def _good_components():
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
    medium = ISM(n_ism=1.0)
    obs = Observer(lumi_dist=1e28, z=0.5, theta_obs=0.0)
    rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.3)
    return jet, medium, obs, rad


def test_model_happy_path():
    Model(*_good_components())


@pytest.mark.parametrize("bad", [0.0, -1e-5, 1.0, 1.1, NAN, INF])
def test_model_bad_rtol(bad):
    jet, medium, obs, rad = _good_components()
    with pytest.raises(ValueError, match="rtol"):
        Model(jet, medium, obs, rad, rtol=bad)


@pytest.mark.parametrize("res", [(-0.1, 0.25, 10), (0.1, 0.0, 10), (0.1, 0.25, NAN), (0.1, 0.25, INF)])
def test_model_bad_resolutions(res):
    jet, medium, obs, rad = _good_components()
    with pytest.raises(ValueError, match="resol"):
        Model(jet, medium, obs, rad, resolutions=res)


# ----------------- flux_density_exposures per-element -----------------

def test_flux_density_exposures_rejects_bad_expo_time():
    model = Model(*_good_components())
    t = np.array([1e3, 1e4, 1e5])
    nu = np.array([5e14, 5e14, 5e14])
    bad_expo = np.array([100.0, -10.0, 100.0])  # negative in middle
    with pytest.raises(ValueError, match=r"expo_time\[1\]"):
        model.flux_density_exposures(t, nu, bad_expo, num_points=5)


def test_flux_density_exposures_rejects_nan_expo_time():
    model = Model(*_good_components())
    t = np.array([1e3, 1e4])
    nu = np.array([5e14, 5e14])
    bad_expo = np.array([100.0, NAN])
    with pytest.raises(ValueError, match=r"expo_time\[1\]"):
        model.flux_density_exposures(t, nu, bad_expo, num_points=5)


# ----------------- generic Ejecta / Medium -----------------

def test_ejecta_happy_path():
    from VegasAfterglow import Ejecta
    Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0)


def test_ejecta_rejects_non_callable():
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="E_iso must be callable"):
        Ejecta(E_iso=1e52, Gamma0=lambda phi, theta: 300.0)
    with pytest.raises(ValueError, match="Gamma0 must be callable"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=300.0)


@pytest.mark.parametrize("bad", [-1.0, 0.0, NAN, INF])
def test_ejecta_bad_duration(bad):
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="duration"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0, duration=bad)


@pytest.mark.parametrize("bad", [NAN, INF, -1e52])
def test_ejecta_rejects_broken_E_iso_profile(bad):
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="E_iso"):
        Ejecta(E_iso=lambda phi, theta: bad, Gamma0=lambda phi, theta: 300.0)


@pytest.mark.parametrize("bad", [NAN, INF, 0.5, -2.0])
def test_ejecta_rejects_broken_Gamma0_profile(bad):
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="Gamma0"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: bad)


def test_ejecta_rejects_broken_sigma0_profile():
    from VegasAfterglow import Ejecta
    with pytest.raises(ValueError, match="sigma0"):
        Ejecta(E_iso=lambda phi, theta: 1e52, Gamma0=lambda phi, theta: 300.0,
               sigma0=lambda phi, theta: NAN)


def test_medium_happy_path():
    from VegasAfterglow import Medium
    Medium(rho=lambda phi, theta, r: 1.67e-24)


def test_medium_rejects_non_callable():
    from VegasAfterglow import Medium
    with pytest.raises(ValueError, match="rho must be callable"):
        Medium(rho=1.0)


@pytest.mark.parametrize("bad", [NAN, INF, -1e-24])
def test_medium_rejects_broken_rho_profile(bad):
    from VegasAfterglow import Medium
    with pytest.raises(ValueError, match="rho"):
        Medium(rho=lambda phi, theta, r: bad)
