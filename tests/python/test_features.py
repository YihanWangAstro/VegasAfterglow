"""
Feature coverage tests for VegasAfterglow v2.0.0.
Tests verify API contracts and output sanity (shapes, finiteness, positivity),
not physics accuracy.
"""

import numpy as np
import pytest

from VegasAfterglow import (
    Ejecta,
    GaussianJet,
    ISM,
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

LOW_RES = (0.3, 1, 10)


# ========================================
#  SSC / Inverse Compton
# ========================================
class TestSSC:
    """Test synchrotron self-Compton emission."""

    @pytest.fixture
    def ssc_model(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, ssc=True)
        return Model(jet, medium, obs, rad, resolutions=LOW_RES)

    def test_ssc_flux_has_components(self, ssc_model):
        """SSC model should produce both sync and ssc flux components."""
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 1e14)
        flux = ssc_model.flux_density(t, nu)
        assert flux.fwd.sync.shape == t.shape
        assert flux.fwd.ssc.shape == t.shape
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_ssc_total_ge_sync(self, ssc_model):
        """Total flux should be >= sync-only flux."""
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 1e14)
        flux = ssc_model.flux_density(t, nu)
        assert np.all(flux.total >= flux.fwd.sync - 1e-30)


# ========================================
#  Klein-Nishina
# ========================================
class TestKN:
    """Test Klein-Nishina corrections."""

    def test_kn_model_runs(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, ssc=True, kn=True)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 1e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)


# ========================================
#  Reverse Shock
# ========================================
class TestReverseShock:
    """Test reverse shock emission."""

    def test_reverse_shock_produces_flux(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        fwd_rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        rvs_rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, fwd_rad, rvs_rad=rvs_rad, resolutions=LOW_RES)
        t = np.logspace(2, 6, 15)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert flux.rvs.sync.shape == t.shape
        assert np.all(np.isfinite(flux.total))
        assert np.any(flux.rvs.sync > 0)

    def test_reverse_shock_details(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        fwd_rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        rvs_rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, fwd_rad, rvs_rad=rvs_rad, resolutions=LOW_RES)
        det = model.details(t_min=1e2, t_max=1e5)
        assert det.rvs.Gamma.size > 0


# ========================================
#  Wind Medium
# ========================================
class TestWindFlux:
    """Test Wind medium produces valid flux."""

    def test_wind_medium_flux(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = Wind(A_star=0.1)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(2, 6, 15)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_wind_with_floor(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = Wind(A_star=0.1, n_ism=0.01)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(2, 6, 10)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)


# ========================================
#  Off-Axis Observer
# ========================================
class TestOffAxis:
    """Test off-axis observer produces valid flux."""

    def test_off_axis_flux(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.3)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(3, 7, 15)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_off_axis_dimmer_at_early_times(self):
        """Off-axis observer should see less flux at early times."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        on_axis = Model(
            jet, medium, Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0),
            rad, resolutions=LOW_RES,
        )
        off_axis = Model(
            jet, medium, Observer(lumi_dist=1e28, z=1.0, theta_obs=0.4),
            rad, resolutions=LOW_RES,
        )
        t = np.logspace(3, 5, 5)
        nu = np.full_like(t, 4.84e14)
        f_on = on_axis.flux_density(t, nu)
        f_off = off_axis.flux_density(t, nu)
        assert f_on.total[0] > f_off.total[0]


# ========================================
#  Missing Jet Types
# ========================================
class TestJetTypeFlux:
    """Test untested jet types produce valid flux."""

    @pytest.fixture(params=["powerlaw_wing", "step_powerlaw"])
    def model(self, request):
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        if request.param == "powerlaw_wing":
            jet = PowerLawWing(
                theta_c=0.05, E_iso_w=1e50, Gamma0_w=50, k_e=2.0, k_g=2.0,
            )
        elif request.param == "step_powerlaw":
            jet = StepPowerLawJet(
                theta_c=0.05, E_iso=1e52, Gamma0=300,
                E_iso_w=1e50, Gamma0_w=50, k_e=2.0, k_g=2.0,
            )
        return Model(jet, medium, obs, rad, resolutions=LOW_RES)

    def test_produces_finite_positive_flux(self, model):
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)


# ========================================
#  CMB Cooling
# ========================================
class TestCMBCooling:
    """Test CMB inverse Compton cooling."""

    def test_cmb_cooling_runs(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=2.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, cmb_cooling=True)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)
