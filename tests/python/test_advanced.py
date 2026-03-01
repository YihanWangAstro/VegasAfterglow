"""
Advanced feature and edge case tests for VegasAfterglow v2.0.0.
"""

import numpy as np
import pytest

from VegasAfterglow import (
    ISM,
    Magnetar,
    Model,
    Observer,
    Radiation,
    TophatJet,
    logscale_screen,
)

LOW_RES = (0.3, 1, 10)


# ========================================
#  Exposure Time Averaging
# ========================================
class TestExposureAveraging:
    """Test flux_density_exposures method."""

    @pytest.fixture
    def model(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        return Model(jet, medium, obs, rad, resolutions=LOW_RES)

    def test_exposure_produces_valid_flux(self, model):
        t = np.logspace(3, 6, 5)
        nu = np.full_like(t, 4.84e14)
        expo = np.full_like(t, 100.0)
        flux = model.flux_density_exposures(t, nu, expo, num_points=5)
        assert flux.total.shape == t.shape
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_exposure_close_to_instantaneous(self, model):
        """Very short exposure should match instantaneous flux."""
        t = np.logspace(4, 5, 3)
        nu = np.full_like(t, 4.84e14)
        flux_inst = model.flux_density(t, nu)
        expo = np.full_like(t, 0.01)
        flux_expo = model.flux_density_exposures(t, nu, expo, num_points=3)
        np.testing.assert_allclose(flux_expo.total, flux_inst.total, rtol=0.1)


# ========================================
#  Magnetar Energy Injection
# ========================================
class TestMagnetar:
    """Test Magnetar energy injection."""

    def test_magnetar_model_runs(self):
        mag = Magnetar(L0=1e47, t0=1e4, q=2)
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300, magnetar=mag)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_magnetar_repr(self):
        mag = Magnetar(L0=1e47, t0=1e4, q=2)
        r = repr(mag)
        assert "Magnetar" in r


# ========================================
#  Model Property Accessors
# ========================================
class TestModelProperties:
    """Test Model read-only property accessors."""

    @pytest.fixture
    def model(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        return Model(jet, medium, obs, rad, resolutions=LOW_RES)

    def test_observer_property(self, model):
        obs = model.observer
        assert obs.z == 1.0
        assert obs.theta_obs == 0.0

    def test_fwd_rad_property(self, model):
        rad = model.fwd_rad
        assert rad.eps_e == pytest.approx(0.1)
        assert rad.p == pytest.approx(2.2)

    def test_rvs_rad_property(self, model):
        assert model.rvs_rad is None

    def test_resolutions_property(self, model):
        r = model.resolutions
        assert len(r) == 3

    def test_rtol_property(self, model):
        assert model.rtol > 0

    def test_axisymmetric_property(self, model):
        assert model.axisymmetric is True

    def test_repr(self, model):
        r = repr(model)
        assert "Model" in r


# ========================================
#  Edge Cases
# ========================================
class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    @pytest.fixture
    def model(self):
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        return Model(jet, medium, obs, rad, resolutions=LOW_RES)

    def test_single_time_point(self, model):
        t = np.array([1e4])
        nu = np.array([4.84e14])
        flux = model.flux_density(t, nu)
        assert flux.total.shape == (1,)
        assert np.isfinite(flux.total[0])

    def test_flux_density_grid_single_freq(self, model):
        t = np.logspace(3, 5, 5)
        nu = np.array([4.84e14])
        flux = model.flux_density_grid(t, nu)
        assert flux.total.shape == (1, 5)

    def test_flux_density_grid_single_time(self, model):
        t = np.array([1e4])
        nu = np.logspace(9, 18, 5)
        flux = model.flux_density_grid(t, nu)
        assert flux.total.shape == (5, 1)

    def test_details_with_ssc(self):
        """Test details() with SSC enabled returns ssc_spectrum."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, ssc=True)
        model = Model(jet, medium, obs, rad, resolutions=LOW_RES)
        det = model.details(t_min=1e3, t_max=1e5)
        assert det.fwd.ssc_spectrum is not None


# ========================================
#  logscale_screen Utility
# ========================================
class TestLogscaleScreen:
    """Test the logscale_screen utility function."""

    def test_basic_screening(self):
        data = np.logspace(1, 5, 1000)
        indices = logscale_screen(data, 10)
        assert len(indices) > 0
        assert len(indices) < 1000

    def test_single_element(self):
        data = np.array([42.0])
        indices = logscale_screen(data, 10)
        assert len(indices) == 1
        assert indices[0] == 0
