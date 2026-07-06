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
        return Model(jet, medium, obs, rad)

    def test_exposure_produces_valid_flux(self, model):
        """Exposure-averaged flux has the same shape as the input times and is finite and positive."""
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
        """A jet with magnetar spin-down energy injection produces finite, positive flux density."""
        mag = Magnetar(L0=1e47, t0=1e4, q=2)
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300, magnetar=mag)
        medium = ISM(n_ism=1.0)
        obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
        model = Model(jet, medium, obs, rad)
        t = np.logspace(3, 6, 10)
        nu = np.full_like(t, 4.84e14)
        flux = model.flux_density(t, nu)
        assert np.all(np.isfinite(flux.total))
        assert np.all(flux.total > 0)

    def test_magnetar_repr(self):
        """repr() of a Magnetar instance contains the class name."""
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
        return Model(jet, medium, obs, rad)

    def test_observer_property(self, model):
        """Model.observer round-trips the redshift and viewing angle given at construction."""
        obs = model.observer
        assert obs.z == 1.0
        assert obs.theta_obs == 0.0

    def test_fwd_rad_property(self, model):
        """Model.fwd_rad round-trips the eps_e and p microphysics parameters given at construction."""
        rad = model.fwd_rad
        assert rad.eps_e == pytest.approx(0.1)
        assert rad.p == pytest.approx(2.2)

    def test_rvs_rad_property(self, model):
        """Model.rvs_rad is None when no reverse-shock radiation is configured."""
        assert model.rvs_rad is None

    def test_resolutions_property(self, model):
        """Model.resolutions returns a length-3 tuple of grid resolutions."""
        r = model.resolutions
        assert len(r) == 3

    def test_rtol_property(self, model):
        """Model.rtol exposes a positive solver relative tolerance."""
        assert model.rtol > 0

    def test_axisymmetric_property(self, model):
        """Model.axisymmetric is True for an on-axis tophat jet configuration."""
        assert model.axisymmetric is True

    def test_repr(self, model):
        """repr() of a Model instance contains the class name."""
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
        return Model(jet, medium, obs, rad)

    def test_single_time_point(self, model):
        """flux_density with a single time-frequency pair returns a length-1 finite array."""
        t = np.array([1e4])
        nu = np.array([4.84e14])
        flux = model.flux_density(t, nu)
        assert flux.total.shape == (1,)
        assert np.isfinite(flux.total[0])

    def test_flux_density_grid_single_freq(self, model):
        """flux_density_grid with one frequency returns a (1, n_times) array."""
        t = np.logspace(3, 5, 5)
        nu = np.array([4.84e14])
        flux = model.flux_density_grid(t, nu)
        assert flux.total.shape == (1, 5)

    def test_flux_density_grid_single_time(self, model):
        """flux_density_grid with one time returns an (n_frequencies, 1) array."""
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
        model = Model(jet, medium, obs, rad)
        det = model.details(t_min=1e3, t_max=1e5)
        assert det.fwd.ssc_spectrum is not None


# ========================================
#  logscale_screen Utility
# ========================================
class TestLogscaleScreen:
    """Test the logscale_screen utility function."""

    def test_basic_screening(self):
        """logscale_screen returns a non-empty proper subset of indices for log-spaced data."""
        data = np.logspace(1, 5, 1000)
        indices = logscale_screen(data, 10)
        assert len(indices) > 0
        assert len(indices) < 1000

    def test_single_element(self):
        """logscale_screen on a single-element array returns exactly index 0."""
        data = np.array([42.0])
        indices = logscale_screen(data, 10)
        assert len(indices) == 1
        assert indices[0] == 0
