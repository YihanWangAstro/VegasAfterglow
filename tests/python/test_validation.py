"""
Input validation tests for VegasAfterglow v2.0.0.
Verifies that invalid inputs raise appropriate errors.
"""

import numpy as np
import pytest

from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet

LOW_RES = (0.3, 1, 10)


@pytest.fixture
def model():
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
    medium = ISM(n_ism=1.0)
    obs = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
    rad = Radiation(eps_e=0.1, eps_B=0.01, p=2.2)
    return Model(jet, medium, obs, rad, resolutions=LOW_RES)


class TestFluxDensityValidation:
    """Test input validation for flux_density."""

    def test_descending_time_raises(self, model):
        """flux_density raises an exception when the time array is in descending order instead of ascending."""
        t = np.array([1e5, 1e4, 1e3])
        nu = np.full_like(t, 4.84e14)
        with pytest.raises(Exception):
            model.flux_density(t, nu)

    def test_mismatched_shapes_raises(self, model):
        """flux_density raises an exception when the time and frequency arrays have different lengths."""
        t = np.logspace(2, 6, 10)
        nu = np.logspace(9, 18, 5)
        with pytest.raises(Exception):
            model.flux_density(t, nu)

    def test_empty_arrays_raises(self, model):
        """flux_density raises an exception when both the time and frequency arrays are empty."""
        with pytest.raises(Exception):
            model.flux_density(np.array([]), np.array([]))


class TestFluxDensityGridValidation:
    """Test input validation for flux_density_grid."""

    def test_empty_time_raises(self, model):
        """flux_density_grid raises an exception when the time array is empty even though the frequency array is non-empty."""
        with pytest.raises(Exception):
            model.flux_density_grid(np.array([]), np.array([1e14]))

    def test_empty_freq_raises(self, model):
        """flux_density_grid raises an exception when the frequency array is empty even though the time array is non-empty."""
        with pytest.raises(Exception):
            model.flux_density_grid(np.array([1e4]), np.array([]))


class TestFluxValidation:
    """Test input validation for flux (bolometric)."""

    def test_empty_time_raises(self, model):
        """flux raises an exception when the time array is empty despite a valid frequency integration range."""
        with pytest.raises(Exception):
            model.flux(np.array([]), 1e9, 1e18, 10)

    def test_negative_nu_min_raises(self, model):
        """flux raises an exception when the lower frequency bound of the integration band is negative."""
        with pytest.raises(Exception):
            model.flux(np.array([1e4]), -1e9, 1e18, 10)

    def test_nu_max_less_than_nu_min_raises(self, model):
        """flux raises an exception when the upper frequency bound of the integration band is below the lower bound."""
        with pytest.raises(Exception):
            model.flux(np.array([1e4]), 1e18, 1e9, 10)

    def test_num_nu_one_raises(self, model):
        """flux raises an exception when only one frequency sample point is requested for the band integration."""
        with pytest.raises(Exception):
            model.flux(np.array([1e4]), 1e9, 1e18, 1)


class TestExposureValidation:
    """Test input validation for flux_density_exposures."""

    def test_mismatched_shapes_raises(self, model):
        """flux_density_exposures raises an exception when the exposure array is shorter than the time and frequency arrays."""
        t = np.logspace(3, 5, 5)
        nu = np.full_like(t, 4.84e14)
        expo = np.full(3, 100.0)
        with pytest.raises(Exception):
            model.flux_density_exposures(t, nu, expo)

    def test_num_points_too_small_raises(self, model):
        """flux_density_exposures raises an exception when num_points is 1, too few samples to average over an exposure window."""
        t = np.logspace(3, 5, 3)
        nu = np.full_like(t, 4.84e14)
        expo = np.full_like(t, 100.0)
        with pytest.raises(Exception):
            model.flux_density_exposures(t, nu, expo, num_points=1)
