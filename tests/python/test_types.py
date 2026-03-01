"""
Tests for Python type definitions and basic API functionality.

These tests verify the Python interface works correctly, including
type definitions, model creation, and basic calculations.
"""

import numpy as np
import pytest

from VegasAfterglow import (
    FitResult,
    GaussianJet,
    ISM,
    Model,
    ModelParams,
    Observer,
    ParamDef,
    PowerLawJet,
    Radiation,
    Scale,
    TophatJet,
    TwoComponentJet,
    Wind,
)


class TestScaleEnum:
    """Test Scale enumeration."""

    def test_scale_values(self):
        """Test that Scale enum has expected values."""
        assert Scale.LINEAR.value == "linear"
        assert Scale.LOG.value == "log"
        assert Scale.FIXED.value == "fixed"

    def test_scale_comparison(self):
        """Test Scale enum comparison."""
        assert Scale.LINEAR == Scale.LINEAR
        assert Scale.LOG != Scale.FIXED


class TestParamDef:
    """Test ParamDef dataclass."""

    def test_paramdef_creation(self):
        """Test creating a ParamDef with all parameters."""
        param = ParamDef(
            name="E_iso",
            lower=1e48,
            upper=1e54,
            scale=Scale.LOG,
            initial=1e52
        )
        assert param.name == "E_iso"
        assert param.lower == 1e48
        assert param.upper == 1e54
        assert param.scale == Scale.LOG
        assert param.initial == 1e52

    def test_paramdef_defaults(self):
        """Test ParamDef default values."""
        param = ParamDef(name="p", lower=2.0, upper=3.0)
        assert param.scale == Scale.LINEAR
        assert param.initial is None

    def test_paramdef_fixed(self):
        """Test fixed parameter definition."""
        param = ParamDef(
            name="theta_v",
            lower=0.0,
            upper=0.0,
            scale=Scale.FIXED
        )
        assert param.scale == Scale.FIXED


class TestJetCreation:
    """Test jet structure creation."""

    def test_tophat_jet(self):
        """Test TophatJet creation."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        # Jet object should be created without error
        assert jet is not None

    def test_gaussian_jet(self):
        """Test GaussianJet creation."""
        jet = GaussianJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        assert jet is not None

    def test_powerlaw_jet(self):
        """Test PowerLawJet creation."""
        jet = PowerLawJet(
            theta_c=0.1, E_iso=1e52, Gamma0=300,
            k_e=2.0, k_g=2.0
        )
        assert jet is not None

    def test_two_component_jet(self):
        """Test TwoComponentJet creation."""
        jet = TwoComponentJet(
            theta_c=0.05, E_iso=1e52, Gamma0=300,
            theta_w=0.2, E_iso_w=1e50, Gamma0_w=50
        )
        assert jet is not None

    def test_jet_with_spreading(self):
        """Test jet with spreading enabled."""
        jet = TophatJet(
            theta_c=0.1, E_iso=1e52, Gamma0=300,
            spreading=True
        )
        assert jet is not None


class TestMediumCreation:
    """Test ambient medium creation."""

    def test_ism_creation(self):
        """Test ISM medium creation."""
        medium = ISM(n_ism=1.0)
        assert medium is not None

    def test_wind_creation(self):
        """Test Wind medium creation."""
        medium = Wind(A_star=0.1)
        assert medium is not None

    def test_wind_with_params(self):
        """Test Wind with additional parameters."""
        medium = Wind(A_star=0.1, n_ism=0.01, n0=1e6, k_m=2)
        assert medium is not None


class TestObserverCreation:
    """Test observer creation."""

    def test_observer_on_axis(self):
        """Test on-axis observer."""
        observer = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        assert observer is not None

    def test_observer_off_axis(self):
        """Test off-axis observer."""
        observer = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.3)
        assert observer is not None


class TestRadiationCreation:
    """Test radiation settings creation."""

    def test_basic_radiation(self):
        """Test basic synchrotron radiation settings."""
        radiation = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, xi_e=1.0)
        assert radiation is not None

    def test_radiation_with_ssc(self):
        """Test radiation with SSC enabled."""
        radiation = Radiation(
            eps_e=0.1, eps_B=0.01, p=2.2, xi_e=1.0,
            ssc=True, kn=True
        )
        assert radiation is not None


class TestModelCreation:
    """Test full model creation."""

    @pytest.fixture
    def basic_components(self):
        """Create basic model components."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        observer = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        radiation = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, xi_e=1.0)
        return jet, medium, observer, radiation

    def test_model_creation(self, basic_components):
        """Test basic model creation."""
        jet, medium, observer, radiation = basic_components
        model = Model(jet, medium, observer, radiation)
        assert model is not None

    def test_model_with_resolution(self, basic_components):
        """Test model with custom resolution."""
        jet, medium, observer, radiation = basic_components
        model = Model(
            jet, medium, observer, radiation,
            resolutions=(0.5, 3, 20)
        )
        assert model is not None


class TestModelCalculations:
    """Test model calculation methods."""

    @pytest.fixture
    def model(self):
        """Create a basic model for testing."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        observer = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        radiation = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, xi_e=1.0)
        return Model(jet, medium, observer, radiation, resolutions=(0.3, 1, 10))

    def test_flux_density(self, model, sample_time_array, sample_frequency):
        """Test flux_density calculation."""
        t = sample_time_array
        nu = np.full_like(t, sample_frequency)

        flux = model.flux_density(t, nu)

        # Check output structure
        assert hasattr(flux, 'total')
        assert hasattr(flux, 'fwd')
        assert hasattr(flux, 'rvs')

        # Check output shape
        assert flux.total.shape == t.shape

        # Check values are positive
        assert np.all(flux.total > 0)

    def test_flux_density_grid(self, model, sample_time_array, sample_frequency_array):
        """Test flux_density_grid calculation."""
        t = sample_time_array[:10]  # Smaller for speed
        nu = sample_frequency_array[:5]

        flux = model.flux_density_grid(t, nu)

        # Check output shape (should be 2D: len(nu) x len(t))
        assert hasattr(flux, 'total')
        assert flux.total.shape == (len(nu), len(t))

    def test_details(self, model):
        """Test details() method returns shock evolution data."""
        details = model.details(t_min=1e2, t_max=1e5)

        # Check structure
        assert hasattr(details, 'fwd')
        assert hasattr(details.fwd, 'r')
        assert hasattr(details.fwd, 'Gamma')
        assert hasattr(details.fwd, 'B_comv')
        assert hasattr(details.fwd, 'gamma_m')
        assert hasattr(details.fwd, 'gamma_c')
        assert hasattr(details.fwd, 'nu_m')
        assert hasattr(details.fwd, 'nu_c')

    def test_jet_profile(self, model):
        """Test jet profile methods."""
        theta = np.linspace(0, 0.5, 20)

        E_iso = model.jet_E_iso(0.0, theta)
        Gamma0 = model.jet_Gamma0(0.0, theta)

        assert E_iso.shape == theta.shape
        assert Gamma0.shape == theta.shape
        assert np.all(E_iso >= 0)
        assert np.all(Gamma0 >= 1)


class TestModelParams:
    """Test ModelParams class."""

    def test_params_creation(self):
        """Test creating ModelParams."""
        params = ModelParams()
        assert params is not None

    def test_params_attributes(self):
        """Test setting ModelParams attributes."""
        params = ModelParams()

        # Set various parameters
        params.E_iso = 1e52
        params.Gamma0 = 300
        params.theta_c = 0.1
        params.p = 2.2
        params.eps_e = 0.1
        params.eps_B = 0.01

        assert params.E_iso == 1e52
        assert params.Gamma0 == 300


class TestFitResult:
    """Test FitResult dataclass."""

    def test_fitresult_creation(self):
        """Test creating FitResult."""
        samples = np.random.randn(100, 5)
        log_probs = np.random.randn(100)
        labels = ["E_iso", "Gamma0", "theta_c", "p", "eps_e"]

        result = FitResult(
            samples=samples,
            log_probs=log_probs,
            labels=labels
        )

        assert result.samples.shape == (100, 5)
        assert result.log_probs.shape == (100,)
        assert len(result.labels) == 5

    def test_fitresult_with_topk(self):
        """Test FitResult with top-k parameters."""
        samples = np.random.randn(100, 5)
        log_probs = np.random.randn(100)
        labels = ["E_iso", "Gamma0", "theta_c", "p", "eps_e"]
        top_k = samples[:10]
        top_k_probs = log_probs[:10]

        result = FitResult(
            samples=samples,
            log_probs=log_probs,
            labels=labels,
            top_k_params=top_k,
            top_k_log_probs=top_k_probs
        )

        assert result.top_k_params.shape == (10, 5)


class TestInputValidation:
    """Test input validation and error handling."""

    @pytest.fixture
    def model(self):
        """Create a model for testing."""
        jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)
        medium = ISM(n_ism=1.0)
        observer = Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0)
        radiation = Radiation(eps_e=0.1, eps_B=0.01, p=2.2, xi_e=1.0)
        return Model(jet, medium, observer, radiation, resolutions=(0.3, 1, 10))

    def test_flux_density_requires_ascending_time(self, model):
        """Test that flux_density requires ascending time array."""
        t = np.array([1e5, 1e4, 1e3])  # Descending - should fail
        nu = np.full_like(t, 4.84e14)

        with pytest.raises(Exception):
            model.flux_density(t, nu)

    def test_flux_density_requires_matching_shapes(self, model):
        """Test that flux_density requires matching array shapes."""
        t = np.logspace(2, 6, 10)
        nu = np.logspace(9, 18, 5)  # Different length

        with pytest.raises(Exception):
            model.flux_density(t, nu)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
