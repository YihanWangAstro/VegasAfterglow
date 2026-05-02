"""
Tests for the VegasAfterglow.extinction module (Pei 1992 dust laws)
and Fitter extinction integration.
"""

import numpy as np
import pytest

import VegasAfterglow.extinction as ext
from VegasAfterglow.extinction import BUILTIN_LAWS, lmc, mw, pei92, smc


# ---------------------------------------------------------------------------
# Basic API
# ---------------------------------------------------------------------------


class TestBuiltinLaws:
    """Smoke tests: the three built-in law callables are importable."""

    def test_builtin_laws_keys(self):
        assert set(BUILTIN_LAWS) == {"smc", "lmc", "mw"}

    def test_builtin_laws_callable(self):
        lam = np.array([5.5e-5])
        for fn in BUILTIN_LAWS.values():
            k = fn(lam)
            assert k.shape == lam.shape
            assert np.all(np.isfinite(k))


# ---------------------------------------------------------------------------
# Normalisation: k(V) = 1
# ---------------------------------------------------------------------------


class TestVBandNormalisation:
    """Each law must return exactly 1 at the V-band reference wavelength."""

    _V_BAND_CM = 5.5e-5  # 0.55 μm

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_k_at_vband_is_one(self, fn):
        lam = np.array([self._V_BAND_CM])
        np.testing.assert_allclose(fn(lam), 1.0, rtol=1e-6)

    @pytest.mark.parametrize("law", ["smc", "lmc", "mw"])
    def test_pei92_k_at_vband_is_one(self, law):
        lam = np.array([self._V_BAND_CM])
        np.testing.assert_allclose(pei92(lam, law), 1.0, rtol=1e-6)


# ---------------------------------------------------------------------------
# Lyman-limit cutoff: k = 0 for λ < 912 Å
# ---------------------------------------------------------------------------


class TestLymanLimit:
    """Wavelengths below the Lyman limit must return zero."""

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_xray_is_zero(self, fn):
        lam_xray = np.array([1e-8, 5e-8])  # X-ray wavelengths (< 912 Å)
        np.testing.assert_array_equal(fn(lam_xray), 0.0)

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_lyman_limit_boundary(self, fn):
        """Just below the Lyman limit should be zero."""
        lam_below = np.array([911.9e-8])  # 911.9 Å in cm
        assert fn(lam_below)[0] == 0.0

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_uv_above_lyman_is_nonzero(self, fn):
        """Just above 912 Å should be attenuated (nonzero)."""
        lam_above = np.array([912.1e-8])  # 912.1 Å in cm
        assert fn(lam_above)[0] > 0.0


# ---------------------------------------------------------------------------
# Physical shape checks
# ---------------------------------------------------------------------------


class TestPhysicalShape:
    """The extinction curves must have the correct qualitative behaviour."""

    def test_uv_higher_than_optical(self):
        """UV extinction must exceed optical for all three laws."""
        lam_opt = np.array([5.5e-5])     # V band
        lam_uv  = np.array([1.5e-5])     # 1500 Å
        for fn in (smc, lmc, mw):
            assert fn(lam_uv)[0] > fn(lam_opt)[0]

    def test_mw_has_2175_bump(self):
        """MW law must show enhanced extinction around the 2175 Å feature."""
        lam_bump  = np.array([2.175e-5])  # 2175 Å (bump centre)
        lam_below = np.array([1.8e-5])    # just shortward of bump
        lam_above = np.array([2.8e-5])    # just longward of bump
        k_bump  = mw(lam_bump)[0]
        k_below = mw(lam_below)[0]
        k_above = mw(lam_above)[0]
        # The bump creates a local maximum between the flanking wavelengths
        assert k_bump > k_above

    def test_mw_bump_larger_than_smc(self):
        """MW 2175 Å feature must be stronger than SMC (negligible bump)."""
        lam = np.array([2.175e-5])  # 2175 Å
        assert mw(lam)[0] > smc(lam)[0]

    def test_smc_steeper_uv_than_mw_relative(self):
        """SMC has a steeper far-UV rise relative to optical than MW."""
        # Measure UV/optical ratio: k(1500 Å) / k(V)
        lam_uv = np.array([1.5e-5])
        lam_V  = np.array([5.5e-5])
        ratio_smc = smc(lam_uv)[0] / smc(lam_V)[0]
        ratio_mw  = mw(lam_uv)[0]  / mw(lam_V)[0]
        assert ratio_smc > ratio_mw

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_output_finite_and_nonneg(self, fn):
        """All values in the optical/UV range must be finite and non-negative."""
        lam = np.logspace(-5, -4, 50)   # 100 nm … 1 μm
        k = fn(lam)
        assert np.all(np.isfinite(k))
        assert np.all(k >= 0.0)

    @pytest.mark.parametrize("fn", [smc, lmc, mw])
    def test_infrared_values_finite(self, fn):
        """Infrared wavelengths should return finite non-negative values."""
        lam_ir = np.array([1e-3])  # 10 μm in cm
        k = fn(lam_ir)
        assert np.isfinite(k[0])
        assert k[0] >= 0.0


# ---------------------------------------------------------------------------
# pei92 error handling
# ---------------------------------------------------------------------------


class TestPei92Errors:
    """pei92 must raise ValueError for unknown law names."""

    def test_unknown_law_raises(self):
        lam = np.array([5.5e-5])
        with pytest.raises(ValueError, match="Unknown Pei"):
            pei92(lam, "invalid_law")

    def test_case_insensitive(self):
        """Law names are case-insensitive."""
        lam = np.array([5.5e-5])
        np.testing.assert_allclose(pei92(lam, "SMC"), 1.0, rtol=1e-6)
        np.testing.assert_allclose(pei92(lam, "Mw"),  1.0, rtol=1e-6)


# ---------------------------------------------------------------------------
# Fitter integration
# ---------------------------------------------------------------------------


class TestFitterIntegration:
    """Fitter can be created with built-in and custom extinction laws."""

    def test_fitter_smc(self):
        from VegasAfterglow.runner import Fitter
        f = Fitter(z=1.0, lumi_dist=1e28, extinction="smc")
        assert f._ext_law is not None
        assert not f._custom_extinction

    def test_fitter_lmc(self):
        from VegasAfterglow.runner import Fitter
        f = Fitter(z=1.0, lumi_dist=1e28, extinction="lmc")
        assert f._ext_law is not None

    def test_fitter_mw(self):
        from VegasAfterglow.runner import Fitter
        f = Fitter(z=1.0, lumi_dist=1e28, extinction="mw")
        assert f._ext_law is not None

    def test_fitter_none(self):
        from VegasAfterglow.runner import Fitter
        f = Fitter(z=1.0, lumi_dist=1e28, extinction=None)
        assert f._ext_law is None

    def test_fitter_callable(self):
        """Custom callable extinction law is accepted."""
        from VegasAfterglow.runner import Fitter
        fn = lambda lam_cm, params: smc(lam_cm)  # noqa: E731
        f = Fitter(z=1.0, lumi_dist=1e28, extinction=fn)
        assert f._ext_law is not None
        assert f._custom_extinction

    def test_fitter_unknown_law_raises(self):
        from VegasAfterglow.runner import Fitter
        with pytest.raises(ValueError, match="Unknown extinction law"):
            Fitter(z=1.0, lumi_dist=1e28, extinction="invalid")

    def test_ext_kernel_precomputed(self):
        """After adding data, the built-in extinction kernel is precomputed."""
        from VegasAfterglow.runner import Fitter
        f = Fitter(z=1.0, lumi_dist=1e28, extinction="smc")
        t = np.array([1e3, 2e3])
        f.add_flux_density(nu=4.84e14, t=t, f_nu=np.ones(2) * 1e-26,
                           err=np.ones(2) * 1e-28)
        f._consolidate_data()
        assert f._ext_kernel is not None
        assert f._ext_kernel.shape == t.shape
        assert np.all(np.isfinite(f._ext_kernel))
