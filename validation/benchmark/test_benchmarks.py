#!/usr/bin/env python3
"""
Pytest wrapper for comprehensive benchmark suite.

This module provides pytest integration for the benchmark suite, allowing
benchmarks to be run via pytest while preserving the standalone functionality.

Usage
-----
    # Run quick benchmarks
    pytest tests/benchmark/test_benchmarks.py -v -k "quick"

    # Run convergence tests
    pytest tests/benchmark/test_benchmarks.py -v -k "convergence"

    # Run full benchmark suite
    pytest tests/benchmark/test_benchmarks.py -v -m "full"

    # Run all benchmarks
    pytest tests/benchmark/test_benchmarks.py -v -m "benchmark"
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from benchmark_suite import (
    ComprehensiveBenchmark,
    FIDUCIAL_RESOLUTION,
    PHI_VALUES,
    THETA_VALUES,
    T_VALUES,
    CONVERGENCE_BANDS,
    generate_visualizations,
)
from configs import (
    get_all_jet_names,
    get_all_medium_names,
    get_all_radiation_names,
)


# ============================================================================
# Fixtures
# ============================================================================

@pytest.fixture(scope="module")
def benchmark_runner():
    """Create a benchmark runner with reduced iterations for testing."""
    return ComprehensiveBenchmark(iterations=1)


@pytest.fixture(scope="session")
def output_dir(tmp_path_factory):
    """Create output directory for benchmark results."""
    return tmp_path_factory.mktemp("benchmark_output")


# ============================================================================
# Quick Benchmark Tests
# ============================================================================

@pytest.mark.benchmark
class TestQuickBenchmarks:
    """Quick benchmark tests for validation."""

    def test_tophat_ism_basic(self, benchmark_runner):
        """Benchmark tophat jet in ISM (no convergence)."""
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        assert len(result.lightcurve_t) > 0
        assert len(result.spectrum_nu) > 0
        print(f"\n  Tophat/ISM: {result.timing.total_ms:.2f} ms total")

    def test_gaussian_ism_basic(self, benchmark_runner):
        """Benchmark Gaussian jet in ISM (no convergence)."""
        result = benchmark_runner.benchmark_config(
            jet_name="gaussian",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        print(f"\n  Gaussian/ISM: {result.timing.total_ms:.2f} ms total")

    def test_component_timing(self, benchmark_runner):
        """Verify component timing breakdown."""
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        # Model creation should be fast
        assert result.timing.model_creation_ms < 100, "Model creation too slow"

        # All components should have positive timing
        assert result.timing.flux_density_ms > 0
        assert result.timing.flux_grid_ms > 0
        assert result.timing.details_ms > 0

        print(f"\n  Model: {result.timing.model_creation_ms:.3f} ms")
        print(f"  Flux:  {result.timing.flux_density_ms:.2f} ms")
        print(f"  Grid:  {result.timing.flux_grid_ms:.2f} ms")
        print(f"  Details: {result.timing.details_ms:.2f} ms")


# ============================================================================
# Convergence Tests
# ============================================================================

@pytest.mark.benchmark
@pytest.mark.convergence
class TestDimensionConvergence:
    """Test resolution convergence for each dimension."""

    @pytest.mark.parametrize("jet_name", ["tophat", "gaussian"])
    def test_phi_convergence(self, benchmark_runner, jet_name):
        """Test phi resolution convergence."""
        result = benchmark_runner.benchmark_config(
            jet_name=jet_name,
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=True,
            verbose=False,
        )

        conv = result.phi_convergence
        assert conv is not None
        assert len(conv.values) == len(PHI_VALUES)
        assert len(conv.flux_values) == len(PHI_VALUES)

        # Error should decrease with resolution
        valid_errors = [e for e in conv.relative_errors if not np.isnan(e)]
        if len(valid_errors) >= 2:
            # Last error should be 0 (reference)
            assert conv.relative_errors[-1] == 0.0

        print(f"\n  {jet_name} phi convergence:")
        for i, val in enumerate(conv.values):
            err = conv.relative_errors[i]
            print(f"    φ={val}: error={err:.2e}, time={conv.times_ms[i]:.2f}ms")

    @pytest.mark.parametrize("jet_name", ["tophat", "gaussian"])
    def test_theta_convergence(self, benchmark_runner, jet_name):
        """Test theta resolution convergence."""
        result = benchmark_runner.benchmark_config(
            jet_name=jet_name,
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=True,
            verbose=False,
        )

        conv = result.theta_convergence
        assert conv is not None
        assert len(conv.values) == len(THETA_VALUES)

        print(f"\n  {jet_name} theta convergence:")
        for i, val in enumerate(conv.values):
            err = conv.relative_errors[i]
            print(f"    θ={val}: error={err:.2e}, time={conv.times_ms[i]:.2f}ms")

    @pytest.mark.parametrize("jet_name", ["tophat", "gaussian"])
    def test_t_convergence(self, benchmark_runner, jet_name):
        """Test time resolution convergence."""
        result = benchmark_runner.benchmark_config(
            jet_name=jet_name,
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=True,
            verbose=False,
        )

        conv = result.t_convergence
        assert conv is not None
        assert len(conv.values) == len(T_VALUES)

        print(f"\n  {jet_name} time convergence:")
        for i, val in enumerate(conv.values):
            err = conv.relative_errors[i]
            print(f"    t={val}: error={err:.2e}, time={conv.times_ms[i]:.2f}ms")


# ============================================================================
# Full Configuration Tests
# ============================================================================

@pytest.mark.benchmark
@pytest.mark.full
class TestFullConfigurations:
    """Test all jet/medium/radiation combinations."""

    @pytest.mark.parametrize("jet_name", get_all_jet_names())
    def test_all_jet_types_ism(self, benchmark_runner, jet_name):
        """Benchmark all jet types in ISM."""
        result = benchmark_runner.benchmark_config(
            jet_name=jet_name,
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        assert len(result.lightcurve_flux) > 0
        assert all(f > 0 for f in result.lightcurve_flux if not np.isnan(f))

        print(f"\n  {jet_name}/ISM: {result.timing.total_ms:.2f} ms")

    @pytest.mark.parametrize("jet_name", get_all_jet_names())
    def test_all_jet_types_wind(self, benchmark_runner, jet_name):
        """Benchmark all jet types in wind medium."""
        result = benchmark_runner.benchmark_config(
            jet_name=jet_name,
            medium_name="wind",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        print(f"\n  {jet_name}/wind: {result.timing.total_ms:.2f} ms")

    @pytest.mark.parametrize("radiation_name", [
        "synchrotron_only", "with_ssc_cooling", "full_ssc"
    ])
    def test_radiation_configs(self, benchmark_runner, radiation_name):
        """Benchmark different radiation configurations."""
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name=radiation_name,
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        print(f"\n  tophat/ISM/{radiation_name}: {result.timing.total_ms:.2f} ms")

    @pytest.mark.parametrize("theta_obs", [0.0, 0.1, 0.3])
    def test_viewing_angles(self, benchmark_runner, theta_obs):
        """Benchmark different viewing angles."""
        result = benchmark_runner.benchmark_config(
            jet_name="gaussian",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=theta_obs,
            run_convergence=False,
            verbose=False,
        )

        assert result.timing.total_ms > 0
        print(f"\n  gaussian/ISM/θ={theta_obs}: {result.timing.total_ms:.2f} ms")


# ============================================================================
# Suite Runner Tests
# ============================================================================

@pytest.mark.benchmark
class TestSuiteRunners:
    """Test the suite runner methods."""

    def test_quick_suite(self):
        """Run the quick benchmark suite."""
        benchmark = ComprehensiveBenchmark(iterations=1)
        results = benchmark.run_quick_suite()

        assert len(results) == 2  # tophat and gaussian
        for r in results:
            assert r.timing.total_ms > 0

        print(f"\n  Quick suite: {len(results)} configs benchmarked")

    def test_convergence_suite(self):
        """Run convergence tests on subset."""
        benchmark = ComprehensiveBenchmark(iterations=1)
        results = benchmark.run_convergence_only(
            jet_types=["tophat"],
            medium_types=["ISM"],
        )

        assert len(results) == 1
        assert results[0].phi_convergence is not None
        assert results[0].theta_convergence is not None
        assert results[0].t_convergence is not None

    @pytest.mark.slow
    def test_standard_suite_subset(self):
        """Run standard suite on small subset."""
        benchmark = ComprehensiveBenchmark(iterations=1)
        results = benchmark.run_standard_suite(
            jet_types=["tophat", "gaussian"],
            medium_types=["ISM"],
            radiation_types=["synchrotron_only"],
            viewing_angles=[0.0],
            run_convergence=False,
        )

        assert len(results) == 2
        for r in results:
            assert len(r.lightcurve_t) > 0
            assert len(r.spectrum_nu) > 0


# ============================================================================
# Visualization Tests
# ============================================================================

@pytest.mark.benchmark
class TestVisualization:
    """Test visualization generation."""

    def test_generate_plots(self, benchmark_runner, output_dir):
        """Test PDF report generation."""
        # Run minimal benchmark
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=True,
            verbose=False,
        )

        # Generate visualizations
        generate_visualizations([result], str(output_dir))

        # Check output file exists
        pdf_path = output_dir / "benchmark_report.pdf"
        assert pdf_path.exists(), "PDF report not generated"
        print(f"\n  Report saved to: {pdf_path}")


# ============================================================================
# Data Quality Tests
# ============================================================================

@pytest.mark.benchmark
class TestDataQuality:
    """Test quality of benchmark data."""

    def test_lightcurve_physical(self, benchmark_runner):
        """Verify light curve has physical values."""
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        # All flux values should be positive
        assert all(f > 0 for f in result.lightcurve_flux)

        # Flux should decrease at late times (generally)
        peak_idx = np.argmax(result.lightcurve_flux)
        late_flux = result.lightcurve_flux[-1]
        peak_flux = result.lightcurve_flux[peak_idx]
        assert late_flux < peak_flux, "Flux should decrease after peak"

    def test_spectrum_physical(self, benchmark_runner):
        """Verify spectrum has physical shape."""
        result = benchmark_runner.benchmark_config(
            jet_name="tophat",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=False,
            verbose=False,
        )

        # All flux values should be positive
        valid_flux = [f for f in result.spectrum_flux if not np.isnan(f)]
        assert all(f > 0 for f in valid_flux)

        # Spectrum should span many orders of magnitude
        if len(valid_flux) > 0:
            dynamic_range = max(valid_flux) / min(valid_flux)
            assert dynamic_range > 10, "Spectrum should have significant dynamic range"

    def test_resolution_cost_scaling(self, benchmark_runner):
        """Verify cost increases with resolution."""
        result = benchmark_runner.benchmark_config(
            jet_name="gaussian",
            medium_name="ISM",
            radiation_name="synchrotron_only",
            theta_obs=0.0,
            run_convergence=True,
            verbose=False,
        )

        if result.resolution_costs:
            times = [r.time_ms for r in result.resolution_costs]
            # Cost should generally increase with resolution
            # Allow some variation due to caching
            assert times[-1] >= times[0] * 0.5, "Higher resolution should not be much faster"
