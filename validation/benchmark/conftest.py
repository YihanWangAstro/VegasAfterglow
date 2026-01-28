"""
Pytest configuration for benchmark tests.

This module registers custom markers and provides shared fixtures
for the benchmark test suite.
"""

import pytest


def pytest_configure(config):
    """Register custom markers."""
    config.addinivalue_line(
        "markers", "benchmark: marks tests as benchmarks (may be slower)"
    )
    config.addinivalue_line(
        "markers", "convergence: marks tests as resolution convergence tests"
    )
    config.addinivalue_line(
        "markers", "full: marks tests as part of the full benchmark suite"
    )
    config.addinivalue_line(
        "markers", "slow: marks tests as slow (deselect with '-m \"not slow\"')"
    )


def pytest_collection_modifyitems(config, items):
    """
    Modify test collection based on markers.

    If running without explicit marker selection, skip slow tests by default.
    """
    # Check if user explicitly selected markers
    markexpr = config.getoption("-m", default="")

    if not markexpr:
        # No marker specified - skip slow tests by default
        skip_slow = pytest.mark.skip(reason="need -m slow or -m full to run")
        for item in items:
            if "slow" in item.keywords:
                item.add_marker(skip_slow)


@pytest.fixture(scope="session")
def benchmark_output_dir(tmp_path_factory):
    """Create a temporary directory for benchmark outputs."""
    return tmp_path_factory.mktemp("benchmark_output")


@pytest.fixture(scope="session")
def va_available():
    """Check if VegasAfterglow is available."""
    try:
        import VegasAfterglow
        return True
    except ImportError:
        return False


@pytest.fixture(autouse=True)
def skip_if_no_va(request, va_available):
    """Skip benchmark tests if VegasAfterglow is not available."""
    if "benchmark" in [m.name for m in request.node.iter_markers()]:
        if not va_available:
            pytest.skip("VegasAfterglow not available")
