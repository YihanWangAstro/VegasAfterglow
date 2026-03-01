"""
Pytest configuration and fixtures for Python API tests.
"""

import pytest
import numpy as np


@pytest.fixture
def sample_time_array():
    """Generate a sample time array in seconds."""
    return np.logspace(2, 6, 30)


@pytest.fixture
def sample_frequency():
    """Return a sample optical frequency in Hz."""
    return 4.84e14  # R-band


@pytest.fixture
def sample_frequency_array():
    """Generate a sample frequency array in Hz."""
    return np.logspace(9, 18, 20)
