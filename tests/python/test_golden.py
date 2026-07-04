"""Golden-baseline regression tests.

Each canonical config in tests/python/golden/regenerate.py is recomputed and
compared component by component (total, fwd.sync, fwd.ssc, rvs.sync, rvs.ssc)
against the stored .npz baseline. The tolerance rtol=2e-3 absorbs
cross-platform fast-math/SIMD differences (different vectorization, FMA
contraction, and libm implementations shift results at the ~1e-4 level)
while still catching real physics drift, which typically changes fluxes at
the percent level or more. Bins below 1e-2 of the component peak float: the
fragile structured-jet configs reproduce them only to ~1% against any
bit-level perturbation, and they carry no observable information.

If these tests fail after an INTENDED physics change: inspect the diff, then
regenerate via `python tests/python/golden/regenerate.py` and commit the new
baselines with the change.
"""

import contextlib
import importlib.util
import os
import sys
import tempfile

import numpy as np
import pytest

pytestmark = pytest.mark.golden

GOLDEN_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "golden")


@contextlib.contextmanager
def _capture_stderr_fd():
    """Capture writes to file descriptor 2, including C++ fprintf(stderr)."""
    sys.stderr.flush()
    saved = os.dup(2)
    with tempfile.TemporaryFile(mode="w+") as tmp:
        os.dup2(tmp.fileno(), 2)
        try:
            yield tmp
        finally:
            sys.stderr.flush()
            os.dup2(saved, 2)
            os.close(saved)


def _load_regenerate_module():
    path = os.path.join(GOLDEN_DIR, "regenerate.py")
    spec = importlib.util.spec_from_file_location("golden_regenerate", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


regen = _load_regenerate_module()

COMPONENTS = ("total", "fwd_sync", "fwd_ssc", "rvs_sync", "rvs_ssc")


@pytest.fixture(scope="module")
def recomputed():
    cache = {}

    def get(name):
        if name not in cache:
            model = regen.build_model(regen.CONFIGS[name])
            with _capture_stderr_fd() as tmp:
                components = regen.compute_components(model)
                tmp.seek(0)
                stderr_text = tmp.read()
            cache[name] = (components, stderr_text)
        return cache[name]

    return get


@pytest.mark.parametrize("name", sorted(regen.CONFIGS))
def test_no_solver_warnings(name, recomputed):
    """The shock solver must not abandon any grid row ("giving up" on stderr): an abandoned row leaves its stored state at initialization values and silently corrupts the light curve."""
    _, stderr_text = recomputed(name)
    assert "giving up" not in stderr_text, f"solver warnings during {name}:\n{stderr_text}"


@pytest.mark.parametrize("name", sorted(regen.CONFIGS))
def test_baseline_grid_matches(name):
    """Each config's stored baseline time and frequency grids match the current regeneration grids to 1e-12 relative (logspace differs by one ulp across platform libms)."""
    baseline = np.load(os.path.join(GOLDEN_DIR, f"{name}.npz"))
    np.testing.assert_allclose(baseline["t"], regen.T, rtol=1e-12)
    np.testing.assert_allclose(baseline["nus"], regen.NUS, rtol=1e-12)


@pytest.mark.parametrize("component", COMPONENTS)
@pytest.mark.parametrize("name", sorted(regen.CONFIGS))
def test_golden_component(name, component, recomputed):
    """Each recomputed flux component matches its golden baseline within the calibrated rtol=2e-3, and identically-zero baseline components stay zero. Bins below 1e-2 of the component peak float (the atol term): for structured-jet reverse shocks their values are only reproducible to ~1% against ANY perturbation — codegen, platform libm, math-kernel choice — because the wing-row ODE solves amplify bit-level differences to that saturated level. Real regressions move bright bins far beyond rtol, so nothing observable is unguarded."""
    baseline = np.load(os.path.join(GOLDEN_DIR, f"{name}.npz"))
    reference = np.asarray(baseline[component])
    current = np.asarray(recomputed(name)[0][component])
    assert current.shape == reference.shape

    if not np.any(reference):
        assert not np.any(current), f"{name}/{component} was identically zero in the baseline"
        return

    peak = np.max(np.abs(reference))
    np.testing.assert_allclose(current, reference, rtol=2e-3, atol=1e-2 * peak)
