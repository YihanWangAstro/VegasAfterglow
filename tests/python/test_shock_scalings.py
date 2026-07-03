"""Fast pytest bridge into the validation suite's shock-dynamics checks.

Imports the SAME scaling tables, model factory, and power-law fitter that
tests/validation/regression uses for the release report, and runs the Blandford-McKee
phase subset on every push. Single source of truth: if the expected exponents
ever change, they change in tests/validation/regression/run_regression.py and this
test follows automatically. The full three-phase sweep (coasting / BM / deep
Newtonian, plus reverse shock and spectral regimes) still runs via
`python tests/validation/run_validation.py`.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from validation.regression.run_regression import (  # noqa: E402
    SHOCK_SCALINGS,
    SLOPE_TOLERANCE,
    TIME_RANGES,
    VIZ_PARAMS,
    RegressionRunner,
    _fwd,
    create_model,
)
from validation.regression.utils import fit_powerlaw  # noqa: E402

pytestmark = pytest.mark.physics

PHASE = "BM"
QUANTITIES = ("u", "r", "B", "N_p")


@pytest.fixture(scope="module")
def models():
    return {m: create_model(m, **VIZ_PARAMS.get(m, {})) for m in ("ISM", "wind")}


@pytest.mark.parametrize("medium", ["ISM", "wind"])
@pytest.mark.parametrize("qty", QUANTITIES)
def test_bm_phase_scaling(models, medium, qty):
    """For each medium (ISM/wind) and forward-shock quantity (four-velocity u, radius r, comoving B, swept-up proton number N_p), the power-law slope fitted over the Blandford-McKee phase time range matches the analytic self-similar exponent from the shared SHOCK_SCALINGS table within SLOPE_TOLERANCE (0.1 in the slope)."""
    expected = float(SHOCK_SCALINGS[PHASE][medium][qty])
    t_range = TIME_RANGES[PHASE][medium]
    details = models[medium].details(t_range[0], t_range[1])
    if qty == "u":
        t, Gamma = _fwd(details, "Gamma")
        y = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
    else:
        t, y = _fwd(details, RegressionRunner._SHOCK_ATTR_MAP[qty])
    mask = (t >= t_range[0]) & (t <= t_range[1]) & (y > 0) & np.isfinite(y)
    assert np.sum(mask) >= 5, f"too few valid points for {medium}/{qty}"
    measured = fit_powerlaw(t[mask], y[mask])
    assert abs(measured - expected) < SLOPE_TOLERANCE, (
        f"{medium} BM-phase {qty}: measured slope {measured:.3f}, "
        f"expected {expected:.3f} (tol {SLOPE_TOLERANCE})"
    )
