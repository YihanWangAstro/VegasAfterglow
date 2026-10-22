"""The blast-wave solution must not depend on where the time lattice starts.

The forward-shock ODE starts at t0 = min(first lattice node, floor, 0.1 t_dec).
Requesting earlier observer times moves t0 earlier; the light curve over the
common window must stay the same. In a wind the mass swept before t0 is a
finite fraction of the deceleration mass, so an init that ignores it shifts
the whole light curve (was ~9% for the wind case below).
"""
import numpy as np
import pytest

from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet, Wind

pytestmark = pytest.mark.physics

T = np.logspace(2, 7, 30)
T_EARLY = np.concatenate([np.logspace(-5, 1.9, 20), T])
NU = np.array([1e9, 5e14, 1e18])
RAD = Radiation(eps_e=0.1, eps_B=1e-3, p=2.3)
# Dense time lattice (lattice-phase interpolation differences < 1e-4) and tight ODE
# tolerance, so only the start time differs.
RES = (0.06, 0.15, 40)


def _flux(medium, t, radiative=True):
    model = Model(jet=TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300), medium=medium,
                  observer=Observer(lumi_dist=1e27, z=0.1, theta_obs=0.0), fwd_rad=RAD,
                  rtol=1e-8, resolutions=RES, radiative_fireball=radiative)
    return np.asarray(model.flux_density_grid(t, NU).total)[:, -len(T):]


@pytest.mark.parametrize("radiative", [True, False])
@pytest.mark.parametrize("medium", [Wind(A_star=0.1), Wind(A_star=1.0), ISM(n_ism=1.0)],
                         ids=["wind_A0.1", "wind_A1", "ism"])
def test_light_curve_independent_of_lattice_start(medium, radiative):
    """Forward-shock light curves over 1e2-1e7 s agree to 2e-3 whether or not the lattice extends down to 1e-5 s."""
    late = _flux(medium, T, radiative)
    early = _flux(medium, T_EARLY, radiative)
    np.testing.assert_allclose(late, early, rtol=2e-3)
