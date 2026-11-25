"""Every cell must contribute at the earliest requested time.

A cell adds no flux before its first time-grid node arrives. A theta row solved
once and reused for every azimuth therefore has to start early enough for the
azimuth farthest from the line of sight. When it started for the nearest one,
off-axis light curves lost 7-40% at their earliest time, and the loss
disappeared as soon as an earlier time shared the call.
"""
import numpy as np
import pytest

from VegasAfterglow import ISM, Medium, Model, Observer, Radiation, TophatJet

pytestmark = pytest.mark.physics

M_P = 1.67262192369e-24
# Late enough that every row starts above the solver's 0.1 s floor in both calls,
# so the dynamics are identical and only the coverage differs.
T = np.logspace(4, 6, 6)
T_EARLIER = np.r_[0.5 * T[0], T]
NU = 5e14
RAD = Radiation(eps_e=0.1, eps_B=1e-3, p=2.3)
# Dense time lattice: moving the lattice start then changes the flux by < 1e-3.
RES = (0.1, 0.3, 20)


def _uniform_medium():
    """ISM(n_ism=1) written as a custom profile, which is solved row by row."""
    return Medium(rho=lambda phi, theta, r: M_P, isotropic=False)


CASES = {
    "non_isotropic_medium": dict(medium=_uniform_medium, theta_obs=0.03),
    "spreading_jet": dict(medium=lambda: ISM(n_ism=1.0), theta_obs=0.03, spreading=True),
    "far_off_axis": dict(medium=lambda: ISM(n_ism=1.0), theta_obs=0.5),
    "reverse_shock": dict(medium=_uniform_medium, theta_obs=0.03, rvs=True),
}


def _components(t, medium, theta_obs, spreading=False, rvs=False):
    """Synchrotron flux of each shock at the last ``len(T)`` times of ``t``."""
    model = Model(jet=TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300, spreading=spreading), medium=medium(),
                  observer=Observer(lumi_dist=1e28, z=0.5, theta_obs=theta_obs), fwd_rad=RAD,
                  rvs_rad=RAD if rvs else None, resolutions=RES)
    flux = model.flux_density(t, np.full_like(t, NU))
    out = {"fwd": np.asarray(flux.fwd.sync)[-len(T):]}
    if rvs:
        out["rvs"] = np.asarray(flux.rvs.sync)[-len(T):]
    return out


@pytest.mark.parametrize("case", CASES)
def test_flux_independent_of_earlier_times_in_the_call(case):
    """Each shock's flux at 1e4-1e6 s is the same to 3e-3 whether or not an earlier time shares the call."""
    alone = _components(T, **CASES[case])
    with_earlier = _components(T_EARLIER, **CASES[case])
    for name in alone:
        np.testing.assert_allclose(alone[name], with_earlier[name], rtol=3e-3, err_msg=f"{case}: {name}")
