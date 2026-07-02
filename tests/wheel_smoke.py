"""Wheel smoke test run by cibuildwheel after each wheel build.

Exercises the full compute path -- imports, C++ extension load, jet/medium/observer
construction, and a real flux_density call -- so a wheel with a broken extension or
missing symbol fails before publish, not after.
"""

import numpy as np

from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet


def main() -> None:
    model = Model(
        TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300),
        ISM(n_ism=1.0),
        Observer(lumi_dist=1e28, z=1.0, theta_obs=0.0),
        Radiation(eps_e=0.1, eps_B=0.01, p=2.2),
        resolutions=(0.3, 1.0, 10),
    )
    t = np.logspace(3, 6, 10)
    nu = np.full_like(t, 1e14)
    flux = model.flux_density(t, nu).total

    assert flux.shape == t.shape, f"shape mismatch: {flux.shape} vs {t.shape}"
    assert np.all(np.isfinite(flux)), f"non-finite flux values: {flux}"
    assert np.all(flux > 0), f"non-positive flux values: {flux}"
    print("OK")


if __name__ == "__main__":
    main()
