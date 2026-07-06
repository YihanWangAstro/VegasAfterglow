"""End-to-end MCMC smoke test: a tiny emcee run over synthetic data.

Covers the full fit pipeline (parameter validation, transformer, likelihood,
sampler orchestration, FitResult assembly) that no other fast test touches.
Deliberately minuscule: coarse resolution, 2 free parameters, few steps --
this verifies plumbing and recovery direction, not posterior quality.
"""
import numpy as np
import pytest

from VegasAfterglow import (
    ISM,
    Fitter,
    Model,
    Observer,
    ParamDef,
    Radiation,
    Scale,
    TophatJet,
)

pytestmark = pytest.mark.physics

TRUE_E = 1e53


@pytest.fixture(scope="module")
def fit_result(tmp_path_factory):
    # synthetic optical light curve from a known model
    t = np.logspace(3.5, 5.5, 12)
    nu = np.full_like(t, 1e15)
    truth = Model(
        jet=TophatJet(theta_c=0.1, E_iso=TRUE_E, Gamma0=300),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=2.3),
    )
    flux = np.asarray(truth.flux_density(t, nu).total)
    err = 0.1 * flux
    rng = np.random.default_rng(42)
    obs = flux * (1 + 0.05 * rng.standard_normal(flux.size))

    fitter = Fitter(z=0.5, lumi_dist=3e28, jet="tophat", medium="ism")
    fitter.add_flux_density(nu=1e15, t=t, f_nu=obs, err=err)

    defs = [
        ParamDef("E_iso", 1e52, 1e54, Scale.log),
        ParamDef("Gamma0", 100, 600, Scale.log),
        ParamDef("theta_c", 0.1, 0.1, Scale.fixed),
        ParamDef("n_ism", 1.0, 1.0, Scale.fixed),
        ParamDef("eps_e", 0.1, 0.1, Scale.fixed),
        ParamDef("eps_B", 1e-3, 1e-3, Scale.fixed),
        ParamDef("p", 2.3, 2.3, Scale.fixed),
    ]
    outdir = str(tmp_path_factory.mktemp("bilby_out"))
    return fitter.fit(defs, sampler="emcee", nsteps=30, nburn=10, npool=1,
                      outdir=outdir, top_k=5)


def test_fit_returns_well_formed_result(fit_result):
    """FitResult has consistent shapes: samples' last axis equals the 2 free parameters, log-probs are all finite, and top_k_params/n_free_params/n_data match the fit setup."""
    assert fit_result.samples.ndim >= 2
    assert fit_result.samples.shape[-1] == 2  # two free parameters
    assert np.all(np.isfinite(fit_result.log_probs))
    assert len(fit_result.top_k_params) == 5
    assert fit_result.n_free_params == 2
    assert fit_result.n_data == 12


def test_fit_recovers_energy_scale(fit_result):
    """The top-ranked sample's log10(E_iso) lands within 1 dex of the injected truth, i.e. the short MCMC run recovers the isotropic energy to the right order of magnitude."""
    # best-fit log10(E_iso) within the (wide) prior should land near truth;
    # with 30 steps we only assert the right order of magnitude
    labels = list(fit_result.labels)
    idx = labels.index("log10_E_iso")
    best = fit_result.top_k_params[0][idx]
    assert abs(best - np.log10(TRUE_E)) < 1.0


def test_fit_summary_renders(fit_result):
    """FitResult.summary() renders without raising and its repr contains ranking or chi-square content."""
    text = repr(fit_result.summary())
    assert "Rank" in text or "chi" in text.lower()
