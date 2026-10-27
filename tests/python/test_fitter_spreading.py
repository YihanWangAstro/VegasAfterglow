"""Lateral jet spreading in fits: ``Fitter(spreading=True)`` (experimental).

The flag reaches the built-in jet constructors, is rejected alongside a custom
jet factory (which owns its own flag), and survives save/load.
"""
import types

import numpy as np
import pytest

from VegasAfterglow import (
    ISM,
    FitResult,
    Fitter,
    Model,
    ModelParams,
    Observer,
    ParamDef,
    Radiation,
    Scale,
    TophatJet,
)
from VegasAfterglow.fitting.config import JETS
from VegasAfterglow.fitting.utils import _build_transformer

pytestmark = pytest.mark.physics

T = np.logspace(4.0, 7.0, 12)
NU = np.full_like(T, 1e15)


@pytest.fixture(scope="module")
def spreading_flux():
    truth = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300, spreading=True),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=2.3),
    )
    return np.asarray(truth.flux_density(T, NU).total)


def _params():
    p = ModelParams()
    p.E_iso, p.Gamma0, p.theta_c, p.theta_v = 1e53, 300, 0.1, 0.0
    p.p, p.eps_e, p.eps_B, p.n_ism = 2.3, 0.1, 1e-3, 1.0
    return p


def _fitter(**kw):
    return Fitter(z=0.5, lumi_dist=3e28, jet="tophat", medium="ism", **kw)


def test_spreading_defaults_off():
    """The Fitter runs without spreading unless asked."""
    assert _fitter().spreading is False


def test_spreading_reaches_the_model(spreading_flux):
    """Data generated with a spreading jet is matched by Fitter(spreading=True) and missed without it."""
    err = 0.05 * spreading_flux
    on, off = _fitter(spreading=True), _fitter()
    for f in (on, off):
        f.add_flux_density(nu=1e15, t=T, f_nu=spreading_flux, err=err)
    chi2_on, chi2_off = on._evaluate(_params()), off._evaluate(_params())
    assert chi2_on < 1e-6
    assert chi2_off > 1.0


@pytest.mark.parametrize("jet_type", sorted(JETS))
def test_every_builtin_jet_accepts_spreading(jet_type):
    """The registry's jet constructors all take the spreading flag the default factory passes."""
    p = types.SimpleNamespace(**vars(ModelParams()))
    p.theta_c, p.E_iso, p.Gamma0 = 0.1, 1e53, 300.0
    p.theta_w, p.E_iso_w, p.Gamma0_w, p.k_e, p.k_g = 0.3, 1e51, 50.0, 2.0, 2.0
    fitter = Fitter(z=0.5, lumi_dist=3e28, jet=jet_type, medium="ism", spreading=True)
    assert fitter.jet_factory(p) is not None


def test_spreading_with_custom_jet_factory_is_rejected():
    """A custom factory owns its spreading flag, so combining it with the Fitter option raises."""
    factory = lambda p: TophatJet(p.theta_c, p.E_iso, p.Gamma0)  # noqa: E731
    with pytest.raises(ValueError, match="spreading"):
        Fitter(z=0.5, lumi_dist=3e28, jet=factory, medium="ism", spreading=True)


def test_spreading_shows_in_repr():
    """repr lists the flag only when it differs from the default."""
    assert "spreading" in repr(_fitter(spreading=True))
    assert "spreading" not in repr(_fitter())


def test_spreading_survives_save_and_load(tmp_path, spreading_flux):
    """The flag round-trips through save/load and the loaded fitter evaluates identically."""
    fitter = _fitter(spreading=True)
    fitter.add_flux_density(nu=1e15, t=T, f_nu=spreading_flux, err=0.05 * spreading_flux)
    defs = [ParamDef("E_iso", 1e50, 1e54, Scale.log)]
    fitter._param_defs = defs
    fitter._to_params = _build_transformer(defs)
    fitter.result = FitResult(
        samples=np.zeros((4, 1, 1)), log_probs=np.zeros((4, 1)), labels=("log10_E_iso",),
        top_k_params=np.array([[53.0]]), top_k_log_probs=np.array([0.0]),
    )
    chi2 = fitter._evaluate(_params())
    path = tmp_path / "fit.h5"
    fitter.save(path)
    loaded = Fitter.load(path)
    assert loaded.spreading is True
    assert loaded._evaluate(_params()) == pytest.approx(chi2, rel=1e-12)
