"""Upper limits in the fitting likelihood.

A non-detection quoted as an n-sigma limit U has noise sigma = U / n_sigma;
its likelihood is the probability that a measurement would have fallen below
U given the model, Phi((U - F_model) / sigma), evaluated in linear flux.
"""
import numpy as np
import pytest
from scipy.special import log_ndtr

from VegasAfterglow import (
    ISM,
    Fitter,
    Model,
    ModelParams,
    Observer,
    Radiation,
    TophatJet,
)

pytestmark = pytest.mark.physics

T = np.logspace(3.5, 5.5, 10)
NU = np.full_like(T, 1e15)
BAND = (2e17, 2e18)


@pytest.fixture(scope="module")
def truth_flux():
    truth = Model(
        jet=TophatJet(theta_c=0.1, E_iso=1e53, Gamma0=300),
        medium=ISM(n_ism=1.0),
        observer=Observer(lumi_dist=3e28, z=0.5, theta_obs=0.0),
        fwd_rad=Radiation(eps_e=0.1, eps_B=1e-3, p=2.3),
    )
    return np.asarray(truth.flux_density(T, NU).total), np.asarray(
        truth.flux(T, BAND[0], BAND[1], 5).total
    )


def _params():
    p = ModelParams()
    p.E_iso, p.Gamma0, p.theta_c, p.theta_v = 1e53, 300, 0.1, 0.0
    p.p, p.eps_e, p.eps_B, p.n_ism = 2.3, 0.1, 1e-3, 1.0
    return p


def _fitter():
    return Fitter(z=0.5, lumi_dist=3e28, jet="tophat", medium="ism")


def _model_flux(fitter, t, nu):
    """The fitter's own model flux at exactly these points (same time lattice as the evaluation)."""
    return np.asarray(fitter._build_model(_params()).flux_density(t, np.full_like(t, nu)).total)


def test_limit_far_above_model_costs_almost_nothing(truth_flux):
    """A 3-sigma limit ten times above the model costs less than 0.01 per limit (the floor is -2 ln Phi(3) ~ 0.003)."""
    flux, _ = truth_flux
    fitter = _fitter()
    fitter.add_flux_density_limits(nu=1e15, t=T, limit=10 * flux)
    assert fitter._evaluate(_params()) < 0.01 * len(T)


def test_limit_equal_to_model_costs_two_ln_two(truth_flux):
    """A limit sitting exactly on the model costs -2 ln(1/2) per limit."""
    fitter = _fitter()
    t = T[:3]
    fitter.add_flux_density_limits(nu=1e15, t=t, limit=_model_flux(fitter, t, 1e15))
    assert fitter._evaluate(_params()) == pytest.approx(3 * 2 * np.log(2), rel=1e-9)


def test_limit_penalty_follows_the_gaussian_tail(truth_flux):
    """A model k sigma above the limit costs -2 ln Phi(-k): with n_sigma=3 and U = F/2, k = 3."""
    fitter = _fitter()
    t = T[:1]
    fitter.add_flux_density_limits(nu=1e15, t=t, limit=_model_flux(fitter, t, 1e15) / 2, n_sigma=3)
    assert fitter._evaluate(_params()) == pytest.approx(-2 * log_ndtr(-3.0), rel=1e-9)


def test_explicit_sigma_overrides_n_sigma(truth_flux):
    """When sigma is given it sets the noise behind the limit and n_sigma is ignored."""
    fitter = _fitter()
    t = T[:1]
    U = _model_flux(fitter, t, 1e15) / 2
    fitter.add_flux_density_limits(nu=1e15, t=t, limit=U, sigma=U / 5, n_sigma=3)
    assert fitter._evaluate(_params()) == pytest.approx(-2 * log_ndtr(-5.0), rel=1e-9)


def test_limits_add_to_detection_chi2(truth_flux):
    """chi2 is the detection sum plus -2 ln Phi((U - F)/sigma) over the limits, from one model evaluation."""
    flux, _ = truth_flux
    obs = 1.3 * flux
    err = 0.1 * obs
    t_lim = T[:2]
    U = flux[:2]
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err)
    fitter.add_flux_density_limits(nu=5e14, t=t_lim, limit=U)
    chi2 = fitter._evaluate(_params())

    model = fitter._build_model(_params())
    eval_t = np.concatenate([T, t_lim])
    eval_nu = np.concatenate([NU, np.full_like(t_lim, 5e14)])
    order = np.argsort(eval_t, kind="stable")
    model_flux = np.asarray(model.flux_density(eval_t[order], eval_nu[order]).total)
    back = np.empty_like(order)
    back[order] = np.arange(len(order))
    f_det = model_flux[back[: len(T)]]
    f_lim = model_flux[back[len(T):]]
    expected = np.sum(((np.log(obs) - np.log(f_det)) / (err / obs)) ** 2) - 2 * np.sum(
        log_ndtr((U - f_lim) / (U / 3))
    )
    assert chi2 == pytest.approx(expected, rel=1e-10)


def test_spectrum_limits(truth_flux):
    """add_spectrum_limits takes one time and many frequencies."""
    fitter = _fitter()
    F = _model_flux(fitter, T[:1], 1e15)[0]
    fitter.add_spectrum_limits(t=T[0], nu=[1e15, 1e15], limit=[F, F])
    assert fitter._evaluate(_params()) == pytest.approx(2 * 2 * np.log(2), rel=1e-9)


def test_band_integrated_limits(truth_flux):
    """add_flux_limits applies the same likelihood to band-integrated limits."""
    fitter = _fitter()
    t = T[:2]
    F = np.asarray(fitter._build_model(_params()).flux(t, BAND[0], BAND[1], 5).total)
    fitter.add_flux_limits(band=BAND, t=t, limit=F / 2, n_sigma=3)
    assert fitter._evaluate(_params()) == pytest.approx(-2 * 2 * log_ndtr(-3.0), rel=1e-9)


def test_limits_use_the_extincted_model_flux(truth_flux):
    """Host extinction attenuates the model before it is compared with an optical limit."""
    flux, _ = truth_flux
    fitter = Fitter(z=0.5, lumi_dist=3e28, jet="tophat", medium="ism", extinction="smc")
    fitter.add_flux_density_limits(nu=1e15, t=T[:1], limit=flux[:1])
    p = _params()
    p.A_V = 1.0
    # the attenuated model sits below the limit, so the cost drops below 2 ln 2
    assert fitter._evaluate(p) < 2 * np.log(2) - 0.1


@pytest.mark.parametrize(
    "kwargs, match",
    [
        (dict(limit=[-1e-28]), "limit"),
        (dict(limit=[1e-28], sigma=[0.0]), "sigma"),
        (dict(limit=[1e-28], n_sigma=0), "n_sigma"),
        (dict(limit=[1e-28, 1e-28]), "shape"),
    ],
)
def test_limit_validation(kwargs, match):
    fitter = _fitter()
    with pytest.raises(ValueError, match=match):
        fitter.add_flux_density_limits(nu=1e15, t=[1e4], **kwargs)


def test_limits_survive_save_and_load(tmp_path, truth_flux):
    """Point and band limits round-trip through save/load and the loaded fitter evaluates identically."""
    from VegasAfterglow import FitResult, ParamDef, Scale
    from VegasAfterglow.fitting.utils import _build_transformer

    flux, band_flux = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux)
    fitter.add_flux_density_limits(nu=5e14, t=T[:2], limit=flux[:2], label="g")
    fitter.add_flux_limits(band=BAND, t=T[:2], limit=band_flux[:2], sigma=band_flux[:2] / 4)
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
    assert loaded._limit_labels == ["g"]
    assert len(loaded._band_limits) == 1
    assert loaded._evaluate(_params()) == pytest.approx(chi2, rel=1e-12)


def test_limits_are_counted_in_the_data_count(truth_flux):
    """n_data counts limits with the detections and n_limits reports them separately."""
    flux, band_flux = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux)
    fitter.add_flux_density_limits(nu=5e14, t=T[:2], limit=flux[:2])
    fitter.add_flux_limits(band=BAND, t=T[:3], limit=band_flux[:3])
    n_data, n_limits = fitter._data_counts()
    assert (n_data, n_limits) == (15, 5)


def _seed(fitter):
    from VegasAfterglow import FitResult, ParamDef, Scale
    from VegasAfterglow.fitting.utils import _build_transformer

    defs = [ParamDef("E_iso", 1e50, 1e54, Scale.log), ParamDef("Gamma0", 10, 1000, Scale.log)]
    fitter._param_defs = defs
    fitter._to_params = _build_transformer(defs)
    fitter.result = FitResult(
        samples=np.zeros((4, 1, 2)), log_probs=np.zeros((4, 1)),
        labels=("log10_E_iso", "log10_Gamma0"),
        top_k_params=np.array([[53.0, np.log10(300.0)]]), top_k_log_probs=np.array([0.0]),
    )


def test_draw_fit_draws_limits(truth_flux):
    """draw_fit draws upper limits as arrows and a model curve for a limit-only frequency."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    flux, band_flux = truth_flux
    with_det = _fitter()
    with_det.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux)
    _seed(with_det)
    fig, (ax, _) = with_det.draw_fit(n_samples=0, show_nu_panel=False)
    n_lines_det_only = len(ax.lines)
    plt.close(fig)

    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux)
    fitter.add_flux_density_limits(nu=1e15, t=T[:2] / 3, limit=flux[:2])  # same band
    fitter.add_flux_density_limits(nu=1e9, t=T[:3], limit=flux[:3], label="VLA")  # limit-only band
    fitter.add_flux_limits(band=BAND, t=T[:2], limit=band_flux[:2])
    _seed(fitter)
    fig, (ax, _) = fitter.draw_fit(n_samples=0, show_nu_panel=False)
    try:
        assert len(ax.lines) > n_lines_det_only  # a curve for the limit-only frequency
        labels = [t.get_text() for t in ax.get_legend().get_texts()]
        assert any("VLA" in lbl for lbl in labels)
    finally:
        plt.close(fig)
