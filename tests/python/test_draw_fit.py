"""Smoke tests for ``Fitter.draw_fit``.

We synthesise a minimal ``FitResult`` directly on the Fitter (no MCMC) so the
tests stay fast and don't depend on a real sampler run. The model is invoked
through the real C++ pipeline -- these tests therefore also serve as a smoke
check that ``flux_density_grid`` / ``flux`` / ``model.details`` survive the
typical post-fit usage pattern.

Layout tests use ``n_samples=0`` to skip the credible-band computation
(the band itself is exercised by the dedicated ``test_credible_band_path``
test with a tight, realistic posterior fixture).
"""

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from VegasAfterglow import FitResult, Fitter, ParamDef, Scale  # noqa: E402
from VegasAfterglow.fitting.utils import _build_transformer  # noqa: E402


def _seed_fit_result(fitter, defs, theta=(52.0, np.log10(300.0), 2.3)):
    """Attach a minimal FitResult so post-fit methods (.flux_density_grid etc.)
    can run without an actual MCMC."""
    fitter._param_defs = list(defs)
    fitter._to_params = _build_transformer(list(defs))
    fitter.result = FitResult(
        samples=np.zeros((10, 1, len(defs))),
        log_probs=np.zeros((10, 1)),
        labels=tuple(d.name for d in defs),
        top_k_params=np.array([list(theta)]),
        top_k_log_probs=np.array([0.0]),
    )


def _default_defs():
    return [
        ParamDef("E_iso", 1e50, 1e54, Scale.log),
        ParamDef("Gamma0", 10, 1000, Scale.log),
        ParamDef("p", 2.0, 3.0, Scale.linear),
    ]


def _mixed_fitter():
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    t = np.logspace(3, 6, 20)
    f.add_flux_density(
        nu=1e9, t=t, f_nu=np.full_like(t, 1e-28), err=np.full_like(t, 1e-29)
    )
    f.add_flux_density(
        nu=1e14, t=t, f_nu=np.full_like(t, 1e-27), err=np.full_like(t, 1e-28)
    )
    f.add_flux(
        band=(1e17, 1e18),
        t=t,
        flux=np.full_like(t, 1e-12),
        err=np.full_like(t, 1e-13),
    )
    _seed_fit_result(f, _default_defs())
    return f


def test_mixed_returns_three_axes():
    """LC + band-integrated -> 3 axes (left, twin, bottom) and equal decade
    span on the dual y-axes."""
    fitter = _mixed_fitter()
    try:
        fig, (ax_top, ax_bot) = fitter.draw_fit(n_samples=0)
        assert ax_bot is not None
        assert len(fig.axes) == 3, f"expected 3 axes, got {len(fig.axes)}"
        # Identify the twin: the other axis on the top panel that is not ax_bot.
        twin_candidates = [a for a in fig.axes if a is not ax_top and a is not ax_bot]
        assert len(twin_candidates) == 1
        ax_R = twin_candidates[0]
        spanL = np.log10(ax_top.get_ylim()[1] / ax_top.get_ylim()[0])
        spanR = np.log10(ax_R.get_ylim()[1] / ax_R.get_ylim()[0])
        assert abs(spanL - spanR) < 1e-6, (
            f"decade spans differ: left={spanL:.6g}, right={spanR:.6g}"
        )
    finally:
        plt.close("all")


def test_lc_only_no_twinx():
    """Light-curve-only path -> 2 axes (top + bottom), no twin axis."""
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    t = np.logspace(3, 6, 10)
    f.add_flux_density(
        nu=1e9, t=t, f_nu=np.full_like(t, 1e-28), err=np.full_like(t, 1e-29)
    )
    _seed_fit_result(f, _default_defs())
    try:
        fig, (ax_top, ax_bot) = f.draw_fit(n_samples=0)
        assert ax_bot is not None
        assert len(fig.axes) == 2, f"expected 2 axes, got {len(fig.axes)}"
    finally:
        plt.close("all")


def test_band_only_no_twinx():
    """Band-only path -> 2 axes (top + bottom), single y-axis on top."""
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    t = np.logspace(3, 6, 10)
    f.add_flux(
        band=(1e17, 1e18),
        t=t,
        flux=np.full_like(t, 1e-12),
        err=np.full_like(t, 1e-13),
    )
    _seed_fit_result(f, _default_defs())
    try:
        fig, (ax_top, ax_bot) = f.draw_fit(n_samples=0)
        assert ax_bot is not None
        assert len(fig.axes) == 2
    finally:
        plt.close("all")


def test_no_nu_panel():
    """show_nu_panel=False -> single panel, ax_bot is None."""
    try:
        fig, (ax_top, ax_bot) = _mixed_fitter().draw_fit(
            n_samples=0, show_nu_panel=False
        )
        assert ax_bot is None
        # Top axis + its twin = 2 axes total (no bottom).
        assert len(fig.axes) == 2
    finally:
        plt.close("all")


def test_no_data_raises():
    """Fitter with no observation data -> ValueError."""
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    _seed_fit_result(f, _default_defs())
    with pytest.raises(ValueError, match="no observation data"):
        f.draw_fit(n_samples=0)


def test_explicit_best_params_works_without_fit():
    """Passing best_params explicitly should bypass the result-requirement."""
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    t = np.logspace(3, 6, 10)
    f.add_flux_density(
        nu=1e14, t=t, f_nu=np.full_like(t, 1e-27), err=np.full_like(t, 1e-28)
    )
    defs = _default_defs()
    f._param_defs = defs
    f._to_params = _build_transformer(defs)
    # No fitter.result attached -- explicit best_params should still work.
    try:
        fig, _ = f.draw_fit(
            best_params=np.array([52.0, np.log10(300.0), 2.3]), n_samples=0
        )
        assert fig is not None
    finally:
        plt.close("all")


def _fitter_with_realistic_posterior():
    f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
    t = np.logspace(3, 6, 10)
    f.add_flux_density(
        nu=1e14, t=t, f_nu=np.full_like(t, 1e-27), err=np.full_like(t, 1e-28)
    )
    defs = _default_defs()
    f._param_defs = defs
    f._to_params = _build_transformer(defs)
    rng = np.random.RandomState(0)
    best = np.array([52.0, np.log10(300.0), 2.3])
    samples = (best[None, :] + rng.normal(0, 0.05, (20, 3)))[:, None, :]
    f.result = FitResult(
        samples=samples,
        log_probs=rng.normal(-30, 0.1, (20, 1)),
        labels=tuple(d.name for d in defs),
        top_k_params=best[None, :],
        top_k_log_probs=np.array([-30.0]),
    )
    return f


def _has_polycollection(ax):
    from matplotlib.collections import PolyCollection
    return any(isinstance(c, PolyCollection) for c in ax.collections)


def test_credible_band_path():
    """With a realistic posterior and n_samples>0, draw_fit renders at least one credible-band fill (PolyCollection) on the top axis."""
    try:
        fig, (ax_top, _) = _fitter_with_realistic_posterior().draw_fit(
            n_samples=8, show_nu_panel=False
        )
        assert _has_polycollection(ax_top)
    finally:
        plt.close("all")


@pytest.mark.parametrize("obs_noise", ["frac", "abs", "none"])
def test_obs_noise_modes_render(obs_noise):
    """All three obs_noise modes render a fill_between band without error."""
    f = _fitter_with_realistic_posterior()
    try:
        fig, (ax, _) = f.draw_fit(
            n_samples=8, show_nu_panel=False, obs_noise=obs_noise
        )
        assert _has_polycollection(ax)
    finally:
        plt.close("all")


def test_obs_noise_invalid_value_raises():
    """Unknown obs_noise string surfaces a clear ValueError listing the valid set."""
    f = _fitter_with_realistic_posterior()
    with pytest.raises(ValueError, match="obs_noise must be one of"):
        f.draw_fit(n_samples=8, show_nu_panel=False, obs_noise="frac_typo")
