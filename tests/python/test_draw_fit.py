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


# ── Calibration systematics: the data are corrected, the model is not ───────

_CAL_DEFS = [
    ParamDef("E_iso", 1e50, 1e54, Scale.log),
    ParamDef("Gamma0", 10, 1000, Scale.log),
    ParamDef("n_ism", 1.0, 1.0, Scale.fixed),
]
_CAL_BEST = np.array([52.0, np.log10(300.0)])
_CAL_T = np.logspace(3, 6, 10)


def _curves(ax, style):
    return [ln for ln in ax.get_lines() if ln.get_linestyle() == style and ln.get_xdata().size > 50]


def _hollow_markers(ax):
    """The markers left where shifted points were measured."""
    return [ln for ln in ax.get_lines() if ln.get_markerfacecolor() == "none"]


def _error_bars(ax):
    """(flux, half-length) of every error-bar set on ``ax``, in drawing order; marker-less sets have flux None."""
    out = []
    for container in ax.containers:
        segments = container.lines[2][0].get_segments()
        half = np.array([(s[1, 1] - s[0, 1]) / 2 for s in segments])
        out.append((None if container.lines[0] is None else container.lines[0].get_ydata(), half))
    return out


def _calibrated_pair(sys_err, kind):
    """Two fitters identical except for ``sys_err``; the data are the model x1.3, so the factor is well defined."""
    from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet

    t = _CAL_T
    truth = Model(TophatJet(0.1, 1e52, 300), ISM(1.0), Observer(1e28, 0.1, 0.0), Radiation(0.1, 0.01, 2.3))
    pair = []
    for s in (0.0, sys_err):
        f = Fitter(z=0.1, lumi_dist=1e28, jet="tophat", medium="ism")
        if kind == "lc":
            flux = 1.3 * np.asarray(truth.flux_density(t, np.full_like(t, 1e14)).total)
            f.add_flux_density(nu=1e14, t=t, f_nu=flux, err=0.05 * flux, sys_err=s)
        else:
            flux = 1.3 * np.asarray(truth.flux(t, 1e17, 1e18, 5).total)
            f.add_flux(band=(1e17, 1e18), t=t, flux=flux, err=0.05 * flux, sys_err=s)
        _seed_fit_result(f, _CAL_DEFS, theta=_CAL_BEST)
        pair.append(f)
    return pair


def _draw_pair(sys_err, kind):
    plain, calib = _calibrated_pair(sys_err, kind)
    _, (ax_plain, _) = plain.draw_fit(n_samples=0, show_nu_panel=False)
    _, (ax_calib, _) = calib.draw_fit(n_samples=0, show_nu_panel=False)
    return calib, ax_plain, ax_calib


@pytest.mark.parametrize("kind,style", [("lc", "-"), ("band", "--")])
def test_calibrated_dataset_is_drawn_corrected_against_the_physical_model(kind, style):
    """With sys_err the model curve is untouched and the data are divided by their best-fit factor, a hollow marker staying where each was measured."""
    try:
        calib, ax_plain, ax_calib = _draw_pair(0.5, kind)
        factor = 1.0 + calib.calibration_offsets(calib._to_params(_CAL_BEST))[0].offset
        assert 1.2 < factor < 1.3  # data at 1.3x the model, pulled slightly toward 1 by the prior
        (plain_curve,) = _curves(ax_plain, style)
        (calib_curve,) = _curves(ax_calib, style)
        np.testing.assert_array_equal(calib_curve.get_ydata(), plain_curve.get_ydata())
        assert not _curves(ax_calib, ":")
        assert f"data $\\times${1 / factor:.3g}" in calib_curve.get_label()
        assert "data" not in plain_curve.get_label()

        ((measured, stat),) = _error_bars(ax_plain)
        ((drawn, drawn_err),) = _error_bars(ax_calib)
        np.testing.assert_allclose(drawn, measured / factor, rtol=1e-12)
        np.testing.assert_allclose(drawn_err, stat / factor, rtol=1e-12)
        (hollow,) = _hollow_markers(ax_calib)
        np.testing.assert_array_equal(hollow.get_ydata(), measured)
        assert not _hollow_markers(ax_plain)
    finally:
        plt.close("all")


def test_calibrated_envelope_stays_on_the_physical_model():
    """The confidence envelope belongs to the model, so it brackets the unscaled curve."""
    _, f = _calibrated_pair(0.5, "lc")
    rng = np.random.RandomState(0)
    f.result = FitResult(
        samples=(_CAL_BEST[None, :] + rng.normal(0, 0.05, (20, 2)))[:, None, :],
        log_probs=rng.normal(-30, 0.1, (20, 1)),
        labels=tuple(d.name for d in _CAL_DEFS if d.scale is not Scale.fixed),
        top_k_params=_CAL_BEST[None, :],
        top_k_log_probs=np.array([-30.0]),
    )
    try:
        _, (ax, _) = f.draw_fit(n_samples=20, ci=0.68, show_nu_panel=False)
        (curve,) = _curves(ax, "-")
        from matplotlib.collections import PolyCollection

        (band,) = [c for c in ax.collections if isinstance(c, PolyCollection)]
        verts = band.get_paths()[0].vertices[:, 1]
        y = curve.get_ydata()
        assert verts.min() <= y.min() * (1 + 1e-9) and verts.max() >= y.max() * (1 - 1e-9)
    finally:
        plt.close("all")


@pytest.mark.parametrize("kind,style", [("lc", "-"), ("band", "--")])
def test_each_sys_err_group_is_corrected_by_its_own_factor(kind, style):
    """Points with different sys_err are divided by their own group's factor and land on the one model curve."""
    sys_err = np.repeat([0.5, 0.1], 5)
    try:
        calib, ax_plain, ax_calib = _draw_pair(sys_err, kind)
        offsets = calib.calibration_offsets(calib._to_params(_CAL_BEST))
        factor = np.repeat([1.0 + o.offset for o in offsets], 5)
        assert factor[0] > factor[-1] > 1.0  # the looser prior lets its group move further
        (plain_curve,) = _curves(ax_plain, style)
        (calib_curve,) = _curves(ax_calib, style)
        np.testing.assert_array_equal(calib_curve.get_ydata(), plain_curve.get_ydata())
        lo, hi = 1 / factor.max(), 1 / factor.min()
        assert f"data $\\times${lo:.3g}–{hi:.3g}" in calib_curve.get_label()

        ((measured, stat),) = _error_bars(ax_plain)
        ((drawn, drawn_err),) = _error_bars(ax_calib)
        np.testing.assert_allclose(drawn, measured / factor, rtol=1e-12)
        np.testing.assert_allclose(drawn_err, stat / factor, rtol=1e-12)
    finally:
        plt.close("all")


def test_point_alone_with_its_sys_err_is_not_moved_and_shows_the_total_error():
    """A value carried by one point cannot be told from its residual: the point stays put and a second bar shows stat and sys combined."""
    sys_err = np.r_[np.full(8, 0.5), 0.2, 0.3]  # one group of eight, two lone points
    try:
        calib, ax_plain, ax_calib = _draw_pair(sys_err, "lc")
        factor = 1.0 + calib.calibration_offsets(calib._to_params(_CAL_BEST))[0].offset
        ((measured, stat),) = _error_bars(ax_plain)
        (no_marker, total), (drawn, drawn_err) = _error_bars(ax_calib)
        assert no_marker is None
        np.testing.assert_allclose(total, np.hypot(stat[8:], sys_err[8:] * measured[8:]), rtol=1e-12)
        np.testing.assert_allclose(drawn[:8], measured[:8] / factor, rtol=1e-12)
        np.testing.assert_array_equal(drawn[8:], measured[8:])
        np.testing.assert_array_equal(drawn_err[8:], stat[8:])
        (hollow,) = _hollow_markers(ax_calib)
        np.testing.assert_array_equal(hollow.get_ydata(), measured[:8])
    finally:
        plt.close("all")


def test_constant_sys_err_array_draws_like_the_scalar():
    """An array holding one value is one calibration factor, drawn exactly as the scalar is."""
    try:
        _, _, ax_scalar = _draw_pair(0.5, "lc")
        _, _, ax_array = _draw_pair(np.full(10, 0.5), "lc")
        (scalar_curve,) = _curves(ax_scalar, "-")
        (array_curve,) = _curves(ax_array, "-")
        assert array_curve.get_label() == scalar_curve.get_label()
        ((scalar_drawn, _),) = _error_bars(ax_scalar)
        ((array_drawn, _),) = _error_bars(ax_array)
        np.testing.assert_array_equal(array_drawn, scalar_drawn)
    finally:
        plt.close("all")
