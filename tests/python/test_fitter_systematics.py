"""Per-dataset calibration systematics in the log-flux likelihood.

``sys_err`` on an ``add_*`` call is the fractional calibration uncertainty of
that dataset: one shared multiplicative factor with a Gaussian prior of that
width in ln-flux, marginalised analytically. Within-dataset precision is kept;
only the dataset's overall level is free to shift.
"""
import numpy as np
import pytest
from scipy.optimize import minimize_scalar

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


def _profiled_chi2(log_obs, log_err, model_flux, sys_err):
    """Reference: minimise over the shared ln-offset c with its Gaussian prior."""
    r = log_obs - np.log(model_flux)
    prec = 1.0 / log_err**2

    def f(c):
        return np.sum(prec * (r - c) ** 2) + c**2 / sys_err**2

    return minimize_scalar(f, bounds=(-2, 2), method="bounded", options={"xatol": 1e-12}).fun


def test_sys_err_zero_is_plain_chi2(truth_flux):
    """With sys_err=0 (the default) the chi2 is the plain per-point sum."""
    flux, _ = truth_flux
    obs = 1.3 * flux
    err = 0.1 * obs
    f0, f1 = _fitter(), _fitter()
    f0.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err)
    f1.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err, sys_err=0.0)
    assert f1._evaluate(_params()) == pytest.approx(f0._evaluate(_params()), rel=1e-12)


def test_sys_err_marginalises_a_shared_offset(truth_flux):
    """chi2 with sys_err equals the profile over one shared ln-flux offset with a Gaussian prior of width sys_err."""
    flux, _ = truth_flux
    obs = 1.3 * flux
    err = 0.1 * obs
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err, sys_err=0.2)
    chi2 = fitter._evaluate(_params())
    model_flux = np.asarray(fitter._build_model(_params()).flux_density(T, NU).total)
    expected = _profiled_chi2(np.log(obs), err / obs, model_flux, 0.2)
    assert chi2 == pytest.approx(expected, rel=1e-8)


def test_large_sys_err_makes_the_dataset_level_free(truth_flux):
    """With a very loose calibration prior, scaling the whole dataset leaves chi2 unchanged."""
    flux, _ = truth_flux
    err = 0.1 * flux
    f1, f2 = _fitter(), _fitter()
    f1.add_flux_density(nu=1e15, t=T, f_nu=flux, err=err, sys_err=100.0)
    f2.add_flux_density(nu=1e15, t=T, f_nu=3.0 * flux, err=3.0 * err, sys_err=100.0)
    # Only the prior's own pull remains: (ln 3)^2 / 100^2 ~ 1e-4 in chi2.
    assert f2._evaluate(_params()) == pytest.approx(f1._evaluate(_params()), abs=1e-3)


def test_sys_err_applies_per_dataset(truth_flux):
    """Each add call is its own dataset: a sys_err on one call leaves the other call's chi2 untouched."""
    flux, _ = truth_flux
    obs = 1.3 * flux
    err = 0.1 * obs
    plain = _fitter()
    plain.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err)
    chi2_plain = plain._evaluate(_params())

    both = _fitter()
    both.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err)
    both.add_flux_density(nu=5e14, t=T, f_nu=obs, err=err, sys_err=0.2)
    marg = _fitter()
    marg.add_flux_density(nu=5e14, t=T, f_nu=obs, err=err, sys_err=0.2)
    assert both._evaluate(_params()) == pytest.approx(
        chi2_plain + marg._evaluate(_params()), rel=1e-10
    )


def test_sys_err_on_band_integrated_flux(truth_flux):
    """add_flux accepts sys_err and marginalises the band dataset's level the same way."""
    _, band_flux = truth_flux
    obs = 1.3 * band_flux
    err = 0.1 * obs
    fitter = _fitter()
    fitter.add_flux(band=BAND, t=T, flux=obs, err=err, sys_err=0.2)
    chi2 = fitter._evaluate(_params())
    model_flux = np.asarray(fitter._build_model(_params()).flux(T, BAND[0], BAND[1], 5).total)
    expected = _profiled_chi2(np.log(obs), err / obs, model_flux, 0.2)
    assert chi2 == pytest.approx(expected, rel=1e-8)


def test_calibration_offsets_reports_the_fitted_shift(truth_flux):
    """calibration_offsets() returns each dataset's best-fit level relative to the model, pulled slightly toward zero by the prior."""
    flux, _ = truth_flux
    obs = 1.05 * flux
    err = 0.1 * obs
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=obs, err=err, sys_err=0.5, label="r")
    fitter.add_flux_density(nu=5e14, t=T, f_nu=flux, err=0.1 * flux)
    offsets = fitter.calibration_offsets(_params())
    assert [o.label for o in offsets] == ["r"]
    prec = len(T) / 0.1**2  # 10 points at 10% each
    expected = np.expm1(np.log(1.05) * prec / (1.0 / 0.5**2 + prec))
    assert offsets[0].offset == pytest.approx(expected, rel=1e-6)
    assert offsets[0].sys_err == 0.5


@pytest.mark.parametrize("bad", [-0.1, np.nan, np.inf])
def test_sys_err_must_be_finite_and_non_negative(bad):
    fitter = _fitter()
    with pytest.raises(ValueError, match="sys_err"):
        fitter.add_flux_density(nu=1e15, t=T, f_nu=np.ones_like(T), err=0.1 * np.ones_like(T), sys_err=bad)


def _seed_result(fitter):
    from VegasAfterglow import FitResult, ParamDef, Scale
    from VegasAfterglow.fitting.utils import _build_transformer

    defs = [ParamDef("E_iso", 1e50, 1e54, Scale.log), ParamDef("Gamma0", 10, 1000, Scale.log)]
    fitter._param_defs = defs
    fitter._to_params = _build_transformer(defs)
    fitter.result = FitResult(
        samples=np.zeros((4, 1, 2)),
        log_probs=np.zeros((4, 1)),
        labels=("log10_E_iso", "log10_Gamma0"),
        top_k_params=np.array([[53.0, np.log10(300.0)]]),
        top_k_log_probs=np.array([0.0]),
    )


def test_sys_err_survives_save_and_load(tmp_path, truth_flux):
    """sys_err of point and band datasets round-trips through save/load, and the loaded fitter evaluates."""
    flux, band_flux = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux, sys_err=0.07, label="r")
    fitter.add_flux(band=BAND, t=T, flux=band_flux, err=0.1 * band_flux, sys_err=0.12)
    _seed_result(fitter)
    chi2 = fitter._evaluate(_params())
    path = tmp_path / "fit.h5"
    fitter.save(path)
    loaded = Fitter.load(path)
    np.testing.assert_array_equal(loaded._point_sys, [np.full(len(T), 0.07)])
    np.testing.assert_array_equal(loaded._band_obs[0].sys_err, np.full(len(T), 0.12))
    assert loaded._evaluate(_params()) == pytest.approx(chi2, rel=1e-12)


def test_summary_reports_calibration_offsets():
    """FitResult.summary() lists each dataset's fitted calibration shift with its prior width."""
    from VegasAfterglow import FitResult

    r = FitResult(
        samples=np.zeros((4, 1, 1)),
        log_probs=np.zeros((4, 1)),
        labels=("log10_E_iso",),
        top_k_params=np.array([[53.0]]),
        top_k_log_probs=np.array([-5.0]),
        n_data=20,
        n_free_params=1,
        calibration_offsets=[{"label": "XRT", "offset": 0.031, "sys_err": 0.1}],
    )
    text = str(r.summary())
    assert "XRT" in text
    assert "+3.1%" in text
    assert "10%" in text


def test_fit_reports_calibration_and_limits(tmp_path, truth_flux):
    """A real fit populates n_limits and calibration_offsets on the result, prints them, and keeps them through save/load."""
    from VegasAfterglow import ParamDef, Scale

    flux, _ = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=1.05 * flux, err=0.1 * flux, sys_err=0.1, label="r")
    fitter.add_flux_density_limits(nu=5e14, t=T[:2], limit=flux[:2])
    defs = [
        ParamDef("E_iso", 1e52, 1e54, Scale.log),
        ParamDef("Gamma0", 300, 300, Scale.fixed),
        ParamDef("theta_c", 0.1, 0.1, Scale.fixed),
        ParamDef("n_ism", 1.0, 1.0, Scale.fixed),
        ParamDef("eps_e", 0.1, 0.1, Scale.fixed),
        ParamDef("eps_B", 1e-3, 1e-3, Scale.fixed),
        ParamDef("p", 2.3, 2.3, Scale.fixed),
    ]
    result = fitter.fit(defs, sampler="emcee", nsteps=12, nburn=4, npool=1,
                        outdir=str(tmp_path), top_k=3)
    assert (result.n_data, result.n_limits) == (12, 2)
    assert len(result.calibration_offsets) == 1
    assert result.calibration_offsets[0]["label"] == "r"
    assert abs(result.calibration_offsets[0]["offset"]) < 0.5
    text = str(result.summary())
    assert "Calibration offsets" in text and "2 upper limits" in text
    path = tmp_path / "fit.h5"
    fitter.save(path)
    loaded = Fitter.load(path).result
    assert loaded.n_limits == 2
    assert loaded.calibration_offsets == result.calibration_offsets


def test_sys_err_on_spectrum(truth_flux):
    """add_spectrum accepts sys_err and marginalises the spectrum's level."""
    flux, _ = truth_flux
    obs = 1.3 * flux[:1]
    err = 0.1 * obs
    fitter = _fitter()
    fitter.add_spectrum(t=T[0], nu=[1e15], f_nu=obs, err=err, sys_err=0.2)
    chi2 = fitter._evaluate(_params())
    model_flux = np.asarray(fitter._build_model(_params()).flux_density(T[:1], NU[:1]).total)
    assert chi2 == pytest.approx(_profiled_chi2(np.log(obs), err / obs, model_flux, 0.2), rel=1e-8)


def test_load_fit_saved_before_sys_err_and_limits(tmp_path, truth_flux):
    """A snapshot written by an older version (no sys_err, no limit entries) loads with defaults and evaluates."""
    import json

    import bilby

    flux, band_flux = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux)
    fitter.add_flux(band=BAND, t=T, flux=band_flux, err=0.1 * band_flux)
    _seed_result(fitter)
    chi2 = fitter._evaluate(_params())
    path = tmp_path / "fit.h5"
    fitter.save(path)

    # Strip everything this release added from the saved snapshot.
    br = bilby.read_in_result(filename=str(path))
    vg = br.meta_data["vegasafterglow"]
    cfg = json.loads(vg["fitter_config_json"])
    cfg["version"] = 1
    for pt in cfg["point_data"]:
        del pt["sys_err"]
    for bd in cfg["band_data"]:
        del bd["sys_err"]
    del cfg["limit_data"], cfg["band_limit_data"]
    vg["fitter_config_json"] = json.dumps(cfg)
    for key in ("n_limits", "calibration_offsets_json"):
        vg.pop(key, None)
    old = tmp_path / "old.h5"
    br.save_to_file(filename=str(old))

    loaded = Fitter.load(old)
    np.testing.assert_array_equal(loaded._point_sys, [np.zeros(len(T))])
    np.testing.assert_array_equal(loaded._band_obs[0].sys_err, np.zeros(len(T)))
    assert loaded._limit_t == [] and loaded._band_limits == []
    assert loaded.result.n_limits is None and loaded.result.calibration_offsets is None
    assert loaded._evaluate(_params()) == pytest.approx(chi2, rel=1e-12)


# ── sys_err given per point: equal values share one factor ──────────────────

NU_SED = np.logspace(9, 18, 10)
T_SED = 1e5
# Interleaved on purpose: grouping is by value, not by position.
SYS = np.array([0.2, 0.05, 0.2, 0.0, 0.05, 0.2, 0.0, 0.05, 0.2, 0.0])
LEVEL = np.where(SYS == 0.2, 1.3, np.where(SYS == 0.05, 0.9, 1.0))


def _dataset(kind):
    """(times, frequencies, model flux) of a 10-point dataset of each kind, from the fitter's own model."""
    model = _fitter()._build_model(_params())
    if kind == "band":
        return T, None, np.asarray(model.flux(T, BAND[0], BAND[1], 5).total)
    t, nu = (T, NU) if kind == "lc" else (np.full_like(NU_SED, T_SED), NU_SED)
    return t, nu, np.asarray(model.flux_density(t, nu).total)


def _add(fitter, kind, obs, err, sys_err, pick=slice(None)):
    if kind == "lc":
        fitter.add_flux_density(nu=1e15, t=T[pick], f_nu=obs[pick], err=err[pick], sys_err=sys_err)
    elif kind == "spectrum":
        fitter.add_spectrum(t=T_SED, nu=NU_SED[pick], f_nu=obs[pick], err=err[pick], sys_err=sys_err)
    else:
        fitter.add_flux(band=BAND, t=T[pick], flux=obs[pick], err=err[pick], sys_err=sys_err)


KINDS = ["lc", "spectrum", "band"]


@pytest.mark.parametrize("kind", KINDS)
def test_constant_sys_err_array_equals_scalar(kind):
    """An array holding one value everywhere is the same likelihood as that scalar."""
    _, _, flux = _dataset(kind)
    obs, err = 1.3 * flux, 0.13 * flux
    scalar, array = _fitter(), _fitter()
    _add(scalar, kind, obs, err, 0.2)
    _add(array, kind, obs, err, np.full(len(obs), 0.2))
    assert array._evaluate(_params()) == pytest.approx(scalar._evaluate(_params()), rel=1e-12)


@pytest.mark.parametrize("kind", KINDS)
def test_sys_err_array_groups_points_by_value(kind):
    """Points given the same value share one marginalised offset; each value is its own offset; zeros are plain."""
    _, _, flux = _dataset(kind)
    obs = LEVEL * flux
    err = 0.1 * obs
    fitter = _fitter()
    _add(fitter, kind, obs, err, SYS)
    expected = sum(
        _profiled_chi2(np.log(obs[SYS == s]), 0.1, flux[SYS == s], s) for s in (0.2, 0.05)
    ) + np.sum((np.log(obs / flux)[SYS == 0.0] / 0.1) ** 2)
    assert fitter._evaluate(_params()) == pytest.approx(expected, rel=1e-8)


@pytest.mark.parametrize("kind", KINDS)
def test_distinct_sys_err_values_add_in_quadrature(kind):
    """When every point has its own value, each one is its own offset: sys_err adds to the ln-flux error in quadrature."""
    _, _, flux = _dataset(kind)
    obs = flux * np.linspace(0.8, 1.3, len(flux))
    err = 0.1 * obs
    sys_err = np.linspace(0.02, 0.3, len(flux))
    fitter = _fitter()
    _add(fitter, kind, obs, err, sys_err)
    expected = np.sum(np.log(obs / flux) ** 2 / (0.1**2 + sys_err**2))
    assert fitter._evaluate(_params()) == pytest.approx(expected, rel=1e-10)


@pytest.mark.parametrize("kind", ["lc", "spectrum"])
def test_sys_err_array_equals_one_call_per_value(kind):
    """One call with an array is the same likelihood as one call per instrument with scalars."""
    _, _, flux = _dataset(kind)
    obs = LEVEL * flux
    err = 0.1 * obs
    array, split = _fitter(), _fitter()
    _add(array, kind, obs, err, SYS)
    for s in (0.2, 0.05, 0.0):
        _add(split, kind, obs, err, s, pick=SYS == s)
    assert array._evaluate(_params()) == pytest.approx(split._evaluate(_params()), rel=1e-12)


def test_equal_sys_err_in_separate_calls_stays_independent():
    """Two instruments with the same sys_err are tied inside one call and independent when added separately."""
    _, _, flux = _dataset("spectrum")
    obs = flux * np.where(np.arange(len(flux)) < 5, 1.3, 0.8)  # the two halves disagree on the level
    err = 0.1 * obs
    tied, apart = _fitter(), _fitter()
    _add(tied, "spectrum", obs, err, 0.2)
    _add(apart, "spectrum", obs, err, 0.2, pick=slice(0, 5))
    _add(apart, "spectrum", obs, err, 0.2, pick=slice(5, None))
    expected = sum(
        _profiled_chi2(np.log(obs[half]), 0.1, flux[half], 0.2) for half in (slice(0, 5), slice(5, None))
    )
    assert apart._evaluate(_params()) == pytest.approx(expected, rel=1e-8)
    assert tied._evaluate(_params()) > 2 * apart._evaluate(_params())


def test_sys_err_array_follows_the_time_sort_of_band_data():
    """add_flux sorts its points by time; each sys_err value stays with its point."""
    _, _, flux = _dataset("band")
    obs = LEVEL * flux
    err = 0.1 * obs
    shuffle = np.random.RandomState(3).permutation(len(T))
    ordered, shuffled = _fitter(), _fitter()
    ordered.add_flux(band=BAND, t=T, flux=obs, err=err, sys_err=SYS)
    shuffled.add_flux(band=BAND, t=T[shuffle], flux=obs[shuffle], err=err[shuffle], sys_err=SYS[shuffle])
    assert shuffled._evaluate(_params()) == pytest.approx(ordered._evaluate(_params()), rel=1e-12)


@pytest.mark.parametrize("kind", KINDS)
def test_calibration_offsets_has_one_entry_per_sys_err_value(kind):
    """Each distinct non-zero value of a call is reported once, in order of first appearance, under the call's label."""
    _, _, flux = _dataset(kind)
    obs = LEVEL * flux
    err = 0.1 * obs
    fitter = _fitter()
    _add(fitter, kind, obs, err, SYS)
    offsets = fitter.calibration_offsets(_params())
    assert [o.sys_err for o in offsets] == [0.2, 0.05]
    assert len({o.label for o in offsets}) == 1
    for o, level in zip(offsets, (1.3, 0.9)):
        prec = np.sum(SYS == o.sys_err) / 0.1**2
        expected = np.expm1(np.log(level) * prec / (1.0 / o.sys_err**2 + prec))
        assert o.offset == pytest.approx(expected, rel=1e-6)


@pytest.mark.parametrize("kind", KINDS)
def test_sys_err_array_must_have_one_value_per_point(kind):
    _, _, flux = _dataset(kind)
    with pytest.raises(ValueError, match="sys_err must be a scalar or one value per point"):
        _add(_fitter(), kind, flux, 0.1 * flux, np.full(len(flux) - 1, 0.1))


@pytest.mark.parametrize("bad", [-0.1, np.nan, np.inf])
@pytest.mark.parametrize("kind", KINDS)
def test_sys_err_array_entries_must_be_finite_and_non_negative(kind, bad):
    _, _, flux = _dataset(kind)
    sys_err = np.full(len(flux), 0.1)
    sys_err[4] = bad
    with pytest.raises(ValueError, match="sys_err must be finite and >= 0"):
        _add(_fitter(), kind, flux, 0.1 * flux, sys_err)


def test_sys_err_array_survives_save_and_load(tmp_path, truth_flux):
    """Per-point sys_err of point and band datasets round-trips through save/load, and the loaded fitter evaluates."""
    flux, band_flux = truth_flux
    fitter = _fitter()
    fitter.add_flux_density(nu=1e15, t=T, f_nu=flux, err=0.1 * flux, sys_err=SYS, label="r")
    fitter.add_flux(band=BAND, t=T, flux=band_flux, err=0.1 * band_flux, sys_err=SYS[::-1])
    _seed_result(fitter)
    chi2 = fitter._evaluate(_params())
    path = tmp_path / "fit.h5"
    fitter.save(path)
    loaded = Fitter.load(path)
    np.testing.assert_array_equal(loaded._point_sys, [SYS])
    np.testing.assert_array_equal(loaded._band_obs[0].sys_err, SYS[::-1])
    assert loaded._evaluate(_params()) == pytest.approx(chi2, rel=1e-12)
