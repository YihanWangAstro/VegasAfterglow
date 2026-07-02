"""FitResult.summary() formatting + Fitter.save / Fitter.load bilby-native roundtrip."""
import os
import tempfile

import numpy as np

from VegasAfterglow import FitResult, Fitter, ParamDef, Scale


def _make_result(top_k=2, latex_labels=("$E$", r"$\Gamma_0$", "$p$")):
    return FitResult(
        samples=np.arange(100 * 4 * 3, dtype=float).reshape(100, 4, 3),
        log_probs=np.linspace(-100, -10, 100 * 4).reshape(100, 4),
        labels=("log10_E_iso", "Gamma0", "p"),
        latex_labels=latex_labels,
        top_k_params=(
            np.array([[52.3, 287.4, 2.31], [52.1, 250.0, 2.30]])[:top_k]
            if top_k else None
        ),
        top_k_log_probs=(
            np.array([-43.21, -45.07])[:top_k] if top_k else None
        ),
    )


def _seed_fitter(result):
    """Build a minimal Fitter pre-populated with a FitResult + ParamDefs so the
    on-disk save/load roundtrip can be exercised without running a real MCMC."""
    fitter = Fitter(z=0.5, lumi_dist=1e28, jet="tophat", medium="ism")
    fitter.add_flux_density(
        nu=5e14,
        t=np.array([1e3, 1e4, 1e5]),
        f_nu=np.array([1e-26, 1e-26, 1e-26]),
        err=np.array([1e-28, 1e-28, 1e-28]),
    )
    fitter._param_defs = [
        ParamDef("E_iso", 1e50, 1e54, Scale.log),
        ParamDef("Gamma0", 10, 1000, Scale.log),
        ParamDef("p", 2.0, 3.0, Scale.linear),
    ]
    fitter.result = result
    return fitter


# ----------------- summary() -----------------


def test_summary_populated_table():
    r = _make_result()
    out = str(r.summary())
    assert "Rank" in out and "chi^2" in out
    for lbl in r.labels:
        assert lbl in out
    # Auto-format: log10_E_iso ~52 -> .2f; Gamma0 ~287 -> .1f; p ~2.31 -> .3f
    assert "52.30" in out and "287.4" in out and "2.310" in out
    # chi^2 = -2 * -43.21 = 86.42
    assert "86.42" in out


def test_summary_title_shows_rows_and_total():
    r = _make_result()
    out = str(r.summary())
    assert "Best-fit summary (top 2 of 2)" in out


def test_summary_includes_fit_quality_header_when_populated():
    """n_data + n_free_params populated -> header line with χ²/DOF, BIC, AIC."""
    r = _make_result()
    r.n_data = 50
    r.n_free_params = 3
    out = str(r.summary())
    # chi²_min = -2 * max(top_k_log_probs) = -2 * -43.21 = 86.42; DOF = 50 - 3 = 47
    assert "χ²_min / DOF = 86.42 / 47" in out
    assert "BIC =" in out and "AIC =" in out


def test_summary_omits_fit_quality_header_when_dof_invalid():
    """Suppress header silently if n_data <= n_free_params (over-parametrized)."""
    r = _make_result()
    r.n_data = 2
    r.n_free_params = 3
    out = str(r.summary())
    assert "DOF" not in out
    assert "BIC" not in out


def test_summary_omits_fit_quality_header_when_unpopulated():
    """Older saved files have n_data=None; header suppressed."""
    r = _make_result()
    # n_data / n_free_params left as default None
    out = str(r.summary())
    assert "DOF" not in out


def test_summary_empty_top_k():
    r = FitResult(
        samples=np.zeros((10, 1, 2)),
        log_probs=np.zeros((10, 1)),
        labels=("a", "b"),
    )
    assert "no top_k_params stored" in str(r.summary())


def test_summary_top_k_kwarg():
    r = _make_result()
    out = str(r.summary(top_k=1))
    # First row present, second pruned
    assert "287.4" in out
    assert "250.0" not in out
    # Title reflects the slice
    assert "top 1 of 2" in out


def test_summary_includes_latex_block_when_labels_set():
    """latex_labels set => single LaTeX block (posterior median, asymmetric 1σ)."""
    r = _make_result()
    out = str(r.summary())
    # Single block, header references median / corner-plot semantics
    assert "LaTeX (median" in out
    # Not per-rank anymore
    assert "LaTeX (Rank" not in out
    # latex_labels = ("$E$", r"$\Gamma_0$", "$p$") — dollars stripped, re-wrapped
    assert "$E = " in out
    assert r"$\Gamma_0 = " in out
    assert "$p = " in out
    # Asymmetric bound markup
    assert "^{+" in out and "}_{-" in out


def test_summary_no_latex_when_labels_missing():
    """Auto mode is silent without latex_labels (no LaTeX block, no note).
    Explicit latex=True surfaces a note explaining how to enable the block."""
    r = FitResult(
        samples=np.arange(100 * 4 * 3, dtype=float).reshape(100, 4, 3),
        log_probs=np.linspace(-100, -10, 100 * 4).reshape(100, 4),
        labels=("a", "b", "c"),
        latex_labels=None,
        top_k_params=np.array([[1.0, 2.0, 3.0]]),
        top_k_log_probs=np.array([-43.21]),
    )
    # Auto: stays quiet
    auto_out = str(r.summary())
    assert "LaTeX" not in auto_out
    assert "no latex_labels" not in auto_out
    # Explicit opt-in with no labels: surface a note
    explicit_out = str(r.summary(latex=True))
    assert "no latex_labels" in explicit_out
    assert "LaTeX (Rank" not in explicit_out


def test_summary_latex_false_suppresses_block():
    """latex=False forces the block off even when latex_labels are set."""
    r = _make_result()
    out = str(r.summary(latex=False))
    assert "LaTeX" not in out
    assert "Rank" in out  # table itself still present


def test_summary_renders_via_repr_for_notebook():
    """The _SummaryTable wrapper's __repr__ returns the text directly so
    Jupyter last-line auto-display renders cleanly (no escaped \\n)."""
    r = _make_result()
    s = r.summary()
    assert repr(s) == str(s)
    assert "\\n" not in repr(s)  # not a Python-escaped string repr
    assert "Rank" in repr(s)


# ----------------- Fitter.save / Fitter.load roundtrip -----------------


def test_save_load_h5_roundtrip():
    fitter = _seed_fitter(_make_result())
    r = fitter.result
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fit.h5")
        fitter.save(path)
        assert os.path.exists(path)

        f2 = Fitter.load(path)
        r2 = f2.result
        np.testing.assert_allclose(r.samples, r2.samples)
        np.testing.assert_allclose(r.log_probs, r2.log_probs)
        np.testing.assert_allclose(r.top_k_params, r2.top_k_params)
        np.testing.assert_allclose(r.top_k_log_probs, r2.top_k_log_probs)
        assert tuple(r.labels) == tuple(r2.labels)
        assert tuple(r.latex_labels) == tuple(r2.latex_labels)
        # bilby_result is populated post-load (was None pre-save on emcee path)
        assert r2.bilby_result is not None
        # summary() output identical
        assert str(r.summary()) == str(r2.summary())
        # Fitter snapshot survived: jet/medium/z/lumi_dist and observation data.
        assert f2.jet == fitter.jet
        assert f2.medium == fitter.medium
        assert f2.z == fitter.z
        assert f2.lumi_dist == fitter.lumi_dist
        np.testing.assert_allclose(f2._point_t[0], fitter._point_t[0])
        # Fit-quality fields default to None in the test fixture; both sides
        # should carry that through.
        assert r2.n_data is None
        assert r2.n_free_params is None


def test_save_load_json_roundtrip():
    fitter = _seed_fitter(_make_result())
    r = fitter.result
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fit.json")
        fitter.save(path)
        assert os.path.exists(path)
        r2 = Fitter.load(path).result
        np.testing.assert_allclose(r.samples, r2.samples)
        np.testing.assert_allclose(r.top_k_params, r2.top_k_params)
        assert tuple(r.labels) == tuple(r2.labels)
        assert r2.n_data is None
        assert r2.n_free_params is None


def test_save_load_roundtrip_preserves_fit_quality_fields():
    """When the FitResult has n_data / n_free_params populated (as it does
    after a real fit() call), both fields must survive HDF5 and JSON
    round-trip and the summary() header must remain identical."""
    r = _make_result()
    r.n_data = 50
    r.n_free_params = 3
    fitter = _seed_fitter(r)
    pre_summary = str(fitter.result.summary())
    assert "DOF" in pre_summary  # sanity: header is present pre-save

    with tempfile.TemporaryDirectory() as d:
        for ext in ("h5", "json"):
            path = os.path.join(d, f"fit.{ext}")
            fitter.save(path)
            r2 = Fitter.load(path).result
            assert r2.n_data == 50, f"{ext}: n_data lost on round-trip"
            assert r2.n_free_params == 3, f"{ext}: n_free_params lost on round-trip"
            # Summary header (and full text) must be byte-identical post-load.
            assert str(r2.summary()) == pre_summary, f"{ext}: summary changed"


def test_save_load_preserves_walker_axis_shape():
    """samples is (N, n_walkers, ndim); shape must survive the bilby flatten."""
    fitter = _seed_fitter(_make_result())
    assert fitter.result.samples.shape == (100, 4, 3)
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fit.h5")
        fitter.save(path)
        r2 = Fitter.load(path).result
        assert r2.samples.shape == (100, 4, 3)
        assert r2.log_probs.shape == (100, 4)


def test_save_requires_completed_fit():
    """Fitter.save before .fit() raises a clear error -- no silent empty file."""
    f = Fitter(z=0.1, lumi_dist=1e28)
    import pytest
    with pytest.raises(RuntimeError, match="requires a completed fit"):
        f.save("/tmp/should_not_exist.h5")


# ----------------- bilby interop -----------------


def test_saved_file_readable_by_bilby_directly():
    """File written by Fitter.save() should be a valid bilby Result file."""
    import bilby

    fitter = _seed_fitter(_make_result())
    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "fit.h5")
        fitter.save(path)
        br = bilby.read_in_result(filename=path)
        assert "log_likelihood" in br.posterior.columns
        assert tuple(br.search_parameter_keys) == tuple(fitter.result.labels)
        # Our snapshot is preserved under the namespaced metadata key
        assert "vegasafterglow" in br.meta_data


def test_plain_bilby_file_rejected_with_clear_error():
    """A bilby Result file with no VegasAfterglow snapshot can't reconstruct a
    Fitter; Fitter.load should raise a clear ValueError directing the user to
    bilby.read_in_result for inspection-only access."""
    import bilby
    import pandas as pd
    import pytest

    labels = ["x", "y"]
    flat = np.random.RandomState(0).randn(200, 2)
    posterior = pd.DataFrame(flat, columns=labels)
    posterior["log_likelihood"] = np.zeros(200)
    priors = bilby.core.prior.PriorDict(
        {name: bilby.core.prior.Uniform(-10, 10, name) for name in labels}
    )
    br = bilby.core.result.Result(
        label="external",
        outdir=".",
        sampler="dynesty",
        posterior=posterior,
        search_parameter_keys=labels,
        priors=priors,
    )

    with tempfile.TemporaryDirectory() as d:
        path = os.path.join(d, "external.h5")
        br.save_to_file(filename=path)
        with pytest.raises(ValueError, match="does not include a Fitter snapshot"):
            Fitter.load(path)
