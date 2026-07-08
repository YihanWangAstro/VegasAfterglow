"""
Pytest configuration and fixtures for Python API tests.
"""

import base64
import importlib.util
import io
import os
import warnings

import numpy as np
import pytest

_GOLDEN_PLOT_CACHE = {}
_GOLDEN_PLOT_MODULE = None


def _golden_plot_module():
    global _GOLDEN_PLOT_MODULE
    if _GOLDEN_PLOT_MODULE is None:
        plot_path = os.path.join(os.path.dirname(__file__), "golden", "plot.py")
        spec = importlib.util.spec_from_file_location("golden_plot", plot_path)
        _GOLDEN_PLOT_MODULE = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_GOLDEN_PLOT_MODULE)
    return _GOLDEN_PLOT_MODULE


def _golden_case_png(item):
    """PNG bytes comparing a golden case's baseline to the run's recomputation.

    Reuses the recomputed components cached by the test fixture, so no extra
    model evaluation happens; the figure itself is rendered once per case and
    also saved to golden/plots/ for tests/report.py to embed.
    """
    name = item.funcargs.get("name")
    recomputed = item.funcargs.get("recomputed")
    if name is None or recomputed is None:
        return None
    if name not in _GOLDEN_PLOT_CACHE:
        try:
            plots_dir = os.path.join(os.path.dirname(__file__), "golden", "plots")
            os.makedirs(plots_dir, exist_ok=True)
            if not _GOLDEN_PLOT_CACHE:
                # first render of this session: drop leftovers from renamed/removed
                # cases so the report never embeds a figure with no matching test
                for stale in os.listdir(plots_dir):
                    if stale.endswith(".png"):
                        os.remove(os.path.join(plots_dir, stale))
            golden_plot = _golden_plot_module()
            fig, _ = golden_plot.render_case(name, recomputed(name)[0])
            buf = io.BytesIO()
            fig.savefig(buf, format="png", dpi=110)
            golden_plot.plt.close(fig)
            _GOLDEN_PLOT_CACHE[name] = buf.getvalue()
            with open(os.path.join(plots_dir, f"{name}.png"), "wb") as f:
                f.write(_GOLDEN_PLOT_CACHE[name])
        except Exception as exc:  # diagnostics only — never take down the test run
            _GOLDEN_PLOT_CACHE[name] = None
            warnings.warn(f"golden comparison plot for {name!r} not rendered: {exc}",
                          stacklevel=2)
    return _GOLDEN_PLOT_CACHE[name]


@pytest.hookimpl(hookwrapper=True)
def pytest_runtest_makereport(item, call):
    """Render the baseline-vs-run comparison figure for golden tests.

    Active only when a report pipeline consumes it — junitxml (tests/report.py
    via run_all.sh) or pytest-html — so bare `pytest -q` loops pay nothing.
    Each case renders once (on its 'total' component test) plus on any golden
    failure; with pytest-html the figure is also attached to the test entry.
    """
    outcome = yield
    report = outcome.get_result()
    if report.when != "call" or "golden" not in item.keywords:
        return
    html_active = getattr(item.config.option, "htmlpath", None)
    junit_active = getattr(item.config.option, "xmlpath", None)
    if not html_active and not junit_active:
        return
    if not (report.failed or item.funcargs.get("component") == "total"):
        return
    png = _golden_case_png(item)
    if png is None or not html_active:
        return
    try:
        from pytest_html import extras as html_extras
    except ImportError:
        return
    report.extras = getattr(report, "extras", []) + [
        html_extras.png(base64.b64encode(png).decode())]


@pytest.fixture
def sample_time_array():
    """Generate a sample time array in seconds."""
    return np.logspace(2, 6, 30)


@pytest.fixture
def sample_frequency():
    """Return a sample optical frequency in Hz."""
    return 4.84e14  # R-band


@pytest.fixture
def sample_frequency_array():
    """Generate a sample frequency array in Hz."""
    return np.logspace(9, 18, 20)
