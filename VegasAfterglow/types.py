import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Sequence, Tuple

import numpy as np

from .VegasAfterglowC import (  # noqa: F401
    ISM,
    Ejecta,
    GaussianJet,
    Magnetar,
    Medium,
    Model,
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
    logscale_screen,
)


class ModelParams:
    """Complete parameter set for GRB afterglow model configuration.

    All values are in user-facing units (CGS where applicable):
    angles in radians, energies in erg, densities in cm^-3, times in seconds.

    Custom attributes can be added freely for use with custom jet/medium factories
    (e.g., ``params.r_scale = 1.5``).
    """

    def __init__(self):
        self.theta_v = 0.0

        # External medium
        self.n_ism = 0.0
        self.n0 = math.inf
        self.A_star = 0.0
        self.k_m = 2.0

        # Jet core
        self.E_iso = 1e52
        self.Gamma0 = 300.0
        self.theta_c = 0.1
        self.k_e = 2.0
        self.k_g = 2.0
        self.tau = 1.0  # Central engine duration [seconds]

        # Jet wing
        self.E_iso_w = 1e52
        self.Gamma0_w = 300.0
        self.theta_w = math.pi / 2

        # Magnetar
        self.L0 = 0.0
        self.t0 = 1.0
        self.q = 2.0

        # Forward shock radiation
        self.p = 2.3
        self.eps_e = 0.1
        self.eps_B = 0.01
        self.xi_e = 1.0

        # Reverse shock radiation
        self.p_r = 2.3
        self.eps_e_r = 0.1
        self.eps_B_r = 0.01
        self.xi_e_r = 1.0

        # Host-galaxy extinction
        self.A_V = 0.0

    def __repr__(self) -> str:
        """Show only fields that differ from the constructor defaults, plus any
        custom attributes the user attached. Empty namespace -> "ModelParams(defaults)".
        """
        defaults = vars(ModelParams())
        items = [
            f"{name}={val!r}"
            for name, val in vars(self).items()
            if name not in defaults or defaults[name] != val
        ]
        return f"ModelParams({', '.join(items)})" if items else "ModelParams(defaults)"


class _SummaryTable:
    """Wraps a formatted multi-line string so it renders cleanly in both
    ``print(...)`` (via ``__str__``) and Jupyter last-line auto-display
    (via ``__repr__``). Test coverage: ``tests/python/test_fit_result_summary.py``
    pins this behavior — do not collapse to a plain ``str``."""

    __slots__ = ("_text",)

    def __init__(self, text: str) -> None:
        self._text = text

    def __repr__(self) -> str:
        return self._text

    __str__ = __repr__


def _auto_format(values: np.ndarray) -> str:
    """Pick a printf-style format string for a column of values, targeting
    ~4 significant figures. Falls back to ``.4g`` (scientific) for values
    outside [1e-2, 1e5)."""
    arr = np.asarray(values, dtype=float).ravel()
    if arr.size == 0:
        return "{:.4f}"
    vmax = float(np.max(np.abs(arr)))
    if vmax == 0 or not np.isfinite(vmax):
        return "{:.4f}"
    magnitude = int(np.floor(np.log10(vmax)))
    if magnitude < -2 or magnitude > 4:
        return "{:.4g}"
    decimals = max(0, 3 - magnitude)
    return f"{{:.{decimals}f}}"


def _strip_dollars(s: str) -> str:
    """Strip a single pair of outer ``$`` delimiters if present, so the label
    can be embedded inside a larger LaTeX equation."""
    s = s.strip()
    if len(s) >= 2 and s.startswith("$") and s.endswith("$"):
        return s[1:-1]
    return s


def _format_value_with_bounds(
    value: float, upper: float, lower: float
) -> Tuple[str, str, str]:
    """Format ``(value, +upper, -lower)`` using the physics-paper convention:
    2 significant figures in the smaller positive bound, with the central
    value rounded to the same decimal place."""
    candidates = [abs(b) for b in (upper, lower) if b > 0 and np.isfinite(b)]
    if not candidates:
        # No usable bound; fall back to 2 sig figs on the value itself.
        ref = abs(value) if value != 0 else 1.0
    else:
        ref = min(candidates)
    if ref == 0 or not np.isfinite(ref):
        fmt = "{:.2f}"
    else:
        magnitude = int(np.floor(np.log10(ref)))
        decimals = max(0, 1 - magnitude)
        fmt = f"{{:.{decimals}f}}"
    return fmt.format(value), fmt.format(upper), fmt.format(lower)


@dataclass(repr=False)
class FitResult:
    """
    The result of a sampling fit.
    """

    samples: np.ndarray
    log_probs: np.ndarray
    labels: Sequence[str]
    latex_labels: Optional[Sequence[str]] = None  # LaTeX labels for plotting
    top_k_params: Optional[np.ndarray] = None
    top_k_log_probs: Optional[np.ndarray] = None
    bilby_result: Optional[object] = None  # Full bilby Result object (for diagnostics)
    # Fit-quality bookkeeping (None for results loaded from older files).
    n_data: Optional[int] = None  # Total observation count across point + band data
    n_free_params: Optional[int] = None  # Sampled (non-fixed) parameter count

    def __repr__(self) -> str:
        """One-line summary; deliberately avoids dumping the (potentially 100k-row)
        ``samples`` / ``log_probs`` arrays the dataclass default would print."""
        n_samples = self.samples.shape[0]
        ndim = self.samples.shape[-1]
        lp_best = float(self.log_probs.max()) if self.log_probs.size else float("nan")
        n_top = len(self.top_k_params) if self.top_k_params is not None else 0
        return (
            f"FitResult(n_samples={n_samples}, ndim={ndim}, "
            f"labels={list(self.labels)}, top_k={n_top}, best_logL={lp_best:.4g})"
        )

    @property
    def flat_samples(self):
        """Posterior samples flattened to shape (n_steps * n_walkers, ndim)."""
        return self.samples.reshape(-1, self.samples.shape[-1])

    def summary(
        self,
        top_k: Optional[int] = None,
        latex: Optional[bool] = None,
    ) -> _SummaryTable:
        """Top-K best-fit table, ranked by log-likelihood.

        Columns: ``Rank``, ``chi^2``, and the sampler-space value of each
        parameter (parameters defined with ``Scale.log`` appear with the
        ``log10_`` prefix matching ``self.labels``).

        The ``chi^2`` column is ``-2 * log_prob`` -- exact for the default
        Gaussian likelihood (``log_likelihood_fn=None``), a generic deviance
        for custom likelihoods.

        When ``latex_labels`` are set on the ``FitResult``, a single
        copy-pasteable LaTeX block is appended showing each parameter's
        posterior median with asymmetric 1σ bounds (16th/50th/84th
        percentiles of ``self.samples``). These values match the conventional
        ``corner`` plot annotations.

        Parameters
        ----------
        top_k : int, optional
            Limit output to the first ``top_k`` rows. If ``None`` (default),
            shows all rows stored in ``self.top_k_params``.
        latex : bool, optional
            Force the LaTeX block on/off. If ``None`` (default), the block is
            included when ``self.latex_labels`` is set.

        Returns
        -------
        _SummaryTable
            Formatted table; renders identically under ``print()`` and
            Jupyter last-line evaluation.
        """
        if (
            self.top_k_params is None
            or self.top_k_log_probs is None
            or len(self.top_k_params) == 0
        ):
            return _SummaryTable("FitResult: no top_k_params stored")

        rows = np.asarray(
            self.top_k_params if top_k is None else self.top_k_params[:top_k],
            dtype=float,
        )
        logp = np.asarray(self.top_k_log_probs[: len(rows)], dtype=float)
        chi2_vals = -2.0 * logp
        n_rows = len(rows)
        n_params = rows.shape[1]
        m_total = len(self.top_k_params)

        chi2_fmt = _auto_format(chi2_vals)
        param_fmts = [_auto_format(rows[:, j]) for j in range(n_params)]

        rank_strs = [str(i + 1) for i in range(n_rows)]
        chi2_strs = [chi2_fmt.format(v) for v in chi2_vals]
        param_strs = [
            [param_fmts[j].format(rows[i, j]) for j in range(n_params)]
            for i in range(n_rows)
        ]

        rank_w = max(len("Rank"), max(len(s) for s in rank_strs))
        chi2_w = max(len("chi^2"), max(len(s) for s in chi2_strs))
        param_ws = [
            max(
                len(self.labels[j]),
                max(len(param_strs[i][j]) for i in range(n_rows)),
            )
            for j in range(n_params)
        ]

        gap = "   "

        def _row(rank_cell: str, chi2_cell: str, param_cells: Sequence[str]) -> str:
            parts = [f"{rank_cell:>{rank_w}}", f"{chi2_cell:>{chi2_w}}"]
            parts.extend(f"{param_cells[j]:>{param_ws[j]}}" for j in range(n_params))
            return gap.join(parts)

        title = f"Best-fit summary (top {n_rows} of {m_total})"
        lines = [title, "=" * len(title)]

        # Suppress the fit-quality header when n_data / n_free_params is
        # missing (older saved fits) or DOF would be <= 0 (over-parametrized).
        if (
            self.n_data is not None
            and self.n_free_params is not None
            and self.n_data > self.n_free_params
        ):
            chi2_min = float(-2.0 * np.max(self.top_k_log_probs))
            n, k = self.n_data, self.n_free_params
            dof = n - k
            bic = chi2_min + k * float(np.log(n))
            aic = chi2_min + 2 * k
            lines.append(
                f"χ²_min / DOF = {chi2_min:.2f} / {dof} = {chi2_min / dof:.3g}"
                f"    BIC = {bic:.3g}    AIC = {aic:.3g}"
            )
            lines.append("")

        lines.append(_row("Rank", "chi^2", list(self.labels)))
        lines.append(
            _row(
                "-" * rank_w,
                "-" * chi2_w,
                ["-" * w for w in param_ws],
            )
        )
        for i in range(n_rows):
            lines.append(_row(rank_strs[i], chi2_strs[i], param_strs[i]))

        show_latex = latex if latex is not None else (self.latex_labels is not None)
        if show_latex:
            if self.latex_labels is None:
                lines.append("")
                lines.append(
                    "# (no latex_labels set; pass latex_labels=[...] to FitResult"
                    " to enable LaTeX block)"
                )
            else:
                flat = self.flat_samples
                p16, p50, p84 = np.percentile(flat, [16, 50, 84], axis=0)
                lines.append("")
                lines.append("LaTeX (median ± 1σ, matches corner plot):")
                for j in range(n_params):
                    label = _strip_dollars(self.latex_labels[j])
                    med_s, up_s, lo_s = _format_value_with_bounds(
                        float(p50[j]),
                        float(p84[j] - p50[j]),
                        float(p50[j] - p16[j]),
                    )
                    lines.append(f"  ${label} = {med_s}^{{+{up_s}}}_{{-{lo_s}}}$")

        return _SummaryTable("\n".join(lines))


class Scale(Enum):
    linear = "linear"
    log = "log"
    fixed = "fixed"

    # Uppercase aliases — supported API; tests/python/test_types.py pins this.
    # Both forms refer to the same enum members (Python Enum value-collision aliasing).
    LINEAR = "linear"
    LOG = "log"
    FIXED = "fixed"

    def __repr__(self) -> str:
        return f"Scale.{self.name}"

    @classmethod
    def _missing_(cls, value) -> Optional["Scale"]:
        if isinstance(value, str):
            for member in cls:
                if member.value == value.lower():
                    return member
        return None


@dataclass
class ParamDef:
    """
    Single-parameter definition for MCMC.
    scale=LOG means we sample log10(x), then transform via 10**v.
    scale=FIXED means this param never appears in the sampler.
    """

    name: str
    lower: float
    upper: float
    scale: Scale = Scale.linear
    initial: Optional[
        float
    ] = None  # Initial value (in linear space, auto-converted for LOG scale)

    def __post_init__(self) -> None:
        """Catch obvious mistakes at construction time so the user sees the
        offending parameter immediately, not 30 seconds into the first MCMC
        step."""
        # Scale.fixed collapses bounds to a single value; ``lower`` serves as
        # the value when ``initial`` is omitted (see fitting/utils.py), so we
        # do not require initial here.
        if self.scale is Scale.fixed:
            return
        if self.lower >= self.upper:
            raise ValueError(
                f"ParamDef(name={self.name!r}): lower={self.lower} must be "
                f"strictly less than upper={self.upper}."
            )
        if self.scale is Scale.log and self.lower <= 0:
            raise ValueError(
                f"ParamDef(name={self.name!r}): scale=Scale.log requires "
                f"lower > 0 (log10 is undefined for non-positive values); "
                f"got lower={self.lower}."
            )
        if self.initial is not None and not (self.lower <= self.initial <= self.upper):
            raise ValueError(
                f"ParamDef(name={self.name!r}): initial={self.initial} must lie "
                f"in [lower={self.lower}, upper={self.upper}]."
            )
