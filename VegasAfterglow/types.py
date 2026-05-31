import math
from dataclasses import dataclass
from enum import Enum
from typing import Optional, Sequence

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
    (via ``__repr__``)."""

    __slots__ = ("_text",)

    def __init__(self, text: str) -> None:
        self._text = text

    def __repr__(self) -> str:
        return self._text

    __str__ = __repr__


@dataclass(repr=False)
class FitResult:
    """
    The result of a sampling fit.
    """

    samples: np.ndarray
    log_probs: np.ndarray
    labels: Sequence[str]
    latex_labels: Sequence[str] = None  # LaTeX labels for plotting
    top_k_params: np.ndarray = None
    top_k_log_probs: np.ndarray = None
    bilby_result: object = None  # Full bilby Result object (for diagnostics)

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

    def summary(self, top_k: Optional[int] = None) -> _SummaryTable:
        """Top-K best-fit table, ranked by log-likelihood.

        Columns: ``Rank``, ``chi^2``, and the sampler-space value of each
        parameter (parameters defined with ``Scale.log`` appear with the
        ``log10_`` prefix matching ``self.labels``).

        The ``chi^2`` column is ``-2 * log_prob`` -- exact for the default
        Gaussian likelihood (``log_likelihood_fn=None``), a generic deviance
        for custom likelihoods.

        Parameters
        ----------
        top_k : int, optional
            Limit output to the first ``top_k`` rows. If ``None`` (default),
            shows all rows stored in ``self.top_k_params``.

        Returns
        -------
        _SummaryTable
            Formatted table; renders identically under ``print()`` and
            Jupyter last-line evaluation.
        """
        if self.top_k_params is None or len(self.top_k_params) == 0:
            return _SummaryTable("FitResult: no top_k_params stored")
        rows = self.top_k_params if top_k is None else self.top_k_params[:top_k]
        logp = self.top_k_log_probs[: len(rows)]
        name_w = max(max(len(s) for s in self.labels), 10)
        header = f"{'Rank':>4}  {'chi^2':>10}  " + "  ".join(
            f"{name:>{name_w}}" for name in self.labels
        )
        lines = [header, "-" * len(header)]
        for i, params in enumerate(rows):
            chi2 = -2.0 * logp[i]
            vals = "  ".join(f"{v:>{name_w}.4f}" for v in params)
            lines.append(f"{i + 1:>4d}  {chi2:>10.2f}  {vals}")
        return _SummaryTable("\n".join(lines))


class Scale(Enum):
    linear = "linear"
    log = "log"
    fixed = "fixed"

    # Uppercase aliases for backward compatibility
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
    scale: Scale = Scale.LINEAR
    initial: float = (
        None  # Initial value (in linear space, auto-converted for LOG scale)
    )
