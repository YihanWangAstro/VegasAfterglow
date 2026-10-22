"""Afterglow model fitting with custom priors, likelihoods, and jet/medium profiles."""

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from types import SimpleNamespace
from typing import Callable, List, Literal, Optional, Sequence, Tuple

import numpy as np
from scipy.special import log_ndtr

from ..types import FitResult, ModelParams, ParamDef, Scale
from ..VegasAfterglowC import Model, Observer, Radiation
from .config import _BandLimits, _BandObs
from .utils import _build_transformer, _default_jet_factory, _default_medium_factory

logger = logging.getLogger(__name__)
# Library convention: attach a NullHandler so the package stays silent unless the
# application configures logging itself. See
# https://docs.python.org/3/howto/logging.html#configuring-logging-for-a-library
if not logger.handlers:
    logger.addHandler(logging.NullHandler())

_C_CGS = 2.99792458e10  # speed of light [cm/s]
_LN10_OVER_2P5 = 0.4 * np.log(10.0)  # 0.4 * ln(10), converts A_V * k -> optical depth


class Fitter:
    """Model-based afterglow fitter with custom priors, likelihood, and jet/medium profiles.

    Uses ``Model`` directly for flux evaluation, with ``ThreadPoolExecutor``
    parallelism (each thread releases the GIL in the C++ computation).

    Resolution lifecycle
    --------------------
    The constructor's ``resolution`` is the session default. ``fit(resolution=…)``
    overrides it *and persists* so post-fit calls (``flux_density_grid``,
    ``flux``) inherit it. Those post-fit methods accept their own one-shot
    ``resolution=`` override that does not mutate the Fitter.

    Likelihood
    ----------
    The chi-squared is evaluated in log-flux space:
    ``chi2 = sum(((ln F_obs - ln F_model) / sigma_ln)**2)`` with
    ``sigma_ln = err / F_obs`` propagated from the linear errors passed to the
    ``add_*`` methods. Scale-invariant, so bright points do not dominate data
    spanning decades; fluxes and errors must be strictly positive.

    Each ``add_*`` call is one dataset. A dataset given with ``sys_err`` shares
    one multiplicative calibration factor with a Gaussian prior of that width
    in ln-flux, marginalised analytically: its shape keeps the full precision
    of the statistical errors while its level may shift by about ``sys_err``.
    Upper limits added with the ``add_*_limits`` methods contribute
    ``-2 ln Phi((U - F_model) / sigma)`` in linear flux. ``log_likelihood_fn``
    (see :meth:`fit`) receives the total of all these terms as one chi-squared.

    Thread safety
    -------------
    Set all state (``add_*``, per-call ``resolution=``) before calling ``fit()``;
    do not mutate the Fitter from another thread while a fit is running.

    Args:
        z: Source redshift
        lumi_dist: Luminosity distance [cm]
        jet: Jet type string ("tophat", "gaussian", "powerlaw", "two_component",
            "step_powerlaw", "powerlaw_wing", "uniform") or a callable
            factory(ModelParams) -> jet object for custom profiles
        medium: Medium type string ("ism", "wind") or a callable
            factory(ModelParams) -> medium object for custom profiles
        fwd_ssc: Enable forward shock inverse Compton
        rvs_ssc: Enable reverse shock inverse Compton
        rvs_shock: Enable reverse shock
        kn: Enable Klein-Nishina corrections
        magnetar: Enable magnetar energy injection
        radiative_fireball: Radiative losses decelerate the blast wave (default True).
            Set False for the adiabatic approximation used by most afterglow codes.
        rtol: Numerical tolerance
        resolution: Grid resolution tuple (phi, theta, t). Defaults to (0.06, 0.15, 5),
            or (0.06, 0.2, 10) when the reverse shock is enabled. Can be overridden per-fit.
        extinction: Host-galaxy dust extinction. ``None`` (default) disables it;
            a string ``"smc" | "lmc" | "mw"`` selects a built-in Pei (1992) law;
            a callable ``f(lam_cm, params) -> k(lam)`` supplies a custom law.
            Requires a fitted ``A_V`` parameter via ``ParamDef``. Galactic
            (Milky Way) extinction along the line of sight should be removed
            from the data before fitting; this layer handles host-galaxy only.
    """

    #: Boolean physics/config flags with their defaults — the single source of truth
    #: for ``__repr__`` and the save/load snapshot in ``fitting/io.py``. Adding a flag:
    #: one entry here plus the constructor keyword.
    CONFIG_FLAGS = {
        "fwd_ssc": False,
        "rvs_ssc": False,
        "rvs_shock": False,
        "kn": False,
        "magnetar": False,
        "radiative_fireball": True,
    }

    # ── Construction & configuration ─────────────────────────────────────────────────────────

    def __init__(
        self,
        *,
        z: float = 0.0,
        lumi_dist: float = 1e26,
        jet="tophat",
        medium="ism",
        fwd_ssc: bool = False,
        rvs_ssc: bool = False,
        rvs_shock: bool = False,
        kn: bool = False,
        magnetar: bool = False,
        radiative_fireball: bool = True,
        rtol: float = 1e-6,
        resolution: Optional[Tuple[float, float, float]] = None,
        extinction=None,
    ):
        self.z = z
        self.lumi_dist = lumi_dist
        self.fwd_ssc = fwd_ssc
        self.rvs_ssc = rvs_ssc
        self.rvs_shock = rvs_shock
        self.kn = kn
        self.magnetar = magnetar
        self.radiative_fireball = radiative_fireball
        self.rtol = rtol
        # None forwards to Model, which resolves the mode-aware default
        # (defaults::grid in the C++ core is the single source of truth).
        self.resolution = tuple(resolution) if resolution is not None else None

        self.jet = jet
        self._custom_jet = callable(jet)
        self.jet_factory = jet if self._custom_jet else _default_jet_factory(self)

        self.medium = medium
        self._custom_medium = callable(medium)
        self.medium_factory = (
            medium if self._custom_medium else _default_medium_factory(self)
        )

        self.extinction = extinction
        self._custom_extinction = callable(extinction)
        self._ext_law = self._resolve_extinction(extinction)
        # Populated lazily by _consolidate_data() once observation data is known.
        self._ext_kernel = None  # 0.4*ln(10) * k(lam_rest), built-in path
        self._lam_rest_cm = None  # rest-frame wavelengths, custom path

        # Set by Fitter.fit() / Fitter.load(); persisted by Fitter.save().
        self._to_params = None
        self._param_defs: Optional[List[ParamDef]] = None
        self.result: Optional[FitResult] = None

        # Observation data (CGS units throughout).
        self._point_t = []
        self._point_nu = []
        self._point_flux = []
        self._point_err = []
        self._point_weights = []
        # Per-entry optional filter/band name for draw_fit legends.
        self._point_labels: List[Optional[str]] = []
        # Per-entry fractional calibration uncertainty (0 = none).
        self._point_sys: List[float] = []
        self._band_obs: List[_BandObs] = []
        # Upper limits (non-detections): limit value and the noise level behind it.
        self._limit_t = []
        self._limit_nu = []
        self._limit_value = []
        self._limit_sigma = []
        self._limit_labels: List[Optional[str]] = []
        self._band_limits: List[_BandLimits] = []

        # Consolidated arrays built lazily by _consolidate_data().
        self._all_t = None
        self._all_nu = None
        self._all_flux = None
        self._all_err = None
        self._all_weights = None
        self._all_log_flux = None
        self._all_log_err = None
        # Calibration-systematics bookkeeping: dataset id per point, per-dataset
        # inverse prior variance (inf when sys_err = 0) and total precision.
        self._all_prec = None
        self._all_dataset = None
        self._dataset_inv_sys2 = None
        self._dataset_prec_sum = None
        # Model-evaluation arrays: detections and limits merged, sorted by time,
        # with the position of each detection / limit in that order.
        self._eval_t = None
        self._eval_nu = None
        self._det_pos = None
        self._lim_pos = None
        self._lim_t = None
        self._lim_value = None
        self._lim_sigma = None

    def __repr__(self) -> str:
        """One-line summary of model selection, physics flags, and loaded data."""
        flags = [
            name if getattr(self, name) else f"{name}=False"
            for name, default in self.CONFIG_FLAGS.items()
            if getattr(self, name) != default
        ]
        flag_str = f", flags=[{', '.join(flags)}]" if flags else ""
        ext_str = (
            f", extinction={self.extinction!r}" if self.extinction is not None else ""
        )
        n_pts = sum(len(t) for t in self._point_t)
        if n_pts or self._band_obs:
            data_str = f", data={n_pts} points + {len(self._band_obs)} bands"
        else:
            data_str = ", data=empty"
        return (
            f"Fitter(z={self.z}, lumi_dist={self.lumi_dist:.3g}, "
            f"jet={self.jet!r}, medium={self.medium!r}"
            f"{ext_str}{flag_str}{data_str})"
        )

    # ── Observation data input ─────────────────────────────────────────────────────────

    def _add_point_data(self, t, nu, f_nu, err, weights, label=None, sys_err=0.0):
        """Append point observation arrays to internal lists.

        ``label`` is an optional filter/band name used by ``draw_fit`` for
        legend display (e.g., ``'r'``, ``'WXT'``, ``'VT_R'``). Not used for
        chi-squared evaluation. ``sys_err`` is the dataset's fractional
        calibration uncertainty (see :meth:`add_flux_density`).
        """
        f_nu = np.asarray(f_nu, dtype=np.float64)
        err = np.asarray(err, dtype=np.float64)
        w = (
            np.asarray(weights, dtype=np.float64)
            if weights is not None
            else np.ones_like(f_nu)
        )
        self._point_t.append(t)
        self._point_nu.append(nu)
        self._point_flux.append(f_nu)
        self._point_err.append(err)
        self._point_weights.append(w)
        self._point_labels.append(label)
        self._point_sys.append(sys_err)
        self._all_t = None

    @staticmethod
    def _validate_sys_err(sys_err, context: str) -> float:
        sys_err = float(sys_err)
        if not np.isfinite(sys_err) or sys_err < 0:
            raise ValueError(
                f"{context}: sys_err must be finite and >= 0 (fractional calibration "
                f"uncertainty of the dataset), got {sys_err}"
            )
        return sys_err

    @staticmethod
    def _validate_observation_arrays(*, t, f_nu, err, weights, context: str):
        """Validate ``t``, ``f_nu``, ``err`` (and optional ``weights``) at the
        boundary of an ``add_*`` method. Catches shape mismatch, empty data,
        non-finite values, and non-positive error bars with a clear ValueError
        naming the caller (``context``).

        Returns the arrays converted to ``float64`` numpy arrays so callers
        don't have to re-cast them.
        """
        t = np.asarray(t, dtype=np.float64)
        f_nu = np.asarray(f_nu, dtype=np.float64)
        err = np.asarray(err, dtype=np.float64)
        if t.size == 0:
            raise ValueError(f"{context}: time array is empty")
        if t.shape != f_nu.shape or t.shape != err.shape:
            raise ValueError(
                f"{context}: t, f_nu, err must have the same shape; got "
                f"t.shape={t.shape}, f_nu.shape={f_nu.shape}, err.shape={err.shape}"
            )
        if not np.isfinite(f_nu).all():
            raise ValueError(
                f"{context}: f_nu contains {int((~np.isfinite(f_nu)).sum())} "
                f"non-finite (NaN or inf) values"
            )
        if not np.isfinite(err).all() or (err <= 0).any():
            raise ValueError(
                f"{context}: err must be finite and > 0 at every point (got "
                f"min={float(err.min())}, max={float(err.max())})"
            )
        w = None
        if weights is not None:
            w = np.asarray(weights, dtype=np.float64)
            if w.shape != t.shape:
                raise ValueError(
                    f"{context}: weights.shape={w.shape} must match t.shape={t.shape}"
                )
            if not np.isfinite(w).all() or (w < 0).any():
                raise ValueError(
                    f"{context}: weights must be finite and >= 0 at every point"
                )
        return t, f_nu, err, w

    def add_flux_density(
        self,
        nu: float,
        t: np.ndarray,
        f_nu: np.ndarray,
        err: np.ndarray,
        weights: Optional[np.ndarray] = None,
        label: Optional[str] = None,
        sys_err: float = 0.0,
    ) -> None:
        """Add light curve data at a single frequency.

        Args:
            nu: Observing frequency [Hz]
            t: Observation times [seconds]
            f_nu: Observed flux densities [erg/cm^2/s/Hz]
            err: Observational (statistical) uncertainties [erg/cm^2/s/Hz]
            weights: Optional statistical weights
            label: Optional filter/band name (e.g., ``'r'``, ``'VT_R'``, ``'WXT'``)
                used by ``draw_fit`` for legend display. Not used for
                chi-squared evaluation. If omitted, the broad-band name
                (Radio / IR / Optical / X-ray / ...) is used as a fallback.
            sys_err: Fractional calibration uncertainty of this dataset
                (e.g. ``0.1`` for a 10% flux-scale uncertainty). All points of
                the call share one multiplicative calibration factor with a
                Gaussian prior of this width in ln-flux; the likelihood
                marginalises it analytically, so the dataset's shape keeps its
                full precision while its overall level may shift by about
                ``sys_err``. Default ``0`` (no calibration freedom). The
                best-fit shift is reported by :meth:`calibration_offsets`.
        """
        if not np.isfinite(nu) or nu <= 0:
            raise ValueError(f"add_flux_density: nu must be finite and > 0, got {nu}")
        sys_err = self._validate_sys_err(sys_err, "add_flux_density")
        t, f_nu, err, weights = self._validate_observation_arrays(
            t=t, f_nu=f_nu, err=err, weights=weights, context="add_flux_density"
        )
        self._add_point_data(
            t, np.full_like(t, nu), f_nu, err, weights, label=label, sys_err=sys_err
        )

    def add_spectrum(
        self,
        t: float,
        nu: np.ndarray,
        f_nu: np.ndarray,
        err: np.ndarray,
        weights: Optional[np.ndarray] = None,
        label: Optional[str] = None,
        sys_err: float = 0.0,
    ) -> None:
        """Add a broadband spectrum at a specific time.

        Args:
            t: Observation time [seconds]
            nu: Observing frequencies [Hz]
            f_nu: Observed flux densities [erg/cm^2/s/Hz]
            err: Observational (statistical) uncertainties [erg/cm^2/s/Hz]
            weights: Optional statistical weights
            label: Optional instrument name used in reports.
            sys_err: Fractional calibration uncertainty shared by the whole
                spectrum (one instrument per call); see :meth:`add_flux_density`.
        """
        if not np.isfinite(t) or t <= 0:
            raise ValueError(f"add_spectrum: t must be finite and > 0, got {t}")
        sys_err = self._validate_sys_err(sys_err, "add_spectrum")
        nu = np.asarray(nu, dtype=np.float64)
        if not np.isfinite(nu).all() or (nu <= 0).any():
            raise ValueError(
                "add_spectrum: nu must be finite and > 0 at every point "
                f"(got min={float(nu.min())}, max={float(nu.max())})"
            )
        # nu plays the role of the "axis" array here (one row per frequency at
        # a fixed time), so validate shapes against it.
        nu, f_nu, err, weights = self._validate_observation_arrays(
            t=nu, f_nu=f_nu, err=err, weights=weights, context="add_spectrum"
        )
        self._add_point_data(
            np.full_like(nu, t), nu, f_nu, err, weights, label=label, sys_err=sys_err
        )

    def add_flux(
        self,
        band,
        t: np.ndarray,
        flux: np.ndarray,
        err: np.ndarray,
        num_points: int = 5,
        weights: Optional[np.ndarray] = None,
        label: Optional[str] = None,
        sys_err: float = 0.0,
    ) -> None:
        """Add band-integrated flux measurements.

        Args:
            band: Frequency band as (nu_min, nu_max) tuple in Hz.
                Use ``units.band("XRT")`` for named instrument bands.
            t: Observation times [seconds]
            flux: Observed integrated fluxes [erg/cm^2/s]
            err: Observational uncertainties [erg/cm^2/s]
            num_points: Frequency sampling points for the band integration
                (Boole's rule on a log grid; best with num_points = 4k+1 —
                5 covers optical/X-ray instrument bands, use 9-13 for radio
                or bands wider than two decades, 17+ for bolometric)
            weights: Optional statistical weights
            label: Optional instrument/band name (e.g., ``'WXT'``, ``'XRT'``)
                used by ``draw_fit`` for legend display.
            sys_err: Fractional calibration uncertainty shared by this
                dataset; see :meth:`add_flux_density`.
        """
        sys_err = self._validate_sys_err(sys_err, "add_flux")
        try:
            nu_min, nu_max = band
        except (TypeError, ValueError):
            raise ValueError(
                f"add_flux: band must be a (nu_min, nu_max) tuple in Hz, got {band!r}"
            ) from None
        if not (np.isfinite(nu_min) and np.isfinite(nu_max) and 0 < nu_min < nu_max):
            raise ValueError(
                f"add_flux: band must satisfy 0 < nu_min < nu_max with both finite; "
                f"got nu_min={nu_min}, nu_max={nu_max}"
            )
        if num_points < 2:
            raise ValueError(
                f"add_flux: num_points must be >= 2 for band integration, got {num_points}"
            )
        t, flux, err, w = self._validate_observation_arrays(
            t=t, f_nu=flux, err=err, weights=weights, context="add_flux"
        )
        if w is None:
            w = np.ones_like(t)

        order = np.argsort(t)
        self._band_obs.append(
            _BandObs(
                nu_min=nu_min,
                nu_max=nu_max,
                num_points=num_points,
                t=t[order],
                flux=flux[order],
                err=err[order],
                weights=w[order],
                name=label,
                sys_err=sys_err,
            )
        )

    # ── Upper limits ─────────────────────────────────────────────────────────

    @staticmethod
    def _validate_axis(axis, name: str, context: str) -> np.ndarray:
        axis = np.atleast_1d(np.asarray(axis, dtype=np.float64))
        if not np.isfinite(axis).all() or (axis <= 0).any():
            raise ValueError(f"{context}: {name} must be finite and > 0 at every point")
        return axis

    @staticmethod
    def _validate_limit_arrays(*, axis, limit, sigma, n_sigma, context: str):
        limit = np.atleast_1d(np.asarray(limit, dtype=np.float64))
        if limit.shape != axis.shape:
            raise ValueError(
                f"{context}: limit must have the same shape as the time/frequency "
                f"array {axis.shape}, got {limit.shape}"
            )
        if not np.isfinite(limit).all() or (limit <= 0).any():
            raise ValueError(f"{context}: limit values must be finite and > 0")
        if sigma is None:
            n_sigma = float(n_sigma)
            if not np.isfinite(n_sigma) or n_sigma <= 0:
                raise ValueError(
                    f"{context}: n_sigma must be finite and > 0, got {n_sigma}"
                )
            sigma = limit / n_sigma
        else:
            sigma = np.atleast_1d(np.asarray(sigma, dtype=np.float64))
            if sigma.shape != limit.shape:
                raise ValueError(
                    f"{context}: sigma must have the same shape as limit {limit.shape}, "
                    f"got {sigma.shape}"
                )
            if not np.isfinite(sigma).all() or (sigma <= 0).any():
                raise ValueError(f"{context}: sigma values must be finite and > 0")
        return limit, sigma

    def _add_limit_data(self, t, nu, limit, sigma, label):
        self._limit_t.append(t)
        self._limit_nu.append(nu)
        self._limit_value.append(limit)
        self._limit_sigma.append(sigma)
        self._limit_labels.append(label)
        self._all_t = None

    def add_flux_density_limits(
        self,
        nu: float,
        t: np.ndarray,
        limit: np.ndarray,
        n_sigma: float = 3.0,
        sigma: Optional[np.ndarray] = None,
        label: Optional[str] = None,
    ) -> None:
        """Add flux-density upper limits (non-detections) at a single frequency.

        A limit ``U`` quoted at ``n_sigma`` has noise ``sigma = U / n_sigma``
        behind it. Its likelihood is the probability that a measurement would
        have come out below ``U`` given the model, ``Phi((U - F_model) / sigma)``:
        no penalty while the model sits well below the limit, ``2 ln 2`` in
        chi2 units when it touches it, and ``((F_model - U) / sigma)^2`` beyond.

        Args:
            nu: Observing frequency [Hz]
            t: Observation times [seconds]
            limit: Upper-limit flux densities [erg/cm^2/s/Hz]
            n_sigma: Significance level the limits are quoted at (default 3).
                Ignored when ``sigma`` is given.
            sigma: Optional noise level behind each limit [erg/cm^2/s/Hz]
                (same shape as ``limit``); overrides ``n_sigma``.
            label: Optional filter/band name for reports and ``draw_fit``.
        """
        if not np.isfinite(nu) or nu <= 0:
            raise ValueError(
                f"add_flux_density_limits: nu must be finite and > 0, got {nu}"
            )
        t = self._validate_axis(t, "t", "add_flux_density_limits")
        limit, sigma = self._validate_limit_arrays(
            axis=t,
            limit=limit,
            sigma=sigma,
            n_sigma=n_sigma,
            context="add_flux_density_limits",
        )
        self._add_limit_data(t, np.full_like(t, nu), limit, sigma, label)

    def add_spectrum_limits(
        self,
        t: float,
        nu: np.ndarray,
        limit: np.ndarray,
        n_sigma: float = 3.0,
        sigma: Optional[np.ndarray] = None,
        label: Optional[str] = None,
    ) -> None:
        """Add flux-density upper limits at one time across frequencies.

        See :meth:`add_flux_density_limits` for the likelihood and arguments.
        """
        if not np.isfinite(t) or t <= 0:
            raise ValueError(f"add_spectrum_limits: t must be finite and > 0, got {t}")
        nu = self._validate_axis(nu, "nu", "add_spectrum_limits")
        limit, sigma = self._validate_limit_arrays(
            axis=nu,
            limit=limit,
            sigma=sigma,
            n_sigma=n_sigma,
            context="add_spectrum_limits",
        )
        self._add_limit_data(np.full_like(nu, t), nu, limit, sigma, label)

    def add_flux_limits(
        self,
        band,
        t: np.ndarray,
        limit: np.ndarray,
        n_sigma: float = 3.0,
        sigma: Optional[np.ndarray] = None,
        num_points: int = 5,
        label: Optional[str] = None,
    ) -> None:
        """Add band-integrated flux upper limits [erg/cm^2/s].

        ``band`` and ``num_points`` are as in :meth:`add_flux`; the limit
        likelihood is as in :meth:`add_flux_density_limits`.
        """
        try:
            nu_min, nu_max = band
        except (TypeError, ValueError):
            raise ValueError(
                f"add_flux_limits: band must be a (nu_min, nu_max) tuple in Hz, got {band!r}"
            ) from None
        if not (np.isfinite(nu_min) and np.isfinite(nu_max) and 0 < nu_min < nu_max):
            raise ValueError(
                f"add_flux_limits: band must satisfy 0 < nu_min < nu_max with both finite; "
                f"got nu_min={nu_min}, nu_max={nu_max}"
            )
        if num_points < 2:
            raise ValueError(
                f"add_flux_limits: num_points must be >= 2 for band integration, got {num_points}"
            )
        t = self._validate_axis(t, "t", "add_flux_limits")
        limit, sigma = self._validate_limit_arrays(
            axis=t, limit=limit, sigma=sigma, n_sigma=n_sigma, context="add_flux_limits"
        )
        order = np.argsort(t)
        self._band_limits.append(
            _BandLimits(
                nu_min=nu_min,
                nu_max=nu_max,
                num_points=num_points,
                t=t[order],
                limit=limit[order],
                sigma=sigma[order],
                name=label,
            )
        )

    def _data_counts(self) -> Tuple[int, int]:
        """``(n_data, n_limits)``: all observations including limits, and limits alone."""
        n_det = sum(len(a) for a in self._point_t) + sum(
            len(b.t) for b in self._band_obs
        )
        n_lim = sum(len(a) for a in self._limit_t) + sum(
            len(b.t) for b in self._band_limits
        )
        return n_det + n_lim, n_lim

    # ── Extinction & data consolidation ─────────────────────────────────────────────────────────

    def _resolve_extinction(self, ext):
        """Resolve the ``extinction`` constructor argument to a callable or None.

        Built-in string profiles are wrapped to discard the unused ``params``
        argument so the call site is uniform regardless of source.
        """
        if ext is None:
            return None
        if callable(ext):
            return ext
        from ..extinction import BUILTIN_LAWS

        if ext not in BUILTIN_LAWS:
            raise ValueError(
                f"Unknown extinction law: {ext!r}. "
                f"Expected one of {sorted(BUILTIN_LAWS)} or a callable."
            )
        law = BUILTIN_LAWS[ext]
        return lambda lam_cm, _params=None: law(lam_cm)

    @staticmethod
    def _require_positive_flux(flux, err) -> None:
        """Guard the log-flux likelihood: ``ln(F)`` and ``err/F`` need positive inputs."""
        if np.any(flux <= 0) or np.any(err <= 0):
            raise ValueError(
                "the log-flux likelihood requires strictly positive fluxes and errors"
            )

    def _consolidate_data(self):
        """Consolidate point observation data into single arrays."""
        if self._all_t is not None:
            return

        for bd in self._band_obs:
            self._require_positive_flux(bd.flux, bd.err)

        if self._point_t:
            self._all_t = np.concatenate(self._point_t)
            self._all_nu = np.concatenate(self._point_nu)
            self._all_flux = np.concatenate(self._point_flux)
            self._all_err = np.concatenate(self._point_err)
            self._all_weights = np.concatenate(self._point_weights)
            self._all_dataset = np.concatenate(
                [np.full(len(t), i, dtype=np.intp) for i, t in enumerate(self._point_t)]
            )

            # Sort by time (flux_density requires ascending time order)
            order = np.argsort(self._all_t)
            self._all_t = self._all_t[order]
            self._all_nu = self._all_nu[order]
            self._all_flux = self._all_flux[order]
            self._all_err = self._all_err[order]
            self._all_weights = self._all_weights[order]
            self._all_dataset = self._all_dataset[order]

            # Normalize weights so they sum to N
            w_sum = self._all_weights.sum()
            if w_sum > 0:
                self._all_weights *= len(self._all_weights) / w_sum

            self._require_positive_flux(self._all_flux, self._all_err)
            self._all_log_flux = np.log(self._all_flux)
            self._all_log_err = self._all_err / self._all_flux
            self._all_prec = self._all_weights / self._all_log_err**2
            self._dataset_inv_sys2 = self._inv_sys2(self._point_sys)
            self._dataset_prec_sum = np.bincount(
                self._all_dataset, weights=self._all_prec, minlength=len(self._point_t)
            )
        else:
            self._all_t = np.array([])
            self._all_nu = np.array([])

        # Upper limits, sorted by time.
        if self._limit_t:
            self._lim_t = np.concatenate(self._limit_t)
            lim_nu = np.concatenate(self._limit_nu)
            self._lim_value = np.concatenate(self._limit_value)
            self._lim_sigma = np.concatenate(self._limit_sigma)
            order = np.argsort(self._lim_t)
            self._lim_t = self._lim_t[order]
            lim_nu = lim_nu[order]
            self._lim_value = self._lim_value[order]
            self._lim_sigma = self._lim_sigma[order]
        else:
            self._lim_t = np.array([])
            lim_nu = np.array([])

        # One model evaluation serves detections and limits: merge, sort by
        # time, and remember where each group landed.
        eval_t = np.concatenate([self._all_t, self._lim_t])
        eval_nu = np.concatenate([self._all_nu, lim_nu])
        order = np.argsort(eval_t, kind="stable")
        self._eval_t = eval_t[order]
        self._eval_nu = eval_nu[order]
        position = np.empty_like(order)
        position[order] = np.arange(len(order))
        n_det = len(self._all_t)
        self._det_pos = position[:n_det]
        self._lim_pos = position[n_det:]

        # Precompute extinction kernel from rest-frame wavelengths.
        # Built-in laws collapse k(lambda) and the 0.4*ln(10) factor into
        # one cached array (one np.exp + one in-place multiply per step).
        # Custom laws may depend on params, so we only cache the wavelengths.
        if self._ext_law is not None and len(self._eval_t) > 0:
            lam_rest_cm = (_C_CGS / self._eval_nu) / (1.0 + self.z)
            if self._custom_extinction:
                self._lam_rest_cm = lam_rest_cm
            else:
                k = self._ext_law(lam_rest_cm, None)
                self._ext_kernel = _LN10_OVER_2P5 * k

    # ── Model evaluation (hot path) ─────────────────────────────────────────────────────────

    def _build_model(self, params: ModelParams) -> Model:
        """Create a Model instance from MCMC parameters."""
        jet = self.jet_factory(params)
        medium = self.medium_factory(params)

        observer = Observer(
            lumi_dist=self.lumi_dist,
            z=self.z,
            theta_obs=params.theta_v,
        )

        fwd_rad = Radiation(
            eps_e=params.eps_e,
            eps_B=params.eps_B,
            p=params.p,
            xi_e=params.xi_e,
            ssc=self.fwd_ssc,
            kn=self.kn,
        )

        rvs_rad = None
        if self.rvs_shock:
            rvs_rad = Radiation(
                eps_e=params.eps_e_r,
                eps_B=params.eps_B_r,
                p=params.p_r,
                xi_e=params.xi_e_r,
                ssc=self.rvs_ssc,
                kn=self.kn,
            )

        return Model(
            jet=jet,
            medium=medium,
            observer=observer,
            fwd_rad=fwd_rad,
            rvs_rad=rvs_rad,
            resolutions=self.resolution,
            rtol=self.rtol,
            radiative_fireball=self.radiative_fireball,
        )

    @staticmethod
    def _inv_sys2(sys_errs) -> np.ndarray:
        """Inverse prior variance of each dataset's ln-calibration factor (inf for none)."""
        s = np.asarray(sys_errs, dtype=np.float64)
        with np.errstate(divide="ignore"):
            return np.where(s > 0, 1.0 / np.maximum(s, 1e-300) ** 2, np.inf)

    @staticmethod
    def _chi2_sum(log_obs, log_err, model_flux, weights) -> float:
        """Weighted chi-squared sum in ln-flux space."""
        diff = log_obs - np.log(np.maximum(model_flux, 1e-300))
        return float(np.sum(weights * (diff / log_err) ** 2))

    @staticmethod
    def _chi2_marginal(log_obs, prec, model_flux, dataset, inv_sys2, prec_sum) -> float:
        """ln-flux chi-squared with each dataset's calibration factor marginalised.

        A shared factor c per dataset with prior N(0, sys_err^2) in ln-flux
        integrates out exactly: chi2 = sum(p r^2) - S1^2 / (1/sys_err^2 + S0),
        with S0 = sum(p), S1 = sum(p r) over the dataset. sys_err = 0 gives
        1/sys_err^2 = inf and the plain sum.
        """
        r = log_obs - np.log(np.maximum(model_flux, 1e-300))
        chi2 = float(np.sum(prec * r * r))
        s1 = np.bincount(dataset, weights=prec * r, minlength=len(inv_sys2))
        return chi2 - float(np.sum(s1 * s1 / (inv_sys2 + prec_sum)))

    @staticmethod
    def _band_dataset_chi2(bd, model_flux) -> float:
        prec = bd.weights / (bd.err / bd.flux) ** 2
        r = np.log(bd.flux) - np.log(np.maximum(model_flux, 1e-300))
        chi2 = float(np.sum(prec * r * r))
        if bd.sys_err > 0:
            s1 = float(np.sum(prec * r))
            chi2 -= s1 * s1 / (1.0 / bd.sys_err**2 + float(np.sum(prec)))
        return chi2

    def _evaluate(self, params: ModelParams) -> float:
        """Compute chi-squared for given parameters using Model directly."""
        self._consolidate_data()

        model = self._build_model(params)

        chi2 = 0.0

        if len(self._eval_t) > 0:
            eval_flux = self._point_model_flux(model, params)
            if len(self._all_t) > 0:
                chi2 += self._chi2_marginal(
                    self._all_log_flux,
                    self._all_prec,
                    eval_flux[self._det_pos],
                    self._all_dataset,
                    self._dataset_inv_sys2,
                    self._dataset_prec_sum,
                )
            if len(self._lim_t) > 0:
                chi2 += self._limits_chi2(
                    self._lim_value, self._lim_sigma, eval_flux[self._lim_pos]
                )

        # Band-integrated flux
        for bd in self._band_obs:
            model_flux = np.asarray(
                model.flux(bd.t, bd.nu_min, bd.nu_max, bd.num_points).total
            )
            chi2 += self._band_dataset_chi2(bd, model_flux)
        for bl in self._band_limits:
            model_flux = np.asarray(
                model.flux(bl.t, bl.nu_min, bl.nu_max, bl.num_points).total
            )
            chi2 += self._limits_chi2(bl.limit, bl.sigma, model_flux)

        return chi2

    @staticmethod
    def _limits_chi2(limit, sigma, model_flux) -> float:
        """-2 ln Phi((U - F_model) / sigma) summed over limits (linear flux)."""
        return float(-2.0 * np.sum(log_ndtr((limit - model_flux) / sigma)))

    def _point_model_flux(self, model, params: ModelParams) -> np.ndarray:
        """Model flux density at every evaluation point (detections and limits), extinction applied."""
        model_flux = np.asarray(model.flux_density(self._eval_t, self._eval_nu).total)
        if self._ext_law is not None and params.A_V != 0.0:
            if self._ext_kernel is not None:
                model_flux = model_flux * np.exp(-params.A_V * self._ext_kernel)
            else:
                k = self._ext_law(self._lam_rest_cm, params)
                model_flux = model_flux * np.exp((-params.A_V * _LN10_OVER_2P5) * k)
        return model_flux

    def _dataset_label(self, i: int) -> str:
        lbl = self._point_labels[i]
        if lbl:
            return lbl
        nu = np.unique(self._point_nu[i])
        if nu.size == 1:
            return f"{float(nu[0]):.3g} Hz"
        return f"spectrum at {float(self._point_t[i][0]):.3g} s"

    def calibration_offsets(self, params: ModelParams) -> List[SimpleNamespace]:
        """Best-fit calibration shift of every dataset given with ``sys_err > 0``.

        Returns one entry per such dataset (in the order the data were
        added, point datasets first, then band datasets) with
        ``label``, ``offset`` = data level relative to the model minus one
        (``+0.03`` means the data sit 3% above the model, i.e. the fit
        prefers scaling that dataset down by 3%), and ``sys_err``.
        """
        self._consolidate_data()
        model = self._build_model(params)
        out: List[SimpleNamespace] = []
        if len(self._all_t) > 0 and np.any(np.isfinite(self._dataset_inv_sys2)):
            model_flux = self._point_model_flux(model, params)[self._det_pos]
            r = self._all_log_flux - np.log(np.maximum(model_flux, 1e-300))
            s1 = np.bincount(
                self._all_dataset,
                weights=self._all_prec * r,
                minlength=len(self._point_t),
            )
            c_hat = s1 / (self._dataset_inv_sys2 + self._dataset_prec_sum)
            for i, sys_err in enumerate(self._point_sys):
                if sys_err > 0:
                    out.append(
                        SimpleNamespace(
                            label=self._dataset_label(i),
                            offset=float(np.expm1(c_hat[i])),
                            sys_err=float(sys_err),
                        )
                    )
        for j, bd in enumerate(self._band_obs):
            if bd.sys_err > 0:
                model_flux = np.asarray(
                    model.flux(bd.t, bd.nu_min, bd.nu_max, bd.num_points).total
                )
                prec = bd.weights / (bd.err / bd.flux) ** 2
                r = np.log(bd.flux) - np.log(np.maximum(model_flux, 1e-300))
                c_hat = float(np.sum(prec * r)) / (
                    1.0 / bd.sys_err**2 + float(np.sum(prec))
                )
                out.append(
                    SimpleNamespace(
                        label=bd.name or f"band {bd.nu_min:.3g}-{bd.nu_max:.3g} Hz",
                        offset=float(np.expm1(c_hat)),
                        sys_err=float(bd.sys_err),
                    )
                )
        return out

    # ── Parameter validation & sampler setup ─────────────────────────────────────────────────────────

    def validate_parameters(self, param_defs: Sequence[ParamDef]) -> None:
        """Validate parameter definitions against the current configuration."""
        from .params import validate_parameters

        return validate_parameters(self, param_defs)

    # ── Sampler dispatch ─────────────────────────────────────────────────────────

    def fit(
        self,
        param_defs: Sequence[ParamDef],
        sampler: str = "emcee",
        resolution: Optional[Tuple[float, float, float]] = None,
        npool: Optional[int] = None,
        top_k: int = 10,
        outdir: str = "bilby_output",
        label: str = "afterglow",
        clean: bool = True,
        resume: bool = False,
        log_likelihood_fn: Optional[Callable] = None,
        priors: Optional[dict] = None,
        **sampler_kwargs,
    ) -> FitResult:
        """Run sampler for parameter estimation.

        Args:
            param_defs: Parameter definitions
            sampler: Sampler name ('emcee', 'dynesty', etc.)
            resolution: Optional (phi, theta, t) override; persists like the
                constructor value (see class docstring "Resolution lifecycle").
            npool: Number of threads for parallelism (default: auto-detect)
            top_k: Number of top samples to return
            outdir: Output directory for bilby
            label: Run label for bilby
            clean: Clean previous runs (bilby)
            resume: Resume from previous run (bilby)
            log_likelihood_fn: Custom likelihood transform.
                Signature: (chi2: float) -> float. Default: lambda chi2: -0.5 * chi2
            priors: Custom prior distributions dict. Keys are parameter labels
                (using ``log10_`` prefix for LOG-scale params), values are
                ``bilby.core.prior.Prior`` objects. Parameters not in this dict
                automatically get Uniform priors from ParamDef bounds.
                Works with all samplers (emcee, dynesty, etc.).
            **sampler_kwargs: Additional sampler arguments
        """
        from .params import build_prior_dict, build_sampler_params
        from .samplers import fit_bilby, fit_emcee

        self.validate_parameters(param_defs)
        defs = list(param_defs)
        self._param_defs = defs

        # Per-fit resolution override persists so post-fit predictions
        # (``flux_density_grid``, ``flux``) render at the SAME resolution used
        # to compute chi-squared.
        if resolution is not None:
            self.resolution = tuple(resolution)

        # Consolidate before dispatching the worker pool so the values worker
        # threads read are stable.
        self._consolidate_data()

        if len(self._eval_t) == 0 and not self._band_obs and not self._band_limits:
            raise ValueError(
                "No observation data. Use add_flux_density(), add_spectrum(), or add_flux()."
            )

        if log_likelihood_fn is None:

            def log_likelihood_fn(chi2):
                return -0.5 * chi2

        labels, pl, pu, ndim = build_sampler_params(self, defs)
        self._to_params = _build_transformer(defs)
        prior_dict = build_prior_dict(labels, pl, pu, defs, priors)

        if sampler.lower() == "emcee":
            result = fit_emcee(
                self,
                defs,
                labels,
                pl,
                pu,
                ndim,
                top_k,
                prior_dict,
                log_likelihood_fn,
                npool,
                **sampler_kwargs,
            )
        else:
            result = fit_bilby(
                self,
                defs,
                labels,
                ndim,
                sampler,
                npool,
                top_k,
                outdir,
                label,
                clean,
                resume,
                log_likelihood_fn,
                prior_dict,
                **sampler_kwargs,
            )
        # Populate fit-quality bookkeeping so FitResult.summary() can surface
        # reduced χ² / BIC / AIC without the user counting data points.
        result.n_data, result.n_limits = self._data_counts()
        result.n_free_params = sum(1 for d in defs if d.scale is not Scale.fixed)
        if any(s > 0 for s in self._point_sys) or any(
            b.sys_err > 0 for b in self._band_obs
        ):
            result.calibration_offsets = [
                vars(o)
                for o in self.calibration_offsets(
                    self._to_params(result.top_k_params[0])
                )
            ]

        # Cache for convenience: matches Fitter.load(...) post-condition.
        self.result = result
        return result

    # ── Lifecycle helpers ─────────────────────────────────────────────────────────

    @contextmanager
    def _override_resolution(self, resolution):
        """Temporarily override resolution, restoring on exit."""
        saved = self.resolution
        if resolution is not None:
            self.resolution = tuple(resolution)
        try:
            yield
        finally:
            self.resolution = saved

    def _require_fitted(self):
        if self._to_params is None:
            raise RuntimeError(
                "Call .fit(...) first, or load a saved fit with Fitter.load(path)."
            )

    # ── Persistence (save / load) ─────────────────────────────────────────────────────────

    def save(self, path) -> None:
        """Persist this fit to a bilby-native HDF5 / JSON file.

        Stores the posterior, top-K best-fit parameters, the original
        ``ParamDef`` list, and a snapshot of the fitter configuration
        (constructor args + added observation data), so that
        :py:meth:`Fitter.load` can rebuild everything in one step. The
        on-disk format remains interoperable with bilby tooling -- the file
        can still be opened directly via ``bilby.read_in_result(path)`` for
        inspection-only access.

        Parameters
        ----------
        path : str | os.PathLike
            Output filename. Extension determines the format: ``.h5`` /
            ``.hdf5`` (recommended; smaller, faster) or ``.json``
            (human-readable). If no extension is provided, HDF5 is used.

        Raises
        ------
        RuntimeError
            If ``.fit(...)`` has not been called yet (nothing to save).
        """
        from .io import save_fitter

        return save_fitter(self, path)

    @classmethod
    def load(cls, path, *, jet=None, medium=None, extinction=None) -> "Fitter":
        """Reload a saved fit as a fully-configured ``Fitter``.

        Symmetric counterpart to :py:meth:`Fitter.save`: reads the file and
        reconstructs the original ``Fitter`` (constructor args, observation
        data, parameter transformer) along with the ``FitResult``. The returned
        fitter is immediately ready for prediction calls -- no manual
        reconfiguration, no MCMC re-run.

        Typical reload workflow::

            from VegasAfterglow import Fitter
            fitter = Fitter.load("vegas_mcmc_fit.h5")     # one step
            result = fitter.result
            lc = fitter.flux_density_grid(result.top_k_params[0], t, nu).total

        Custom (callable) ``jet`` / ``medium`` / ``extinction`` can't round-trip
        through serialisation; pass the same callable via
        ``Fitter.load(path, jet=...)`` in that case.

        Parameters
        ----------
        path : str | os.PathLike
            Path to a previously-saved HDF5 / JSON file.
        jet, medium, extinction : optional
            Overrides; required if the original used a custom callable for
            these, ignored otherwise.

        Returns
        -------
        Fitter
            Fully-configured fitter, with ``.result`` set to the loaded
            ``FitResult``.

        Raises
        ------
        ValueError
            If the file was written by an older release that doesn't include
            the fitter snapshot, or if the original used a custom callable
            that wasn't provided as an override.
        """
        from .io import load_fitter

        return load_fitter(path, jet=jet, medium=medium, extinction=extinction)

    # ── Post-fit prediction ─────────────────────────────────────────────────────────

    def _model_from_fit(self, best_params):
        """Build a Model from sampler-space parameters (requires prior fit)."""
        self._require_fitted()
        return self._build_model(self._to_params(best_params))

    def _apply_extinction_to_grid(self, result, nu, params):
        """Multiply each component of ``result`` by ``exp(-A_V*0.4*ln10*k(lam_rest))``
        along the nu axis; pass-through when extinction is disabled or ``A_V == 0``."""
        if self._ext_law is None or params.A_V == 0.0:
            return result

        lam_rest = (_C_CGS / np.asarray(nu)) / (1.0 + self.z)
        k = self._ext_law(lam_rest, params)
        attn = np.exp(-params.A_V * _LN10_OVER_2P5 * k)  # shape (n_nu,)

        def _attn(arr):
            a = np.asarray(arr)
            if a.size == 0:
                return a
            # PyFlux components are shape (n_nu, n_t); broadcast attn along axis 0
            return a * attn.reshape((-1,) + (1,) * (a.ndim - 1))

        return SimpleNamespace(
            total=_attn(result.total),
            fwd=SimpleNamespace(sync=_attn(result.fwd.sync), ssc=_attn(result.fwd.ssc)),
            rvs=SimpleNamespace(sync=_attn(result.rvs.sync), ssc=_attn(result.rvs.ssc)),
        )

    def flux_density_grid(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu: np.ndarray,
        resolution: Optional[Tuple[float, float, float]] = None,
    ):
        """Compute flux density grid at best-fit parameters.

        Host-galaxy dust extinction (when ``Fitter(extinction=...)`` was
        configured and the best-fit ``A_V != 0``) is applied per-nu, matching
        the fitter's per-frequency chi-squared path. The return type then
        becomes a ``SimpleNamespace`` with the same ``.total / .fwd / .rvs``
        layout as ``PyFlux``; in the no-extinction case the underlying
        ``PyFlux`` is returned directly.

        Args:
            best_params: Best-fit parameter array (in sampler space)
            t: Time array [seconds]
            nu: Frequency array [Hz]
            resolution: Optional resolution override (phi, theta, t)
        """
        with self._override_resolution(resolution):
            params = self._to_params(best_params)
            result = self._build_model(params).flux_density_grid(t, nu)
            return self._apply_extinction_to_grid(result, nu, params)

    def flux(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        band,
        num_points: int = 5,
        resolution: Optional[Tuple[float, float, float]] = None,
    ):
        """Compute integrated flux at best-fit parameters.

        Extinction is not applied; use ``flux_density_grid`` for extincted
        predictions (matches the band-integrated chi-squared path).

        Args:
            best_params: Best-fit parameter array (in sampler space)
            t: Time array [seconds]
            band: Frequency band as (nu_min, nu_max) tuple in Hz.
                Use ``units.band("XRT")`` for named instrument bands.
            num_points: Frequency sampling points (Boole's rule; best with
                4k+1 — 5 covers optical/X-ray bands, 9-17 for radio or
                wider than two decades)
            resolution: Optional resolution override (phi, theta, t)
        """
        nu_min, nu_max = band
        with self._override_resolution(resolution):
            return self._model_from_fit(best_params).flux(t, nu_min, nu_max, num_points)

    def _has_posterior_samples(self) -> bool:
        """True when ``self.result`` carries usable posterior samples (set by
        ``.fit()`` or loaded via ``Fitter.load``)."""
        return (
            self.result is not None
            and self.result.samples is not None
            and self.result.samples.size > 0
        )

    def _draw_posterior(
        self, n_samples: int, rng: Optional[np.random.Generator]
    ) -> np.ndarray:
        """Sample ``n_samples`` rows from the flat posterior. Sampling is with
        replacement only when ``n_samples`` exceeds the posterior size."""
        if not self._has_posterior_samples():
            raise RuntimeError(
                "credible bands require posterior samples; call .fit() first"
            )
        if n_samples < 2:
            raise ValueError(f"n_samples must be >= 2, got {n_samples}")
        rng = np.random.default_rng() if rng is None else rng
        flat = self.result.flat_samples
        idx = rng.choice(
            flat.shape[0], size=n_samples, replace=(n_samples > flat.shape[0])
        )
        return flat[idx]

    def _draw_chi2_region(
        self, cl: float, n_samples: int, rng: Optional[np.random.Generator]
    ) -> np.ndarray:
        """Rows of the flat posterior inside the joint chi-squared confidence
        region ``chi2 <= chi2_min + delta_chi2(cl, n_free)``. The offset form
        makes the selection independent of any constant in the stored
        log-probabilities. The minimum-chi2 sample is always included so the
        best-fit curve lies inside the resulting envelope."""
        if not self._has_posterior_samples():
            raise RuntimeError(
                "confidence bands require posterior samples; call .fit() first"
            )
        if not (0.0 < cl < 1.0):
            raise ValueError(f"cl must be in (0, 1), got {cl}")
        from scipy.stats import chi2 as _chi2_dist

        flat = self.result.flat_samples
        chi2 = -2.0 * self.result.log_probs.ravel()
        n_free = self.result.n_free_params or flat.shape[1]
        thresh = float(chi2.min()) + float(_chi2_dist.ppf(cl, df=n_free))
        region = flat[chi2 <= thresh]
        best = flat[int(np.argmin(chi2))]
        if len(region) > n_samples:
            rng = np.random.default_rng() if rng is None else rng
            idx = rng.choice(region.shape[0], size=n_samples, replace=False)
            region = region[idx]
        return np.vstack([best[None, :], region])

    def _credible_band(
        self,
        draws: np.ndarray,
        ci: float,
        evaluate: Callable[[np.ndarray], np.ndarray],
        n_workers: Optional[int],
        resolution: Optional[Tuple[float, float, float]],
    ) -> SimpleNamespace:
        """Apply ``evaluate(params)`` to each posterior draw, return the
        central-``ci`` percentile envelope from those draws plus the
        pointwise median curve of the same draws. Threads share a single
        ``_override_resolution`` context so the per-call resolution
        mutation doesn't race between draws."""
        if not (0.0 < ci <= 1.0):
            raise ValueError(f"ci must be in (0, 1], got {ci}")
        workers = n_workers if n_workers is not None else max(1, os.cpu_count() or 1)
        with self._override_resolution(resolution):
            if workers > 1:
                with ThreadPoolExecutor(max_workers=workers) as pool:
                    results = list(pool.map(evaluate, draws))
            else:
                results = [evaluate(p) for p in draws]

        stack = np.stack([np.asarray(r) for r in results], axis=0)
        lower, median, upper = np.percentile(
            stack, [50.0 * (1.0 - ci), 50.0, 50.0 * (1.0 + ci)], axis=0
        )
        return SimpleNamespace(lower=lower, median=median, upper=upper)

    def flux_density_credible(
        self,
        t: np.ndarray,
        nu: np.ndarray,
        ci: float = 0.68,
        n_samples: int = 200,
        n_workers: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> SimpleNamespace:
        """Posterior credible band on the total flux density.

        Draws ``n_samples`` random posterior samples, evaluates
        ``flux_density_grid`` at each, and returns the central ``ci``
        credible interval per ``(nu, t)`` cell. Extinction is applied
        per draw (matching ``flux_density_grid``) so the band reflects
        ``A_V`` uncertainty when ``A_V`` is sampled.

        Args:
            t: Time array [seconds]
            nu: Frequency array [Hz]
            ci: Credible interval (e.g. 0.68 for ~1σ, 0.95 for ~2σ)
            n_samples: Posterior draws to evaluate (larger = smoother)
            n_workers: Thread count for parallel evaluation. ``None`` -> ``os.cpu_count()``.
            rng: ``numpy.random.Generator``; defaults to ``np.random.default_rng()``
            resolution: Optional ``(phi, theta, t)`` override

        Returns:
            ``SimpleNamespace(lower, median, upper)`` with arrays of shape
            ``(len(nu), len(t))``. All three are per-cell percentiles of the
            ``n_samples`` posterior draws (``median`` is the 50th), so the
            median curve always lies inside the band. It represents the
            typical posterior prediction, not any single parameter set; for
            a best-fit overlay evaluate ``flux_density_grid`` at
            ``result.top_k_params[0]``.
            Usage: ``ax.fill_between(t, out.lower[i_nu], out.upper[i_nu])``
            followed by ``ax.plot(t, out.median[i_nu])``.
        """
        draws = self._draw_posterior(n_samples, rng)
        t_arr = np.atleast_1d(np.asarray(t, dtype=np.float64))
        nu_arr = np.atleast_1d(np.asarray(nu, dtype=np.float64))
        return self._credible_band(
            draws,
            ci,
            lambda p: np.asarray(self.flux_density_grid(p, t_arr, nu_arr).total),
            n_workers,
            resolution,
        )

    def flux_credible(
        self,
        t: np.ndarray,
        band,
        ci: float = 0.68,
        n_samples: int = 200,
        num_points: int = 5,
        n_workers: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> SimpleNamespace:
        """Posterior credible band on band-integrated flux.

        Same semantics as :meth:`flux_density_credible` but mirrors
        :meth:`flux` (band-integrated, no extinction applied).

        Returns:
            ``SimpleNamespace(lower, median, upper)`` with arrays of shape ``(len(t),)``.
        """
        try:
            _, _ = band
        except (TypeError, ValueError):
            raise ValueError(
                f"flux_credible: band must be a (nu_min, nu_max) tuple in Hz, got {band!r}"
            ) from None
        draws = self._draw_posterior(n_samples, rng)
        t_arr = np.atleast_1d(np.asarray(t, dtype=np.float64))
        return self._credible_band(
            draws,
            ci,
            lambda p: np.asarray(self.flux(p, t_arr, band, num_points).total),
            n_workers,
            resolution,
        )

    def flux_density_confidence(
        self,
        t: np.ndarray,
        nu: np.ndarray,
        cl: float = 0.68,
        n_samples: int = 200,
        n_workers: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> SimpleNamespace:
        """Joint chi-squared confidence envelope on the total flux density.

        The likelihood-consistent companion of the best-fit curve: selects
        every posterior sample with ``chi2 <= chi2_min + delta_chi2(cl,
        n_free)`` (the joint ``cl`` confidence region), evaluates the model
        at up to ``n_samples`` of them, and returns the pointwise min/max
        envelope of those curves. The minimum-chi2 sample is always
        included, so the best-fit curve lies inside the envelope by
        construction. The selection depends only on chi-squared
        *differences*, so it is insensitive to any constant offset in the
        stored log-probabilities.

        Compared to :meth:`flux_density_credible` (pointwise posterior
        percentiles), this band is a *joint* region — for many free
        parameters it is wider, and its edges are set by the extreme
        curves inside the region, so larger ``n_samples`` gives a more
        complete envelope.

        Args:
            t: Time array [seconds]
            nu: Frequency array [Hz]
            cl: Joint confidence level (e.g. 0.68 for ~1σ)
            n_samples: Maximum region samples to evaluate
            n_workers: Thread count for parallel evaluation. ``None`` -> ``os.cpu_count()``.
            rng: ``numpy.random.Generator``; defaults to ``np.random.default_rng()``
            resolution: Optional ``(phi, theta, t)`` override

        Returns:
            ``SimpleNamespace(lower, median, upper)`` with arrays of shape
            ``(len(nu), len(t))``. ``lower`` / ``upper`` are the envelope of
            region curves; ``median`` is their pointwise median.
        """
        draws = self._draw_chi2_region(cl, n_samples, rng)
        t_arr = np.atleast_1d(np.asarray(t, dtype=np.float64))
        nu_arr = np.atleast_1d(np.asarray(nu, dtype=np.float64))
        return self._credible_band(
            draws,
            1.0,
            lambda p: np.asarray(self.flux_density_grid(p, t_arr, nu_arr).total),
            n_workers,
            resolution,
        )

    def flux_confidence(
        self,
        t: np.ndarray,
        band,
        cl: float = 0.68,
        n_samples: int = 200,
        num_points: int = 5,
        n_workers: Optional[int] = None,
        rng: Optional[np.random.Generator] = None,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> SimpleNamespace:
        """Joint chi-squared confidence envelope on band-integrated flux.

        Same semantics as :meth:`flux_density_confidence` but mirrors
        :meth:`flux` (band-integrated, no extinction applied).

        Returns:
            ``SimpleNamespace(lower, median, upper)`` with arrays of shape ``(len(t),)``.
        """
        try:
            _, _ = band
        except (TypeError, ValueError):
            raise ValueError(
                f"flux_confidence: band must be a (nu_min, nu_max) tuple in Hz, got {band!r}"
            ) from None
        draws = self._draw_chi2_region(cl, n_samples, rng)
        t_arr = np.atleast_1d(np.asarray(t, dtype=np.float64))
        return self._credible_band(
            draws,
            1.0,
            lambda p: np.asarray(self.flux(p, t_arr, band, num_points).total),
            n_workers,
            resolution,
        )

    def model(
        self,
        best_params: np.ndarray,
    ) -> Model:
        """Return the raw underlying ``Model`` at best-fit parameters.

        No extinction is applied -- this is the escape hatch for users who
        want to drive the C++ flux interfaces directly. To get extincted
        predictions, use ``flux_density_grid`` instead.
        """
        return self._model_from_fit(best_params)

    def draw_fit(
        self,
        best_params: Optional[np.ndarray] = None,
        *,
        ci: float = 0.68,
        n_samples: int = 100,
        obs_noise: Literal["frac", "abs", "none"] = "none",
        t_range: Optional[Tuple[float, float]] = None,
        n_t: int = 500,
        shifts: Optional[dict] = None,
        auto_shift_gap: float = 1.0,
        show_nu_panel: bool = True,
        resolution: Optional[Tuple[float, float, float]] = None,
        fig=None,
        axes=None,
    ):
        """Diagnostic plot of observation data overlaid with the fitted model.

        Builds a two-panel figure:

        * **Top panel** -- data + model curves. Each band's central line is
          the **best-fit model** (``best_params``, defaulting to
          ``result.top_k_params[0]``); when posterior samples are available
          (post-``fit()``) it is surrounded by the joint chi-squared
          confidence envelope — curves from samples with
          ``chi2 <= chi2_min + delta_chi2(ci, n_free)`` — which is
          likelihood-consistent with the best-fit line and contains it by
          construction (see :meth:`flux_density_confidence`).
          By default (``obs_noise='none'``) the band is the pure
          model-curve envelope (pinches at the data, no observation
          noise added). Switch to ``obs_noise='frac'`` to broaden into a
          *posterior predictive* envelope where σ scales with the model
          flux as ``central * median(err/flux)`` — represents *"where
          would a new observation fall?"* and stays visually well-behaved
          on log-y plots that span many decades. Use ``obs_noise='abs'``
          for a constant absolute σ instead (suitable when noise is
          detector-limited).
          Single-frequency light curves (added via ``add_flux_density``) and
          band-integrated fluxes (added via ``add_flux``) get a dual y-axis
          if both kinds are present: left is ``F_nu`` in erg/cm^2/s/Hz,
          right is ``F`` in erg/cm^2/s. Multiple bands are shifted vertically
          (in log) so they don't overlap; the legend reports the shift factors.

        * **Bottom panel** (optional, ``show_nu_panel=True``) -- evolution of
          the characteristic synchrotron break frequencies
          (``nu_a``, ``nu_m``, ``nu_c``) in the **observer frame**, with the
          observed frequencies overlaid as horizontal lines and bands as shaded
          regions for visual cross-reference. Frequencies are evaluated at the
          jet theta-column closest to ``theta_v`` so the Doppler boost reflects
          what the observer actually sees off-axis (recovers the on-axis cell
          when ``theta_v == 0``).

        Parameters
        ----------
        best_params
            Sampler-space parameter array. Defaults to ``self.result.top_k_params[0]``.
            Sets the central model trajectory in the top panel and the
            break-frequency evolution in the bottom panel.
        ci
            Credible interval for the shaded band (e.g. 0.68 for ~1σ,
            0.95 for ~2σ). Ignored when ``n_samples == 0``.
        n_samples
            Posterior draws for the shaded band. ``0`` skips the band and
            plots the MAP trajectory only. Default ``100``.
        obs_noise
            Noise model used to broaden the confidence envelope into a
            posterior predictive band:

            * ``'none'`` (default) -- skip observation noise; the band is
              the pure model-curve envelope and pinches at the data.
            * ``'frac'`` -- σ scales with the model flux as
              ``central * median(err/flux)``; visually well-behaved on
              log-y plots over many decades. Use this for a goodness-of-fit
              display where ~``ci`` of the data points should fall inside
              the band.
            * ``'abs'`` -- constant absolute σ equal to the per-band median
              ``err``. Appropriate when the noise really is a constant floor
              (detector-noise-limited measurements).
        t_range
            (t_min, t_max) for the model curves in seconds. Defaults to one
            decade below ``tmin`` and two decades above ``tmax`` across all
            observation data.
        n_t
            Number of points in the model time grid.
        shifts
            Optional dict overriding auto-computed shifts. Keys are
            ``("lc", nu_hz)`` for light-curve bands or
            ``("band", (nu_min, nu_max))`` for band-integrated entries; values
            are ``log10`` multipliers.
        auto_shift_gap
            Decades added between consecutive bands (sorted by frequency).
            Each band is shifted by ``(rank - (n-1)/2) * auto_shift_gap`` so
            the legend shows small, symmetric ``× 10^{±k}`` factors. Bands
            with wide internal flux ranges can still overlap visually --
            raise this value to push them apart.
        show_nu_panel
            Whether to include the bottom break-frequency panel.
        resolution
            Optional ``(phi, theta, t)`` resolution override.
        fig, axes
            Optional existing figure / axes to draw into. If ``axes`` is
            provided it must be a tuple ``(ax_top, ax_bot)`` (``ax_bot`` may be
            ``None`` when ``show_nu_panel=False``).

        Returns
        -------
        fig, (ax_top, ax_bot)
            The matplotlib figure and a 2-tuple of axes. ``ax_bot`` is ``None``
            when ``show_nu_panel=False``.

        Raises
        ------
        ValueError
            If the fitter has no observation data, or if no fit has been run
            and ``best_params`` is None.

        Notes
        -----
        For ``add_spectrum`` data (broadband spectra at fixed times), v1 logs a
        warning and skips them -- the natural model curve geometry differs from
        light curves. Reverse-shock break frequencies are not drawn in v1
        (forward-shock only).
        """
        from ..plotting.diagnostic import draw_fit as _impl

        return _impl(
            self,
            best_params=best_params,
            ci=ci,
            n_samples=n_samples,
            obs_noise=obs_noise,
            t_range=t_range,
            n_t=n_t,
            shifts=shifts,
            auto_shift_gap=auto_shift_gap,
            show_nu_panel=show_nu_panel,
            resolution=resolution,
            fig=fig,
            axes=axes,
        )
