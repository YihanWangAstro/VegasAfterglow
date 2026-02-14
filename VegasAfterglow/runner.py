"""Afterglow model fitting with custom priors, likelihoods, and jet/medium profiles."""

import logging
import math
import os
import types
from concurrent.futures import ThreadPoolExecutor
from contextlib import contextmanager
from dataclasses import dataclass
from typing import Callable, List, Optional, Sequence, Tuple

import bilby
import emcee
import numpy as np
from bilby.core.sampler.emcee import Emcee as _BilbyEmcee

from .types import FitResult, ModelParams, ParamDef, Scale
from .VegasAfterglowC import (
    ISM,
    GaussianJet,
    Magnetar,
    Model,
    Observer,
    PowerLawJet,
    PowerLawWing,
    Radiation,
    StepPowerLawJet,
    TophatJet,
    TwoComponentJet,
    Wind,
)

# Patch bilby to accept 'moves' kwarg for emcee sampler
_BilbyEmcee.default_kwargs["moves"] = None


class ThreadPoolWithClose(ThreadPoolExecutor):
    def close(self):
        self.shutdown(wait=True)

    def join(self):
        pass


logger = logging.getLogger(__name__)

LATEX_LABELS = {
    "E_iso": r"$E_{\rm iso}$",
    "Gamma0": r"$\Gamma_0$",
    "theta_c": r"$\theta_c$",
    "theta_v": r"$\theta_v$",
    "theta_w": r"$\theta_w$",
    "k_e": r"$k_E$",
    "k_g": r"$k_\Gamma$",
    "E_iso_w": r"$E_{\rm iso,w}$",
    "Gamma0_w": r"$\Gamma_{0,w}$",
    "n_ism": r"$n_{\rm ISM}$",
    "A_star": r"$A_*$",
    "n0": r"$n_0$",
    "k_m": r"$k_m$",
    "p": r"$p$",
    "eps_e": r"$\epsilon_e$",
    "eps_B": r"$\epsilon_B$",
    "xi_e": r"$\xi_e$",
    "p_r": r"$p_r$",
    "eps_e_r": r"$\epsilon_{e,r}$",
    "eps_B_r": r"$\epsilon_{B,r}$",
    "xi_e_r": r"$\xi_{e,r}$",
    "tau": r"$\tau$",
    "L0": r"$L_0$",
    "t0": r"$t_0$",
    "q": r"$q$",
}

MEDIUM_RULES = {
    "ism": ({"n_ism"}, {"A_star", "n0", "k_m"}),
    "wind": ({"A_star"}, set()),
}

JET_RULES = {
    "tophat": ({"theta_c", "E_iso", "Gamma0"}, {"k_e", "k_g", "E_iso_w", "Gamma0_w"}),
    "gaussian": ({"theta_c", "E_iso", "Gamma0"}, {"k_e", "k_g", "E_iso_w", "Gamma0_w"}),
    "powerlaw": ({"theta_c", "E_iso", "Gamma0", "k_e", "k_g"}, {"E_iso_w", "Gamma0_w"}),
    "two_component": (
        {"theta_c", "E_iso", "Gamma0", "theta_w", "E_iso_w", "Gamma0_w"},
        {"k_e", "k_g"},
    ),
    "step_powerlaw": (
        {"theta_c", "E_iso", "Gamma0", "E_iso_w", "Gamma0_w", "k_e", "k_g"},
        set(),
    ),
    "powerlaw_wing": (
        {"theta_c", "E_iso_w", "Gamma0_w", "k_e", "k_g"},
        {"E_iso", "Gamma0"},
    ),
}

TOGGLE_RULES = {
    "forward_shock": ({"eps_e", "eps_B", "p"}, set()),
    "rvs_shock": (
        {"p_r", "eps_e_r", "eps_B_r", "tau"},
        {"p_r", "eps_e_r", "eps_B_r", "xi_e_r"},
    ),
    "magnetar": ({"L0", "t0", "q"}, {"L0", "t0", "q"}),
}

SAMPLER_DEFAULTS = {
    "dynesty": {
        "nlive": 500,
        "dlogz": 0.1,
        "sample": "rslice",
        "maxmcmc": 5000,
    },
    "emcee": {
        "nsteps": 5000,
        "nburn": 1000,
        "thin": 1,
        "moves": [
            (emcee.moves.DEMove(), 0.7),
            (emcee.moves.DESnookerMove(), 0.3),
        ],
    },
}


# --- Default Jet/Medium Factories ---


def _default_jet_factory(fitter: "Fitter"):
    """Build a jet factory from fitter's jet type string."""

    def factory(params: ModelParams):
        jet_type = fitter.jet
        spreading = False
        duration = params.tau if hasattr(params, "tau") else 1.0

        if fitter.magnetar:
            magnetar = Magnetar(L0=params.L0, t0=params.t0, q=params.q)
        else:
            magnetar = None

        if jet_type == "tophat":
            return TophatJet(
                theta_c=params.theta_c,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        elif jet_type == "gaussian":
            return GaussianJet(
                theta_c=params.theta_c,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        elif jet_type == "powerlaw":
            return PowerLawJet(
                theta_c=params.theta_c,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                k_e=params.k_e,
                k_g=params.k_g,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        elif jet_type == "two_component":
            return TwoComponentJet(
                theta_c=params.theta_c,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                theta_w=params.theta_w,
                E_iso_w=params.E_iso_w,
                Gamma0_w=params.Gamma0_w,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        elif jet_type == "step_powerlaw":
            return StepPowerLawJet(
                theta_c=params.theta_c,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                E_iso_w=params.E_iso_w,
                Gamma0_w=params.Gamma0_w,
                k_e=params.k_e,
                k_g=params.k_g,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        elif jet_type == "powerlaw_wing":
            return PowerLawWing(
                theta_c=params.theta_c,
                E_iso_w=params.E_iso_w,
                Gamma0_w=params.Gamma0_w,
                k_e=params.k_e,
                k_g=params.k_g,
                spreading=spreading,
                duration=duration,
            )
        elif jet_type == "uniform":
            return TophatJet(
                theta_c=np.pi / 2,
                E_iso=params.E_iso,
                Gamma0=params.Gamma0,
                spreading=spreading,
                duration=duration,
                magnetar=magnetar,
            )
        else:
            raise ValueError(f"Unknown jet type: {jet_type}")

    return factory


def _default_medium_factory(fitter: "Fitter"):
    """Build a medium factory from fitter's medium type string."""

    def factory(params: ModelParams):
        if fitter.medium == "ism":
            return ISM(n_ism=params.n_ism)
        elif fitter.medium == "wind":
            return Wind(
                A_star=params.A_star, n_ism=params.n_ism, n0=params.n0, k=params.k_m
            )
        else:
            raise ValueError(f"Unknown medium type: {fitter.medium}")

    return factory


@dataclass
class _BandObs:
    """Band-integrated flux observation group."""

    nu_min: float
    nu_max: float
    num_points: int
    t: np.ndarray
    flux: np.ndarray
    err: np.ndarray
    weights: np.ndarray


# --- Helper Functions ---


def get_optimal_nwalkers(ndim: int, ncpu: Optional[int] = None) -> int:
    """Compute optimal nwalkers for emcee.

    Ensures nwalkers >= 4 * ndim and is a multiple of 2 * ncpu
    (emcee splits walkers in half, so each half should fill CPUs evenly).
    """
    if ncpu is None:
        ncpu = os.cpu_count() or 1

    math_floor = 4 * ndim
    align_unit = 2 * ncpu
    units_needed = math.ceil(math_floor / align_unit)
    units_needed = max(1, units_needed)

    return units_needed * align_unit


def get_optimal_queue_size(ncpu, nlive):
    """Compute optimal queue_size for dynesty based on hardware and sampling parameters."""
    target_size = ncpu * 2
    max_safe_size = int(nlive * 0.15)
    optimal_size = min(target_size, max_safe_size)
    aligned_size = (optimal_size // ncpu) * ncpu

    if aligned_size < ncpu:
        aligned_size = ncpu
    return aligned_size


def _get_latex_label(param_def: ParamDef) -> str:
    """Generate LaTeX label for parameter (with log10 wrapper if needed)."""
    base_latex = LATEX_LABELS.get(param_def.name, param_def.name)
    if param_def.scale is Scale.LOG:
        return rf"$\log_{{10}}({base_latex.strip('$')})$"
    return base_latex


def _get_model_params_defaults() -> dict:
    """Get default values for all ModelParams fields."""
    mp = ModelParams()
    defaults = {}
    for attr in dir(mp):
        if not attr.startswith("_"):
            try:
                defaults[attr] = getattr(mp, attr)
            except Exception:
                pass
    return defaults


_MODEL_PARAMS_DEFAULTS = None


def _build_transformer(param_defs: List[ParamDef]):
    """Build a parameter transformer from sampler array to parameter namespace.

    Standard ModelParams fields get their defaults; ParamDef values override them.
    Custom parameters (e.g., r_scale) are also supported.
    """
    global _MODEL_PARAMS_DEFAULTS
    if _MODEL_PARAMS_DEFAULTS is None:
        _MODEL_PARAMS_DEFAULTS = _get_model_params_defaults()

    free_mappings = []  # (name, is_log)
    fixed_values = []  # (name, value)
    for pd in param_defs:
        if pd.scale is Scale.FIXED:
            val = pd.initial if pd.initial is not None else pd.lower
            fixed_values.append((pd.name, val))
        else:
            free_mappings.append((pd.name, pd.scale is Scale.LOG))

    defaults = dict(_MODEL_PARAMS_DEFAULTS)

    def transformer(theta):
        params = types.SimpleNamespace(**defaults)
        for name, val in fixed_values:
            setattr(params, name, val)
        for i, (name, is_log) in enumerate(free_mappings):
            setattr(params, name, 10 ** theta[i] if is_log else theta[i])
        return params

    return transformer


# --- Bilby Likelihood ---


class AfterglowLikelihood(bilby.Likelihood):
    """Bilby-compatible likelihood using Model directly."""

    __slots__ = (
        "parameters",
        "param_keys",
        "_fitter",
        "_theta",
        "_transformer",
        "_log_likelihood_fn",
    )

    def __init__(
        self,
        fitter: "Fitter",
        param_defs: List[ParamDef],
        log_likelihood_fn,
        transformer,
    ):
        param_keys = tuple(
            f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name
            for pd in param_defs
            if pd.scale is not Scale.FIXED
        )
        super().__init__(parameters={key: None for key in param_keys})
        self.param_keys = param_keys
        self._fitter = fitter
        self._theta = np.empty(len(param_keys), dtype=np.float64)
        self._log_likelihood_fn = log_likelihood_fn
        self._transformer = transformer

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.__slots__}

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)

    def log_likelihood(self, parameters=None) -> float:
        if parameters is not None:
            self.parameters.update(parameters)
        for i, key in enumerate(self.param_keys):
            self._theta[i] = self.parameters[key]

        try:
            params = self._transformer(self._theta)
            chi2 = self._fitter._evaluate(params)
            if not np.isfinite(chi2):
                return -np.inf
            return self._log_likelihood_fn(chi2)
        except Exception:
            return -np.inf


# --- Fitter ---


class Fitter:
    """Model-based afterglow fitter with custom priors, likelihood, and jet/medium profiles.

    Uses Model directly for flux evaluation, with ThreadPoolExecutor for parallelism.
    Each thread creates a Model instance (brief GIL hold), then calls flux_density
    (GIL released during C++ computation).

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
        ssc_cooling: Enable IC cooling effects
        kn: Enable Klein-Nishina corrections
        magnetar: Enable magnetar energy injection
        rtol: Numerical tolerance
        resolution: Grid resolution tuple (phi, theta, t). Can be overridden per-fit.
    """

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
        ssc_cooling: bool = False,
        kn: bool = False,
        magnetar: bool = False,
        rtol: float = 1e-6,
        resolution: Tuple[float, float, float] = (0.15, 0.5, 10),
    ):
        self.z = z
        self.lumi_dist = lumi_dist
        self.fwd_ssc = fwd_ssc
        self.rvs_ssc = rvs_ssc
        self.rvs_shock = rvs_shock
        self.ssc_cooling = ssc_cooling
        self.kn = kn
        self.magnetar = magnetar
        self.rtol = rtol

        # Resolution (can be overridden per-fit)
        self.phi_resol, self.theta_resol, self.t_resol = resolution

        # Jet: string (built-in type) or callable (custom factory)
        self.jet = jet
        self._custom_jet = callable(jet)
        if self._custom_jet:
            self.jet_factory = jet
        else:
            self.jet_factory = _default_jet_factory(self)

        # Medium: string (built-in type) or callable (custom factory)
        self.medium = medium
        self._custom_medium = callable(medium)
        if self._custom_medium:
            self.medium_factory = medium
        else:
            self.medium_factory = _default_medium_factory(self)

        self._to_params = None

        # Observation data in user units (seconds, Hz, erg/cm^2/s/Hz)
        self._point_t = []
        self._point_nu = []
        self._point_flux = []
        self._point_err = []
        self._point_weights = []
        self._band_obs: List[_BandObs] = []

        # Consolidated arrays (built lazily)
        self._all_t = None
        self._all_nu = None
        self._all_flux = None
        self._all_err = None
        self._all_weights = None

    def add_flux_density(
        self,
        nu: float,
        t: np.ndarray,
        f_nu: np.ndarray,
        err: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ):
        """Add light curve data at a single frequency.

        Args:
            nu: Observing frequency [Hz]
            t: Observation times [seconds]
            f_nu: Observed flux densities [erg/cm^2/s/Hz]
            err: Observational uncertainties [erg/cm^2/s/Hz]
            weights: Optional statistical weights
        """
        t = np.asarray(t, dtype=np.float64)
        f_nu = np.asarray(f_nu, dtype=np.float64)
        err = np.asarray(err, dtype=np.float64)
        w = (
            np.asarray(weights, dtype=np.float64)
            if weights is not None
            else np.ones_like(t)
        )

        self._point_t.append(t)
        self._point_nu.append(np.full_like(t, nu))
        self._point_flux.append(f_nu)
        self._point_err.append(err)
        self._point_weights.append(w)
        self._all_t = None  # Invalidate cache

    def add_spectrum(
        self,
        t: float,
        nu: np.ndarray,
        f_nu: np.ndarray,
        err: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ):
        """Add a broadband spectrum at a specific time.

        Args:
            t: Observation time [seconds]
            nu: Observing frequencies [Hz]
            f_nu: Observed flux densities [erg/cm^2/s/Hz]
            err: Observational uncertainties [erg/cm^2/s/Hz]
            weights: Optional statistical weights
        """
        nu = np.asarray(nu, dtype=np.float64)
        f_nu = np.asarray(f_nu, dtype=np.float64)
        err = np.asarray(err, dtype=np.float64)
        w = (
            np.asarray(weights, dtype=np.float64)
            if weights is not None
            else np.ones_like(nu)
        )

        self._point_t.append(np.full_like(nu, t))
        self._point_nu.append(nu)
        self._point_flux.append(f_nu)
        self._point_err.append(err)
        self._point_weights.append(w)
        self._all_t = None

    def add_flux(
        self,
        nu_min: float,
        nu_max: float,
        num_points: int,
        t: np.ndarray,
        flux: np.ndarray,
        err: np.ndarray,
        weights: Optional[np.ndarray] = None,
    ):
        """Add band-integrated flux measurements.

        Args:
            nu_min: Lower frequency bound [Hz]
            nu_max: Upper frequency bound [Hz]
            num_points: Number of frequency sampling points for integration
            t: Observation times [seconds]
            flux: Observed integrated fluxes [erg/cm^2/s]
            err: Observational uncertainties [erg/cm^2/s]
            weights: Optional statistical weights
        """
        t = np.asarray(t, dtype=np.float64)
        flux = np.asarray(flux, dtype=np.float64)
        err = np.asarray(err, dtype=np.float64)
        w = (
            np.asarray(weights, dtype=np.float64)
            if weights is not None
            else np.ones_like(t)
        )

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
            )
        )

    def _consolidate_data(self):
        """Consolidate point observation data into single arrays."""
        if self._all_t is not None:
            return

        if self._point_t:
            self._all_t = np.concatenate(self._point_t)
            self._all_nu = np.concatenate(self._point_nu)
            self._all_flux = np.concatenate(self._point_flux)
            self._all_err = np.concatenate(self._point_err)
            self._all_weights = np.concatenate(self._point_weights)

            # Sort by time (flux_density requires ascending time order)
            order = np.argsort(self._all_t)
            self._all_t = self._all_t[order]
            self._all_nu = self._all_nu[order]
            self._all_flux = self._all_flux[order]
            self._all_err = self._all_err[order]
            self._all_weights = self._all_weights[order]

            # Normalize weights so they sum to N
            w_sum = self._all_weights.sum()
            if w_sum > 0:
                self._all_weights *= len(self._all_weights) / w_sum
        else:
            self._all_t = np.array([])

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
            ssc_cooling=self.ssc_cooling,
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
                ssc_cooling=self.ssc_cooling,
                ssc=self.rvs_ssc,
                kn=self.kn,
            )

        return Model(
            jet=jet,
            medium=medium,
            observer=observer,
            fwd_rad=fwd_rad,
            rvs_rad=rvs_rad,
            resolutions=(self.phi_resol, self.theta_resol, self.t_resol),
            rtol=self.rtol,
        )

    def _evaluate(self, params: ModelParams) -> float:
        """Compute chi-squared for given parameters using Model directly."""
        self._consolidate_data()

        model = self._build_model(params)

        chi2 = 0.0

        if len(self._all_t) > 0:
            flux_result = model.flux_density(self._all_t, self._all_nu)
            model_flux = np.asarray(flux_result.total)
            diff = self._all_flux - model_flux
            chi2 += float(np.sum(self._all_weights * (diff / self._all_err) ** 2))

        # Band-integrated flux
        for bd in self._band_obs:
            model_flux = np.asarray(
                model.flux(bd.t, bd.nu_min, bd.nu_max, bd.num_points).total
            )
            diff = bd.flux - model_flux
            chi2 += float(np.sum(bd.weights * (diff / bd.err) ** 2))

        return chi2

    # --- Parameter Validation ---

    def validate_parameters(self, param_defs: Sequence[ParamDef]) -> None:
        """Validate parameter definitions against the current configuration."""
        # Skip validation when using custom factories (user manages their own params)
        if self._custom_jet or self._custom_medium:
            return

        param_names = {pd.name for pd in param_defs}
        missing, incompatible = [], []

        def add_violations(required: set, forbidden: set, context: str):
            missing.extend(f"{p} ({context})" for p in required - param_names)
            incompatible.extend(
                f"{p} (not used with {context})" for p in forbidden & param_names
            )

        for rules, config_attr, label in [
            (MEDIUM_RULES, self.medium, "medium"),
            (JET_RULES, self.jet, "jet"),
        ]:
            if config_attr in rules:
                required, forbidden = rules[config_attr]
                add_violations(required, forbidden, f"{config_attr} {label}")

        for toggle, (required_on, forbidden_off) in TOGGLE_RULES.items():
            enabled = toggle == "forward_shock" or getattr(self, toggle, False)
            if enabled:
                add_violations(required_on, set(), toggle.replace("_", " "))
            else:
                incompatible.extend(
                    f"{p} ({toggle} disabled)" for p in forbidden_off & param_names
                )

        if missing or incompatible:
            msg = "Parameter validation failed:\n"
            if missing:
                msg += "Missing:\n  - " + "\n  - ".join(missing) + "\n"
            if incompatible:
                msg += "Incompatible:\n  - " + "\n  - ".join(incompatible) + "\n"
            msg += f"\nConfig: medium='{self.medium}', jet='{self.jet}', "
            msg += f"rvs_shock={self.rvs_shock}, magnetar={self.magnetar}"
            raise ValueError(msg)

    # --- Sampler Parameter Setup ---

    def _build_sampler_params(
        self, param_defs: List[ParamDef]
    ) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray, int]:
        """Build parameter labels and bounds for sampler."""
        if not (self._custom_jet or self._custom_medium):
            p_test = ModelParams()
            for pd in param_defs:
                if not hasattr(p_test, pd.name):
                    raise AttributeError(f"'{pd.name}' is not a valid MCMC parameter")

        labels, lowers, uppers = zip(
            *(
                (
                    f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name,
                    np.log10(pd.lower) if pd.scale is Scale.LOG else pd.lower,
                    np.log10(pd.upper) if pd.scale is Scale.LOG else pd.upper,
                )
                for pd in param_defs
                if pd.scale is not Scale.FIXED
            )
        )
        return labels, np.array(lowers), np.array(uppers), len(labels)

    def _build_prior_dict(
        self,
        labels: Tuple[str, ...],
        pl: np.ndarray,
        pu: np.ndarray,
        defs: List[ParamDef],
        user_priors: Optional[dict],
    ):
        """Build bilby PriorDict: user priors where provided, Uniform for the rest."""
        label_to_def = {
            (f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name): pd
            for pd in defs
            if pd.scale is not Scale.FIXED
        }

        priors_dict = {}
        for i, name in enumerate(labels):
            if user_priors and name in user_priors:
                priors_dict[name] = user_priors[name]
            else:
                priors_dict[name] = bilby.core.prior.Uniform(
                    pl[i], pu[i], name, _get_latex_label(label_to_def[name])
                )
        return bilby.core.prior.PriorDict(priors_dict)

    def _generate_initial_positions(
        self,
        param_defs: List[ParamDef],
        lower_bounds: np.ndarray,
        upper_bounds: np.ndarray,
        nwalkers: int,
        ndim: int,
    ) -> np.ndarray:
        """Generate initial walker positions centered around parameter initial values."""
        initial_vals = []
        idx = 0
        for pd in param_defs:
            if pd.scale is Scale.FIXED:
                continue
            if pd.initial is not None:
                val = np.log10(pd.initial) if pd.scale is Scale.LOG else pd.initial
            else:
                val = 0.5 * (lower_bounds[idx] + upper_bounds[idx])
            initial_vals.append(val)
            idx += 1

        initial_vals = np.array(initial_vals)
        spread = 0.1 * (upper_bounds - lower_bounds)
        pos0 = initial_vals + spread * np.random.randn(nwalkers, ndim)

        eps = 1e-6 * (upper_bounds - lower_bounds)
        pos0 = np.clip(pos0, lower_bounds + eps, upper_bounds - eps)
        return pos0

    # --- Main Fit Interface ---

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
            resolution: Optional (phi, theta, t) override. If None, uses the
                values set at construction time.
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
        self.validate_parameters(param_defs)
        defs = list(param_defs)

        # Override resolution if provided
        if resolution is not None:
            self.phi_resol, self.theta_resol, self.t_resol = resolution

        # Consolidate observation data
        self._consolidate_data()

        if len(self._all_t) == 0 and not self._band_obs:
            raise ValueError(
                "No observation data. Use add_flux_density(), add_spectrum(), or add_flux()."
            )

        # Default likelihood transform
        if log_likelihood_fn is None:

            def log_likelihood_fn(chi2):
                return -0.5 * chi2

        # Build sampler params
        labels, pl, pu, ndim = self._build_sampler_params(defs)

        self._to_params = _build_transformer(defs)

        # Build unified prior dict (used by both emcee and bilby)
        prior_dict = self._build_prior_dict(labels, pl, pu, defs, priors)

        if sampler.lower() == "emcee":
            return self._fit_emcee(
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
            return self._fit_bilby(
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

    # --- Emcee (ThreadPoolExecutor) ---

    def _fit_emcee(
        self,
        defs: List[ParamDef],
        labels: Tuple[str, ...],
        pl: np.ndarray,
        pu: np.ndarray,
        ndim: int,
        top_k: int,
        prior_dict,
        log_likelihood_fn: Callable,
        npool: Optional[int],
        **sampler_kwargs,
    ) -> FitResult:
        """Run emcee with ThreadPoolExecutor parallelism."""
        defaults = dict(SAMPLER_DEFAULTS.get("emcee", {}))
        transformer = self._to_params

        nwalkers = sampler_kwargs.pop("nwalkers", None)
        if nwalkers is None:
            nwalkers = get_optimal_nwalkers(ndim)

        if npool is None:
            npool = os.cpu_count() or 1

        logger.info(
            "Using emcee: nwalkers=%d, ndim=%d, npool=%d", nwalkers, ndim, npool
        )

        pool = ThreadPoolExecutor(max_workers=npool)

        def eval_one(theta):
            """Evaluate log-probability for a single walker."""
            params = transformer(theta)
            chi2 = self._evaluate(params)
            if not np.isfinite(chi2):
                return -np.inf
            return log_likelihood_fn(chi2)

        def log_prob_batch(samples: np.ndarray) -> np.ndarray:
            in_bounds = np.all((samples >= pl) & (samples <= pu), axis=1)
            log_probs = np.full(samples.shape[0], -np.inf)

            valid_indices = np.where(in_bounds)[0]
            if len(valid_indices) > 0:
                valid_samples = [samples[i] for i in valid_indices]
                results = list(pool.map(eval_one, valid_samples))
                log_likes = np.array(results)
                log_likes[~np.isfinite(log_likes)] = -np.inf

                valid_array = np.array(valid_samples)
                log_prior = np.zeros(len(valid_samples))
                for i, name in enumerate(labels):
                    log_prior += prior_dict[name].ln_prob(valid_array[:, i])
                log_likes += log_prior

                log_probs[valid_indices] = log_likes

            return log_probs

        pos0 = self._generate_initial_positions(defs, pl, pu, nwalkers, ndim)

        nsteps = sampler_kwargs.pop("nsteps", defaults.get("nsteps", 5000))
        nburn = sampler_kwargs.pop("nburn", defaults.get("nburn", 1000))
        thin = sampler_kwargs.pop("thin", defaults.get("thin", 1))
        moves = sampler_kwargs.pop("moves", defaults.get("moves", None))

        sampler_obj = emcee.EnsembleSampler(
            nwalkers,
            ndim,
            log_prob_batch,
            vectorize=True,
            moves=moves,
        )

        logger.info("Running emcee: nsteps=%d, nburn=%d", nsteps, nburn)

        try:
            sampler_obj.run_mcmc(pos0, nsteps, progress=True)
        finally:
            pool.shutdown(wait=False)

        chain = sampler_obj.get_chain(discard=nburn, thin=thin, flat=True)
        log_probs_flat = sampler_obj.get_log_prob(discard=nburn, thin=thin, flat=True)

        return self._process_samples(chain, log_probs_flat, labels, defs, ndim, top_k)

    # --- Bilby (dynesty, etc.) ---

    def _fit_bilby(
        self,
        defs: List[ParamDef],
        labels: Tuple[str, ...],
        ndim: int,
        sampler: str,
        npool: Optional[int],
        top_k: int,
        outdir: str,
        label: str,
        clean: bool,
        resume: bool,
        log_likelihood_fn: Callable,
        prior_dict,
        **sampler_kwargs,
    ) -> FitResult:
        """Run bilby sampler (dynesty, etc.)."""
        likelihood = AfterglowLikelihood(
            fitter=self,
            param_defs=defs,
            log_likelihood_fn=log_likelihood_fn,
            transformer=self._to_params,
        )

        if npool is None:
            npool = os.cpu_count() or 1
        else:
            npool = max(1, npool)

        pool = ThreadPoolWithClose(max_workers=npool)

        run_kwargs = {
            "likelihood": likelihood,
            "priors": prior_dict,
            "sampler": sampler,
            "outdir": outdir,
            "label": label,
            "clean": clean,
            "resume": resume,
            "pool": pool,
        }

        logger.info(
            "Running %s sampler (npool=%d, using threads)",
            sampler,
            npool,
        )

        defaults = dict(SAMPLER_DEFAULTS.get(sampler.lower(), {}))
        run_kwargs.update({**defaults, **sampler_kwargs})

        if sampler.lower() == "dynesty":
            run_kwargs["use_pool"] = {"loglikelihood": True}
            run_kwargs["queue_size"] = get_optimal_queue_size(
                npool, run_kwargs.get("nlive", 500)
            )

        try:
            result = bilby.run_sampler(**run_kwargs)
        finally:
            pool.shutdown(wait=True)

        samples = result.posterior[list(labels)].values
        log_likelihoods = result.posterior["log_likelihood"].values

        fit_result = self._process_samples(
            samples, log_likelihoods, labels, defs, ndim, top_k
        )
        fit_result.bilby_result = result
        return fit_result

    # --- Sample Processing ---

    def _process_samples(
        self,
        samples: np.ndarray,
        log_likelihoods: np.ndarray,
        labels: Tuple[str, ...],
        defs: List[ParamDef],
        ndim: int,
        top_k: int,
    ) -> FitResult:
        """Process MCMC samples to find top-k fits."""
        latex_labels = [
            _get_latex_label(pd) for pd in defs if pd.scale is not Scale.FIXED
        ]

        medians = np.median(samples, axis=0)
        stds = np.std(samples, axis=0)
        within_1sigma = np.all(np.abs(samples - medians) <= stds, axis=1)

        if np.sum(within_1sigma) >= top_k:
            candidate_samples = samples[within_1sigma]
            candidate_logp = log_likelihoods[within_1sigma]
        else:
            logger.warning(
                "Only %d samples within 1-sigma, using all samples for top-k selection",
                np.sum(within_1sigma),
            )
            candidate_samples = samples
            candidate_logp = log_likelihoods

        sorted_idx = np.argsort(candidate_logp)[::-1]
        sorted_samples = candidate_samples[sorted_idx]
        sorted_logp = candidate_logp[sorted_idx]

        seen = set()
        unique_indices = []
        for i, sample in enumerate(np.round(sorted_samples, 12)):
            key = tuple(sample)
            if key not in seen:
                seen.add(key)
                unique_indices.append(i)
                if len(unique_indices) >= top_k:
                    break

        top_k_params = sorted_samples[unique_indices]
        top_k_log_probs = sorted_logp[unique_indices]

        logger.info(
            "Found %d unique fits within 1-sigma (log-L: %.2f to %.2f)",
            len(unique_indices),
            top_k_log_probs[0],
            top_k_log_probs[-1],
        )

        return FitResult(
            samples=samples.reshape(-1, 1, ndim),
            log_probs=log_likelihoods.reshape(-1, 1),
            labels=labels,
            latex_labels=latex_labels,
            top_k_params=top_k_params,
            top_k_log_probs=top_k_log_probs,
            bilby_result=None,
        )

    # --- Post-fit Model Predictions ---

    @contextmanager
    def _override_resolution(self, resolution):
        """Temporarily override resolution, restoring on exit."""
        if resolution is not None:
            saved = (self.phi_resol, self.theta_resol, self.t_resol)
            self.phi_resol, self.theta_resol, self.t_resol = resolution
            try:
                yield
            finally:
                self.phi_resol, self.theta_resol, self.t_resol = saved
        else:
            yield

    def _require_fitted(self):
        if self._to_params is None:
            raise RuntimeError("Call .fit(...) first")

    def flux_density_grid(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu: np.ndarray,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> np.ndarray:
        """Compute flux density grid at best-fit parameters.

        Args:
            best_params: Best-fit parameter array (in sampler space)
            t: Time array [seconds]
            nu: Frequency array [Hz]
            resolution: Optional resolution override (phi, theta, t)
        """
        self._require_fitted()
        with self._override_resolution(resolution):
            params = self._to_params(best_params)
            model = self._build_model(params)
            return model.flux_density_grid(t, nu)

    def flux(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu_min: float,
        nu_max: float,
        num_points: int,
        resolution: Optional[Tuple[float, float, float]] = None,
    ) -> np.ndarray:
        """Compute integrated flux at best-fit parameters.

        Args:
            best_params: Best-fit parameter array (in sampler space)
            t: Time array [seconds]
            nu_min: Lower frequency bound [Hz]
            nu_max: Upper frequency bound [Hz]
            num_points: Number of frequency sampling points
            resolution: Optional resolution override (phi, theta, t)
        """
        self._require_fitted()
        with self._override_resolution(resolution):
            params = self._to_params(best_params)
            model = self._build_model(params)
            return model.flux(t, nu_min, nu_max, num_points)
