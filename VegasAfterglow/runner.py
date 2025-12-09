# src/vegasglow/runner.py
"""Afterglow model fitting with bilby samplers."""

import logging
from typing import Callable, List, Sequence, Tuple, Type

import bilby
import emcee
import numpy as np

from .types import FitResult, ModelParams, ObsData, ParamDef, Scale, Setups, VegasMC

logger = logging.getLogger(__name__)

# =============================================================================
# Configuration Constants
# =============================================================================

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

# Validation rules: {key: (required_params, incompatible_params)}
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

# Toggle rules: {(config_attr, enabled_value): (required, incompatible)}
TOGGLE_RULES = {
    "forward_shock": ({"eps_e", "eps_B", "p"}, set()),  # Always required
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
        "sample": "rwalk",
        "walks": 100,
        "nact": 10,
        "maxmcmc": 5000,
    },
    "emcee": {
        "nsteps": 5000,
        "nburn": 1000,
        "thin": 1,
        "moves": [
            (emcee.moves.StretchMove(), 0.5),
            (emcee.moves.DEMove(), 0.4),
            (emcee.moves.DESnookerMove(), 0.1),
        ],
    },
}


# =============================================================================
# Utility Functions
# =============================================================================


def clone_config(base_cfg: Setups, resolution: Tuple[float, float, float]) -> Setups:
    """Clone config and override resolution (phi, theta, t)."""
    cfg = type(base_cfg)()
    for attr in dir(base_cfg):
        if not attr.startswith("_") and hasattr(cfg, attr):
            try:
                setattr(cfg, attr, getattr(base_cfg, attr))
            except Exception:
                pass
    cfg.phi_resol, cfg.theta_resol, cfg.t_resol = resolution
    return cfg


# =============================================================================
# Core Classes
# =============================================================================


class AfterglowLikelihood(bilby.Likelihood):
    """Bilby-compatible likelihood with lazy model creation for multiprocessing."""

    __slots__ = (
        "parameters",
        "to_params",
        "param_keys",
        "_data",
        "_config",
        "_model_cls",
        "_model",
        "_theta",
    )

    def __init__(
        self,
        data: ObsData,
        config: Setups,
        to_params: Callable[[np.ndarray], ModelParams],
        param_keys: Sequence[str],
        model_cls: Type[VegasMC],
    ):
        super().__init__(parameters={key: None for key in param_keys})
        self.to_params = to_params
        self.param_keys = tuple(param_keys)
        self._data, self._config, self._model_cls = data, config, model_cls
        self._model = None
        self._theta = np.empty(len(param_keys), dtype=np.float64)

    def __getstate__(self):
        return {k: getattr(self, k) for k in self.__slots__ if k != "_model"}

    def __setstate__(self, state):
        for k, v in state.items():
            setattr(self, k, v)
        self._model = None

    def _get_model(self) -> VegasMC:
        if self._model is None:
            self._model = self._model_cls(self._data)
            self._model.set(self._config)
        return self._model

    def log_likelihood(self) -> float:
        for i, key in enumerate(self.param_keys):
            self._theta[i] = self.parameters[key]
        try:
            chi2 = self._get_model().estimate_chi2(self.to_params(self._theta))
            return -0.5 * chi2 if chi2 == chi2 else -np.inf
        except Exception:
            return -np.inf


class ParamTransformer:
    """
    Picklable parameter transformer for converting sampler arrays to ModelParams.

    This class replaces a local closure to enable multiprocessing with pickle.
    Optimized with pre-computed indices and vectorized log10 conversion.
    """

    __slots__ = ("_names", "_log_mask", "_fixed_names", "_fixed_values")

    def __init__(self, param_defs: List[ParamDef]):
        names = []
        log_mask = []
        fixed_names = []
        fixed_values = []

        for pd in param_defs:
            if pd.scale is Scale.FIXED:
                fixed_names.append(pd.name)
                fixed_values.append(0.5 * (pd.lower + pd.upper))
            else:
                names.append(pd.name)
                log_mask.append(pd.scale is Scale.LOG)

        self._names = tuple(names)
        self._log_mask = np.array(log_mask, dtype=bool)
        self._fixed_names = tuple(fixed_names)
        self._fixed_values = tuple(fixed_values)

    def __call__(self, x: np.ndarray) -> ModelParams:
        p = ModelParams()
        values = np.where(self._log_mask, 10.0**x, x)
        for name, val in zip(self._names, values):
            setattr(p, name, val)
        for name, val in zip(self._fixed_names, self._fixed_values):
            setattr(p, name, val)
        return p


class Fitter:
    def __init__(self, data: ObsData, config: Setups):
        """
        Parameters
        ----------
        data : ObsData
            Observed light curves and spectra.
        config : Setups
            Model configuration (grids, environment, etc).
        """
        self.data = data
        self.config = config
        self._param_defs = None
        self._to_params = None

    def validate_parameters(self, param_defs: Sequence[ParamDef]) -> None:
        """Validate parameter definitions against the current configuration."""
        param_names = {pd.name for pd in param_defs}
        missing, incompatible = [], []

        def check(required: set, incompat: set, context: str):
            missing.extend(f"{p} ({context})" for p in required - param_names)
            incompatible.extend(
                f"{p} (not used with {context})" for p in incompat & param_names
            )

        def apply_rules(rules: dict, config_value: str, context_name: str):
            if config_value in rules:
                req, incom = rules[config_value]
                check(req, incom, f"{config_value} {context_name}")

        # Apply lookup-based rules
        apply_rules(MEDIUM_RULES, self.config.medium, "medium")
        apply_rules(JET_RULES, self.config.jet, "jet")

        # Apply toggle-based rules
        for toggle, (required_on, incompatible_off) in TOGGLE_RULES.items():
            enabled = toggle == "forward_shock" or getattr(self.config, toggle, False)
            if enabled:
                check(required_on, set(), toggle.replace("_", " "))
            else:
                incompatible.extend(
                    f"{p} ({toggle} disabled)" for p in incompatible_off & param_names
                )

        if missing or incompatible:
            msg = "Parameter validation failed:\n"
            if missing:
                msg += "Missing:\n  - " + "\n  - ".join(missing) + "\n"
            if incompatible:
                msg += "Incompatible:\n  - " + "\n  - ".join(incompatible) + "\n"
            msg += f"\nConfig: medium='{self.config.medium}', jet='{self.config.jet}', "
            msg += f"rvs_shock={self.config.rvs_shock}, magnetar={self.config.magnetar}"
            raise ValueError(msg)

    def fit(
        self,
        param_defs: Sequence[ParamDef],
        resolution: Tuple[float, float, float] = (0.3, 1, 10),
        sampler: str = "dynesty",
        npool: int = 1,
        top_k: int = 10,
        outdir: str = "bilby_output",
        label: str = "afterglow",
        clean: bool = True,
        resume: bool = False,
        **sampler_kwargs,
    ) -> FitResult:
        """
        Run bilby sampler for parameter estimation.

        Parameters
        ----------
        param_defs : Sequence[ParamDef]
            Parameter definitions for fitting.
        resolution : Tuple[float, float, float]
            Grid resolution (phi, theta, t).
        sampler : str
            'dynesty' (nested sampling) or 'emcee' (MCMC).
        npool : int
            Number of parallel processes.
        top_k : int
            Number of top fits to return.
        **sampler_kwargs
            Sampler-specific options (nlive, dlogz, nwalkers, iterations, etc.)
        """
        self.validate_parameters(param_defs)
        defs = list(param_defs)
        self._param_defs = defs

        # Extract labels and bounds (log-transform if needed)
        # For LOG scale: prepend 'log10_' to name and transform bounds to log10 space
        labels, lowers, uppers = zip(
            *(
                (
                    f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name,
                    np.log10(pd.lower) if pd.scale is Scale.LOG else pd.lower,
                    np.log10(pd.upper) if pd.scale is Scale.LOG else pd.upper,
                )
                for pd in defs
                if pd.scale is not Scale.FIXED
            )
        )
        pl, pu = np.array(lowers), np.array(uppers)
        ndim = len(labels)

        # Validate parameter names
        p_test = ModelParams()
        for pd in defs:
            if not hasattr(p_test, pd.name):
                raise AttributeError(f"'{pd.name}' is not a valid MCMC parameter")

        self._to_params = ParamTransformer(defs)

        # Create likelihood and priors
        likelihood = AfterglowLikelihood(
            data=self.data,
            config=clone_config(self.config, resolution),
            to_params=self._to_params,
            param_keys=labels,
            model_cls=VegasMC,
        )

        # Build LaTeX labels: wrap LOG-scale params in \log_{10}(...)
        def get_latex_label(param_name: str, param_def: ParamDef) -> str:
            base_latex = LATEX_LABELS.get(param_def.name, param_def.name)
            if param_def.scale is Scale.LOG:
                return rf"$\log_{{10}}({base_latex.strip('$')})$"
            return base_latex

        # Map labels back to their param defs for LaTeX generation
        label_to_def = {}
        for pd in defs:
            if pd.scale is not Scale.FIXED:
                label_name = f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name
                label_to_def[label_name] = pd

        priors = bilby.core.prior.PriorDict(
            {
                name: bilby.core.prior.Uniform(
                    pl[i], pu[i], name, get_latex_label(name, label_to_def[name])
                )
                for i, name in enumerate(labels)
            }
        )

        # Build run kwargs with defaults
        run_kwargs = {
            "likelihood": likelihood,
            "priors": priors,
            "sampler": sampler,
            "outdir": outdir,
            "label": label,
            "clean": clean,
            "resume": resume,
        }

        # Apply sampler-specific defaults and handle parallelization
        defaults = dict(SAMPLER_DEFAULTS.get(sampler.lower(), {}))
        sampler_lower = sampler.lower()
        run_kwargs["npool"] = npool

        if sampler_lower == "emcee":
            defaults.setdefault("nwalkers", 2 * ndim)

        run_kwargs.update({**defaults, **sampler_kwargs})

        logger.info(
            "Running %s sampler at resolution %s (npool=%d)", sampler, resolution, npool
        )
        result = bilby.run_sampler(**run_kwargs)

        # Extract and process results
        samples = result.posterior[list(labels)].values
        log_likelihoods = result.posterior["log_likelihood"].values

        # Find top-k unique fits
        sorted_idx = np.argsort(log_likelihoods)[::-1]
        _, unique_idx = np.unique(
            np.round(samples[sorted_idx], 12), axis=0, return_index=True
        )
        final_idx = sorted_idx[np.sort(unique_idx)[:top_k]]

        logger.info(
            "Found %d unique fits (log-L: %.2f to %.2f)",
            len(final_idx),
            log_likelihoods[final_idx[0]],
            log_likelihoods[final_idx[-1]],
        )

        # Extract LaTeX labels from priors for corner plots
        latex_labels = [priors[name].latex_label for name in labels]

        return FitResult(
            samples=samples.reshape(-1, 1, ndim),
            log_probs=log_likelihoods.reshape(-1, 1),
            labels=labels,
            latex_labels=latex_labels,
            top_k_params=samples[final_idx],
            top_k_log_probs=log_likelihoods[final_idx],
            bilby_result=result,
        )

    def _prepare_model(
        self, best_params: np.ndarray, resolution: Tuple[float, float, float]
    ):
        """Prepare model with given parameters and resolution."""
        if self._to_params is None:
            raise RuntimeError("Call .fit(...) first")
        model = VegasMC(self.data)
        model.set(clone_config(self.config, resolution))
        return model, self._to_params(best_params)

    def flux_density_grid(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu: np.ndarray,
        resolution: Tuple[float, float, float] = (0.3, 1, 10),
    ) -> np.ndarray:
        """Compute flux density grid at best-fit parameters.
        Returns shape (len(nu), len(t)).
        """
        model, p = self._prepare_model(best_params, resolution)
        return model.flux_density_grid(p, t, nu)

    def flux(
        self,
        best_params: np.ndarray,
        t: np.ndarray,
        nu_min: float,
        nu_max: float,
        num_points: int,
        resolution: Tuple[float, float, float] = (0.3, 1, 10),
    ) -> np.ndarray:
        """Compute integrated flux at best-fit parameters. Returns shape (len(t),)."""
        model, p = self._prepare_model(best_params, resolution)
        return model.flux(p, t, nu_min, nu_max, num_points)
