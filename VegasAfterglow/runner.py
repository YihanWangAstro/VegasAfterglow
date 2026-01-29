"""Afterglow model fitting with bilby samplers."""

import logging
import math
import os
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional, Sequence, Tuple, Type

import bilby
import emcee
import numpy as np
from bilby.core.sampler.emcee import Emcee as _BilbyEmcee

from .types import FitResult, ModelParams, ObsData, ParamDef, Scale, Setups, VegasMC

# Import C++ classes for MCMC
try:
    from .VegasAfterglowC import ParamDef as CppParamDef
    from .VegasAfterglowC import ParamTransformer as CppParamTransformer
except ImportError:
    CppParamDef = None
    CppParamTransformer = None

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
        # "sample": "rwalk",
        # "walks": 25,
        "sample": "rslice",
        # "bootstrap": 5,
        # "slices": 3,
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


def get_optimal_nwalkers(ndim: int, ncpu: Optional[int] = None) -> int:
    """Compute optimal nwalkers for emcee with OpenMP batch processing.

    The formula ensures:
    1. nwalkers >= 4 * ndim (minimum for emcee algorithm)
    2. nwalkers is a multiple of 4 * ncpu (emcee splits walkers in half,
       so the half-batch should fill all CPUs evenly)

    Args:
        ndim: Number of dimensions (sampled parameters)
        ncpu: Number of CPU cores. If None, auto-detect.

    Returns:
        Optimal number of walkers
    """
    if ncpu is None:
        ncpu = os.cpu_count() or 1

    math_floor = 4 * ndim
    align_unit = 2 * ncpu
    units_needed = math.ceil(math_floor / align_unit)
    units_needed = max(1, units_needed)

    return units_needed * align_unit


def get_optimal_queue_size(ncpu, nlive):
    """
    Calculates the optimal queue_size for dynesty based on hardware and sampling parameters.
    Parameters:
    -----------
    ncpu : int
        Number of CPUs or Threads available.
    nlive : int
        Number of live points used in the Nested Sampler.

    Returns:
    --------
    int
        The recommended queue_size.
    """
    target_multiplier = 2
    target_size = ncpu * target_multiplier

    max_safe_size = int(nlive * 0.15)
    optimal_size = min(target_size, max_safe_size)
    aligned_size = (optimal_size // ncpu) * ncpu

    if aligned_size < ncpu:
        aligned_size = ncpu
    return aligned_size


class AfterglowLikelihood(bilby.Likelihood):
    """Bilby-compatible likelihood using C++ ParamTransformer for efficiency."""

    __slots__ = (
        "parameters",
        "param_keys",
        "_data",
        "_config",
        "_model_cls",
        "_model",
        "_theta",
        "_transformer",
    )

    def __init__(
        self,
        data: ObsData,
        config: Setups,
        param_defs: List[ParamDef],
        model_cls: Type[VegasMC],
    ):
        if CppParamTransformer is None:
            raise RuntimeError(
                "VegasAfterglowC extension not loaded. "
                "AfterglowLikelihood requires the C++ ParamTransformer."
            )

        # Build parameter keys from defs (only non-fixed params)
        param_keys = tuple(
            f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name
            for pd in param_defs
            if pd.scale is not Scale.FIXED
        )
        super().__init__(parameters={key: None for key in param_keys})
        self.param_keys = param_keys
        self._data, self._config, self._model_cls = data, config, model_cls
        self._model = None
        self._theta = np.empty(len(param_keys), dtype=np.float64)

        # Create C++ ParamTransformer (handles all transformation logic in C++)
        cpp_param_defs = [Fitter._to_cpp_param_def(pd) for pd in param_defs]
        self._transformer = CppParamTransformer(cpp_param_defs)

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
        # Gather parameter values into array
        for i, key in enumerate(self.param_keys):
            self._theta[i] = self.parameters[key]

        BAD_LOGL_FLOOR = -1e100

        try:
            # Use batch interface with single sample - C++ handles all transformations
            samples = self._theta.reshape(1, -1)
            chi2_arr = self._get_model().batch_estimate_chi2(samples, self._transformer)
            chi2 = chi2_arr[0]
            loglike = -0.5 * chi2 if np.isfinite(chi2) else BAD_LOGL_FLOOR
            # print(f"[DEBUG] log_likelihood call: params={self._theta}, chi2={chi2}, loglike={loglike}")
            return loglike
        except Exception as e:
            print(f"[DEBUG] log_likelihood exception: params={self._theta}, error={e}")
            return BAD_LOGL_FLOOR


class Fitter:
    def __init__(self, data: ObsData, config: Setups):
        self.data = data
        self.config = config
        self._param_defs = None
        self._to_params = None

    @staticmethod
    def _to_cpp_param_def(pd: ParamDef):
        """Convert Python ParamDef to C++ ParamDef.

        Args:
            pd: Python parameter definition

        Returns:
            C++ parameter definition with proper scale encoding and initial values
        """
        cpp_pd = CppParamDef()
        cpp_pd.name = pd.name
        cpp_pd.lower = pd.lower
        cpp_pd.upper = pd.upper
        # Convert enum to int: LINEAR=0, LOG=1, FIXED=2
        if pd.scale == Scale.LINEAR:
            cpp_pd.scale = 0
        elif pd.scale == Scale.LOG:
            cpp_pd.scale = 1
        elif pd.scale == Scale.FIXED:
            cpp_pd.scale = 2
        else:
            raise ValueError(f"Unknown scale: {pd.scale}")

        # Set initial value (optional for sampled parameters)
        if pd.initial is not None:
            cpp_pd.initial = pd.initial
        else:
            # Compute sensible default based on scale
            if pd.scale == Scale.LOG:
                # Geometric mean for log-scale parameters (centered in log space)
                cpp_pd.initial = np.sqrt(pd.lower * pd.upper)
            else:
                # Arithmetic mean for linear-scale and fixed parameters
                cpp_pd.initial = 0.5 * (pd.lower + pd.upper)
        return cpp_pd

    def validate_parameters(self, param_defs: Sequence[ParamDef]) -> None:
        """Validate parameter definitions against the current configuration."""
        param_names = {pd.name for pd in param_defs}
        missing, incompatible = [], []

        def add_violations(required: set, forbidden: set, context: str):
            """Add missing and incompatible parameters for given context."""
            missing.extend(f"{p} ({context})" for p in required - param_names)
            incompatible.extend(
                f"{p} (not used with {context})" for p in forbidden & param_names
            )

        # Check medium and jet configuration rules
        for rules, config_attr, label in [
            (MEDIUM_RULES, self.config.medium, "medium"),
            (JET_RULES, self.config.jet, "jet"),
        ]:
            if config_attr in rules:
                required, forbidden = rules[config_attr]
                add_violations(required, forbidden, f"{config_attr} {label}")

        # Check physics toggle rules (forward_shock, rvs_shock, magnetar)
        for toggle, (required_on, forbidden_off) in TOGGLE_RULES.items():
            enabled = toggle == "forward_shock" or getattr(self.config, toggle, False)
            if enabled:
                add_violations(required_on, set(), toggle.replace("_", " "))
            else:
                incompatible.extend(
                    f"{p} ({toggle} disabled)" for p in forbidden_off & param_names
                )

        # Raise error if any violations found
        if missing or incompatible:
            self._raise_validation_error(missing, incompatible)

    def _raise_validation_error(
        self, missing: List[str], incompatible: List[str]
    ) -> None:
        """Format and raise validation error with diagnostic information."""
        msg = "Parameter validation failed:\n"
        if missing:
            msg += "Missing:\n  - " + "\n  - ".join(missing) + "\n"
        if incompatible:
            msg += "Incompatible:\n  - " + "\n  - ".join(incompatible) + "\n"
        msg += f"\nConfig: medium='{self.config.medium}', jet='{self.config.jet}', "
        msg += f"rvs_shock={self.config.rvs_shock}, magnetar={self.config.magnetar}"
        raise ValueError(msg)

    def _build_sampler_params(
        self, param_defs: List[ParamDef]
    ) -> Tuple[Tuple[str, ...], np.ndarray, np.ndarray, int]:
        """Build parameter labels and bounds for sampler from definitions.

        Returns:
            labels: Parameter names (with log10_ prefix for log-scale params)
            lower_bounds: Lower bounds array (in log10 space for log-scale params)
            upper_bounds: Upper bounds array (in log10 space for log-scale params)
            ndim: Number of free (sampled) parameters
        """
        # Validate parameter names exist in ModelParams
        p_test = ModelParams()
        for pd in param_defs:
            if not hasattr(p_test, pd.name):
                raise AttributeError(f"'{pd.name}' is not a valid MCMC parameter")

        # Extract labels and bounds for non-fixed parameters
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

    def _generate_initial_positions(
        self,
        param_defs: List[ParamDef],
        lower_bounds: np.ndarray,
        upper_bounds: np.ndarray,
        nwalkers: int,
        ndim: int,
    ) -> np.ndarray:
        """Generate initial walker positions centered around parameter initial values.

        Args:
            param_defs: Parameter definitions with initial values
            lower_bounds: Lower bounds in sampler space (log10 for log-scale params)
            upper_bounds: Upper bounds in sampler space
            nwalkers: Number of walkers
            ndim: Number of dimensions

        Returns:
            Initial positions array of shape (nwalkers, ndim)
        """
        # Extract initial values for non-fixed parameters
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

        # Add small random scatter around initial values
        initial_vals = np.array(initial_vals)
        spread = 0.1 * (upper_bounds - lower_bounds)
        pos0 = initial_vals + spread * np.random.randn(nwalkers, ndim)

        # Ensure positions are within bounds
        eps = 1e-6 * (upper_bounds - lower_bounds)
        pos0 = np.clip(pos0, lower_bounds + eps, upper_bounds - eps)
        return pos0

    @staticmethod
    def _get_latex_label(param_def: ParamDef) -> str:
        """Generate LaTeX label for parameter (with log10 wrapper if needed)."""
        base_latex = LATEX_LABELS.get(param_def.name, param_def.name)
        if param_def.scale is Scale.LOG:
            return rf"$\log_{{10}}({base_latex.strip('$')})$"
        return base_latex

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
        vectorize: bool = True,
        **sampler_kwargs,
    ) -> FitResult:
        """Run sampler for parameter estimation.

        Args:
            param_defs: Parameter definitions
            resolution: (phi_resol, theta_resol, t_resol) tuple
            sampler: Sampler name ('emcee', 'dynesty', etc.)
            npool: Number of processes (ignored for vectorized emcee)
            top_k: Number of top samples to return
            outdir: Output directory for bilby
            label: Run label for bilby
            clean: Clean previous runs (bilby)
            resume: Resume from previous run (bilby)
            vectorize: Use OpenMP batch processing for emcee (default True)
            **sampler_kwargs: Additional sampler arguments (nsteps, nburn, nwalkers, etc.)

        Returns:
            FitResult with samples and best-fit parameters
        """
        self.validate_parameters(param_defs)
        defs = list(param_defs)
        self._param_defs = defs

        # Extract parameter labels and bounds for sampler
        labels, pl, pu, ndim = self._build_sampler_params(defs)

        # Create C++ ParamTransformer for flux predictions
        cpp_param_defs = [self._to_cpp_param_def(pd) for pd in defs]
        self._to_params = CppParamTransformer(cpp_param_defs)

        # Use vectorized emcee if available and requested
        sampler_lower = sampler.lower()
        use_vectorized = (
            sampler_lower in ("emcee") and vectorize and CppParamTransformer is not None
        )

        # Force npool=None for batch evaluation and warn user
        if use_vectorized:
            # Always check both positional and keyword npool
            npool_effective = npool
            if "npool" in sampler_kwargs:
                npool_effective = sampler_kwargs["npool"]
            if npool_effective not in (None, 1):
                import warnings

                warnings.warn(
                    f"npool={npool_effective} is ignored: batch likelihood uses multi-threaded parallelism. "
                    "Set npool=None or 1 for vectorized emcee"
                )
            # Always force npool=None for batch
            sampler_kwargs["npool"] = None
            return self._fit_vectorized_sampler(
                sampler_lower,
                defs,
                labels,
                pl,
                pu,
                ndim,
                resolution,
                top_k,
                **sampler_kwargs,
            )

        # Fall back to bilby for other samplers or non-vectorized emcee
        return self._fit_bilby(
            defs,
            labels,
            pl,
            pu,
            ndim,
            resolution,
            sampler,
            npool,
            top_k,
            outdir,
            label,
            clean,
            resume,
            **sampler_kwargs,
        )

    def _fit_vectorized_sampler(
        self,
        sampler: str,
        defs: List[ParamDef],
        labels: Tuple[str, ...],
        pl: np.ndarray,
        pu: np.ndarray,
        ndim: int,
        resolution: Tuple[float, float, float],
        top_k: int,
        **sampler_kwargs,
    ) -> FitResult:
        """Run emcee with OpenMP-vectorized batch chi2 computation."""
        defaults = dict(SAMPLER_DEFAULTS.get(sampler, {}))

        # Reuse C++ ParamTransformer created in fit()
        transformer = self._to_params

        # Compute optimal nwalkers for CPU utilization
        nwalkers = sampler_kwargs.pop("nwalkers", None)
        if nwalkers is None:
            nwalkers = get_optimal_nwalkers(ndim)
        logger.info(f"Using vectorized {sampler}: nwalkers=%d, ndim=%d", nwalkers, ndim)

        # Create model (used by all OpenMP threads internally)
        model = VegasMC(self.data)
        model.set(clone_config(self.config, resolution))

        # Vectorized log-probability function (clean and simple!)
        def log_prob_batch(samples: np.ndarray) -> np.ndarray:
            in_bounds = np.all((samples >= pl) & (samples <= pu), axis=1)
            log_probs = np.full(samples.shape[0], -np.inf)
            if np.any(in_bounds):
                chi2_arr = model.batch_estimate_chi2(samples[in_bounds], transformer)
                log_likes = -0.5 * chi2_arr
                log_likes[~np.isfinite(log_likes)] = -np.inf
                log_probs[in_bounds] = log_likes
            return log_probs

        # Generate initial walker positions
        pos0 = self._generate_initial_positions(defs, pl, pu, nwalkers, ndim)

        # Get sampler parameters
        nsteps = sampler_kwargs.pop("nsteps", defaults.get("nsteps", 5000))
        nburn = sampler_kwargs.pop("nburn", defaults.get("nburn", 1000))
        thin = sampler_kwargs.pop("thin", defaults.get("thin", 1))
        moves = sampler_kwargs.pop("moves", defaults.get("moves", None))

        # Create sampler with vectorized log_prob
        if sampler == "emcee":
            sampler_obj = emcee.EnsembleSampler(
                nwalkers, ndim, log_prob_batch, vectorize=True, moves=moves
            )
        else:
            raise ValueError(f"Unknown vectorized sampler: {sampler}")

        logger.info(f"Running {sampler}: nsteps=%d, nburn=%d", nsteps, nburn)
        sampler_obj.run_mcmc(pos0, nsteps, progress=True)

        # Extract samples
        chain = sampler_obj.get_chain(discard=nburn, thin=thin, flat=True)
        log_probs_flat = sampler_obj.get_log_prob(discard=nburn, thin=thin, flat=True)

        return self._process_samples(chain, log_probs_flat, labels, defs, ndim, top_k)

    def _fit_bilby(
        self,
        defs: List[ParamDef],
        labels: Tuple[str, ...],
        pl: np.ndarray,
        pu: np.ndarray,
        ndim: int,
        resolution: Tuple[float, float, float],
        sampler: str,
        npool: int,
        top_k: int,
        outdir: str,
        label: str,
        clean: bool,
        resume: bool,
        **sampler_kwargs,
    ) -> FitResult:
        """Run bilby sampler (non-vectorized path)."""
        # Build mapping from label to parameter definition
        label_to_def = {
            (f"log10_{pd.name}" if pd.scale is Scale.LOG else pd.name): pd
            for pd in defs
            if pd.scale is not Scale.FIXED
        }

        likelihood = AfterglowLikelihood(
            data=self.data,
            config=clone_config(self.config, resolution),
            param_defs=defs,
            model_cls=VegasMC,
        )

        priors = bilby.core.prior.PriorDict(
            {
                name: bilby.core.prior.Uniform(
                    pl[i], pu[i], name, self._get_latex_label(label_to_def[name])
                )
                for i, name in enumerate(labels)
            }
        )

        run_kwargs = {
            "likelihood": likelihood,
            "priors": priors,
            "sampler": sampler,
            "outdir": outdir,
            "label": label,
            "clean": clean,
            "resume": resume,
        }

        # Use ThreadPoolExecutor instead of multiprocessing for better performance
        # (GIL is released in C++ code via py::gil_scoped_release)
        pool = None
        if npool is None:
            npool = os.cpu_count() or 1
        else:
            npool = max(1, npool)

        pool = ThreadPoolWithClose(max_workers=npool)
        run_kwargs["pool"] = pool
        logger.info(
            "Running %s sampler at resolution %s (npool=%d, using threads)",
            sampler,
            resolution,
            npool,
        )

        defaults = dict(SAMPLER_DEFAULTS.get(sampler.lower(), {}))
        run_kwargs.update({**defaults, **sampler_kwargs})

        # Explicitly tell dynesty to use pool for likelihood evaluations
        if sampler.lower() == "dynesty":
            run_kwargs["use_pool"] = {"loglikelihood": True}
            run_kwargs["queue_size"] = get_optimal_queue_size(
                npool, run_kwargs.get("nlive", 500)
            )

        try:
            result = bilby.run_sampler(**run_kwargs)
        finally:
            # Clean up thread pool
            if pool is not None:
                pool.shutdown(wait=True)

        samples = result.posterior[list(labels)].values
        log_likelihoods = result.posterior["log_likelihood"].values

        fit_result = self._process_samples(
            samples, log_likelihoods, labels, defs, ndim, top_k
        )
        # Add bilby result
        return FitResult(
            samples=fit_result.samples,
            log_probs=fit_result.log_probs,
            labels=fit_result.labels,
            latex_labels=fit_result.latex_labels,
            top_k_params=fit_result.top_k_params,
            top_k_log_probs=fit_result.top_k_log_probs,
            bilby_result=result,
        )

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
            self._get_latex_label(pd) for pd in defs if pd.scale is not Scale.FIXED
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
        """Compute flux density grid at best-fit parameters."""
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
        """Compute integrated flux at best-fit parameters."""
        model, p = self._prepare_model(best_params, resolution)
        return model.flux(p, t, nu_min, nu_max, num_points)
