"""Utility functions and classes for afterglow model fitting."""

import math
import os
import types
from concurrent.futures import ThreadPoolExecutor
from typing import List, Optional

import bilby
import numpy as np
from bilby.core.sampler.emcee import Emcee as _BilbyEmcee

from ._fitting_config import _JET_CONSTRUCTORS, LATEX_LABELS
from .types import ModelParams, ParamDef, Scale
from .VegasAfterglowC import ISM, Magnetar, Wind

# Patch bilby to accept 'moves' kwarg for emcee sampler
_BilbyEmcee.default_kwargs["moves"] = None


class ThreadPoolWithClose(ThreadPoolExecutor):
    def close(self):
        self.shutdown(wait=True)

    def join(self):
        pass


# --- Default Jet/Medium Factories ---


def _default_jet_factory(fitter):
    """Build a jet factory from fitter's jet type string."""

    def factory(params: ModelParams):
        jet_type = fitter.jet
        if jet_type not in _JET_CONSTRUCTORS:
            raise ValueError(f"Unknown jet type: {jet_type}")

        cls, param_names = _JET_CONSTRUCTORS[jet_type]
        kwargs = {name: getattr(params, name) for name in param_names}

        if jet_type == "uniform":
            kwargs["theta_c"] = np.pi / 2

        kwargs["spreading"] = False
        kwargs["duration"] = params.tau if hasattr(params, "tau") else 1.0

        if jet_type != "powerlaw_wing":
            kwargs["magnetar"] = (
                Magnetar(L0=params.L0, t0=params.t0, q=params.q)
                if fitter.magnetar
                else None
            )

        return cls(**kwargs)

    return factory


def _default_medium_factory(fitter):
    """Build a medium factory from fitter's medium type string."""

    def factory(params: ModelParams):
        if fitter.medium == "ism":
            return ISM(n_ism=params.n_ism)
        elif fitter.medium == "wind":
            return Wind(
                A_star=params.A_star, n_ism=params.n_ism, n0=params.n0, k_m=params.k_m
            )
        else:
            raise ValueError(f"Unknown medium type: {fitter.medium}")

    return factory


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
    return vars(ModelParams())


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
        fitter,
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
