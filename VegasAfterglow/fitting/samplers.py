"""Sampler back-ends called from :meth:`Fitter.fit`.

* ``fit_emcee`` — emcee with a ``ThreadPoolExecutor`` and vectorised
  log-probability evaluation.
* ``fit_bilby`` — bilby front-end for dynesty (and other bilby samplers).
* ``process_samples`` — turn a flat ``(samples, log_likelihoods)`` pair into a
  ``FitResult`` with the top-K best-fit picks.
"""

import logging
import os
from concurrent.futures import ThreadPoolExecutor
from typing import Callable, List, Optional, Tuple

import bilby
import emcee
import numpy as np

from ..types import FitResult, ParamDef, Scale
from .config import SAMPLER_DEFAULTS
from .params import generate_initial_positions
from .utils import (
    AfterglowLikelihood,
    ThreadPoolWithClose,
    _get_latex_label,
    get_optimal_nwalkers,
    get_optimal_queue_size,
)

logger = logging.getLogger(__name__)


def fit_emcee(
    fitter,
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
    transformer = fitter._to_params

    nwalkers = sampler_kwargs.pop("nwalkers", None)
    if nwalkers is None:
        nwalkers = get_optimal_nwalkers(ndim)

    if npool is None:
        npool = os.cpu_count() or 1

    logger.info("Using emcee: nwalkers=%d, ndim=%d, npool=%d", nwalkers, ndim, npool)

    pool = ThreadPoolExecutor(max_workers=npool)

    def eval_one(theta):
        """Evaluate log-probability for a single walker."""
        try:
            params = transformer(theta)
            chi2 = fitter._evaluate(params)
            if not np.isfinite(chi2):
                return -np.inf
            return log_likelihood_fn(chi2)
        except Exception:
            return -np.inf

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

    pos0 = generate_initial_positions(defs, pl, pu, nwalkers, ndim)

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
        pool.shutdown(wait=True)

    chain = sampler_obj.get_chain(discard=nburn, thin=thin, flat=True)
    log_probs_flat = sampler_obj.get_log_prob(discard=nburn, thin=thin, flat=True)

    return process_samples(chain, log_probs_flat, labels, defs, ndim, top_k)


def fit_bilby(
    fitter,
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
        fitter=fitter,
        param_defs=defs,
        log_likelihood_fn=log_likelihood_fn,
        transformer=fitter._to_params,
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

    fit_result = process_samples(samples, log_likelihoods, labels, defs, ndim, top_k)
    fit_result.bilby_result = result
    return fit_result


def process_samples(
    samples: np.ndarray,
    log_likelihoods: np.ndarray,
    labels: Tuple[str, ...],
    defs: List[ParamDef],
    ndim: int,
    top_k: int,
) -> FitResult:
    """Pick top-K best-fit parameter sets and wrap everything in a ``FitResult``.

    Selection rule:
      1. Keep only samples inside the per-parameter 1-sigma band around the
         median (so the top-K reflects the bulk of the posterior, not a stray
         high-likelihood outlier). If fewer than ``top_k`` survive, fall back
         to the full set with a warning.
      2. Sort surviving samples by log-likelihood, descending.
      3. Pick the first ``top_k`` *unique* parameter vectors (rounded to 12
         decimals) so repeated identical proposals don't all count.
    """
    latex_labels = [_get_latex_label(pd) for pd in defs if pd.scale is not Scale.fixed]

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
