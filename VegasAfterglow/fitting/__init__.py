"""Afterglow model fitting: data setup, sampler orchestration, and post-fit prediction.

Public API:
    Fitter              — model-based fitter (one class, in :mod:`.fitter`)
    AfterglowLikelihood — bilby-compatible likelihood object (in :mod:`.utils`)

Internal helpers live in :mod:`.config` (jet/medium rules, parameter latex
labels) and :mod:`.utils` (transformer, factories, optimal worker counts).
"""

from .fitter import Fitter
from .utils import AfterglowLikelihood

__all__ = ["Fitter", "AfterglowLikelihood"]
