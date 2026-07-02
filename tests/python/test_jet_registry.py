"""The JETS/MEDIA registry: derived validation rules and constructibility.

The required/forbidden sets were historically hand-maintained in separate
JET_RULES/MEDIUM_RULES dicts; they are now derived from each JetSpec. These
tests pin the derivation to the historical values so registry edits that
would change validation behavior are caught.
"""
import math

import pytest

from VegasAfterglow.fitting.config import JETS, MEDIA
from VegasAfterglow.types import ModelParams

# Historical JET_RULES (required, forbidden), verbatim from the pre-registry config.
HISTORICAL_JET_RULES = {
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

HISTORICAL_MEDIUM_RULES = {
    "ism": ({"n_ism"}, {"A_star", "n0", "k_m"}),
    "wind": ({"A_star"}, set()),
}


@pytest.mark.parametrize("jet_type", sorted(HISTORICAL_JET_RULES))
def test_derived_jet_rules_match_historical(jet_type):
    required, forbidden = HISTORICAL_JET_RULES[jet_type]
    spec = JETS[jet_type]
    assert spec.required == required
    assert spec.forbidden == forbidden


def test_uniform_jet_rules():
    # "uniform" had no JET_RULES entry historically (it went unvalidated);
    # the registry now gives it sensible derived rules.
    spec = JETS["uniform"]
    assert spec.required == {"E_iso", "Gamma0"}
    assert spec.forbidden == {"k_e", "k_g", "E_iso_w", "Gamma0_w"}
    assert spec.fixed_kwargs == {"theta_c": math.pi / 2}


@pytest.mark.parametrize("medium_type", sorted(HISTORICAL_MEDIUM_RULES))
def test_medium_rules_match_historical(medium_type):
    required, forbidden = HISTORICAL_MEDIUM_RULES[medium_type]
    spec = MEDIA[medium_type]
    assert spec.required == required
    assert spec.forbidden == forbidden


@pytest.mark.parametrize("jet_type", sorted(JETS))
def test_every_jet_constructs_from_defaults(jet_type):
    spec = JETS[jet_type]
    params = ModelParams()
    kwargs = {name: getattr(params, name) for name in spec.params}
    kwargs.update(spec.fixed_kwargs)
    jet = spec.constructor(**kwargs)
    assert jet is not None


@pytest.mark.parametrize("medium_type", sorted(MEDIA))
def test_every_medium_constructs_from_defaults(medium_type):
    spec = MEDIA[medium_type]
    params = ModelParams()
    if medium_type == "ism":
        params.n_ism = 1.0  # default 0 is rejected by ISM validation
    if medium_type == "wind":
        params.A_star = 0.1  # default 0 is rejected by Wind validation
    kwargs = {name: getattr(params, name) for name in spec.params}
    medium = spec.constructor(**kwargs)
    assert medium is not None


@pytest.mark.parametrize("jet_type", sorted(JETS))
def test_registry_self_consistency(jet_type):
    spec = JETS[jet_type]
    assert not (spec.required & spec.forbidden)
    assert not (set(spec.params) & set(spec.fixed_kwargs)), (
        "a constructor arg must be either sampled or fixed, not both"
    )
