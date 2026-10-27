"""Bad jet/medium arguments fail at construction, and swallowed evaluation errors are reported.

A jet or medium that was neither a registered name nor a callable used to pass
validation, then fail inside every likelihood call; the sampler saw only -inf and
hung while seeding its initial points.
"""
import logging

import numpy as np
import pytest

from VegasAfterglow import ISM, Fitter, ParamDef, Scale, TophatJet
from VegasAfterglow.fitting import utils

pytestmark = pytest.mark.physics


@pytest.mark.parametrize("kwargs", [{"jet": "tophatt"}, {"medium": "wnd"}])
def test_unknown_name_is_rejected(kwargs):
    """A misspelt jet or medium name raises ValueError listing the valid names."""
    with pytest.raises(ValueError, match="choose from"):
        Fitter(z=0.5, lumi_dist=3e28, **kwargs)


@pytest.mark.parametrize(
    "kwargs",
    [{"jet": TophatJet(0.1, 1e52, 100.0)}, {"medium": ISM(1.0)}, {"jet": None}],
)
def test_non_name_non_callable_is_rejected(kwargs):
    """A built jet/medium object (or None) raises TypeError pointing at the factory form."""
    with pytest.raises(TypeError, match="callable factory"):
        Fitter(z=0.5, lumi_dist=3e28, **kwargs)


def test_registered_names_and_factories_still_construct():
    """Every registered name and a callable factory remain valid."""
    from VegasAfterglow.fitting.config import JETS, MEDIA

    for name in JETS:
        Fitter(z=0.5, lumi_dist=3e28, jet=name)
    for name in MEDIA:
        Fitter(z=0.5, lumi_dist=3e28, medium=name)
    Fitter(z=0.5, lumi_dist=3e28, jet=lambda p: TophatJet(p.theta_c, p.E_iso, p.Gamma0))


def test_swallowed_likelihood_error_is_logged(caplog):
    """A failing evaluation still returns -inf to the sampler but logs the cause once."""
    utils._reported_failures.clear()
    fitter = Fitter(z=0.5, lumi_dist=3e28, jet=lambda p: 1 / 0, medium="ism")
    fitter.add_flux_density(nu=1e15, t=[1e4, 1e5], f_nu=[1e-26, 1e-27], err=[1e-28, 1e-28])
    defs = [ParamDef("n_ism", 1e-3, 10, Scale.log)]
    like = utils.AfterglowLikelihood(
        fitter, defs, lambda c: -0.5 * c, utils._build_transformer(defs)
    )
    with caplog.at_level(logging.WARNING, logger=utils.logger.name):
        assert like.log_likelihood({"log10_n_ism": 0.0}) == -np.inf
        assert like.log_likelihood({"log10_n_ism": 0.5}) == -np.inf
    msgs = [r.getMessage() for r in caplog.records if "ZeroDivisionError" in r.getMessage()]
    assert len(msgs) == 1
