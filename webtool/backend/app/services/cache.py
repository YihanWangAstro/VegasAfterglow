from __future__ import annotations

from functools import lru_cache
from typing import Any

import orjson

from ..compute import compute_model, compute_sed, compute_skymap


def canonical_json(value: Any) -> str:
    return orjson.dumps(value, option=orjson.OPT_SORT_KEYS | orjson.OPT_NON_STR_KEYS).decode()


@lru_cache(maxsize=64)
def cached_compute_model(params_json: str):
    return compute_model(orjson.loads(params_json))


@lru_cache(maxsize=64)
def cached_compute_sed(params_json: str):
    return compute_sed(orjson.loads(params_json))


@lru_cache(maxsize=32)
def cached_compute_skymap(params_json: str):
    return compute_skymap(orjson.loads(params_json))
