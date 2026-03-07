#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


def post_json(base_url: str, path: str, payload: dict) -> dict:
    req = Request(
        f"{base_url}{path}",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    with urlopen(req, timeout=20) as resp:
        return json.loads(resp.read().decode("utf-8"))


def expect_http_400(base_url: str, path: str, payload: dict) -> None:
    req = Request(
        f"{base_url}{path}",
        data=json.dumps(payload).encode("utf-8"),
        headers={"Content-Type": "application/json"},
        method="POST",
    )
    try:
        with urlopen(req, timeout=20):
            raise AssertionError(f"expected HTTP 400 for {path}")
    except HTTPError as exc:
        if exc.code != 400:
            raise AssertionError(f"expected HTTP 400 for {path}, got {exc.code}") from exc


def main() -> int:
    base = os.getenv("BASE_URL", "http://127.0.0.1:8000").rstrip("/")

    lc = post_json(
        base,
        "/api/lightcurve",
        {
            "frequencies_input": "1e9",
            "include_figure": False,
            "include_exports": False,
        },
    )
    assert "meta" in lc and "compute_seconds" in lc["meta"], "lightcurve meta.compute_seconds missing"
    assert "figure" not in lc and "exports" not in lc, "lightcurve include flags broken"

    sed = post_json(
        base,
        "/api/spectrum",
        {
            "t_snapshots_input": "1e4",
            "include_figure": False,
            "include_exports": False,
        },
    )
    assert "meta" in sed and "compute_seconds" in sed["meta"], "spectrum meta.compute_seconds missing"
    assert "figure" not in sed and "exports" not in sed, "spectrum include flags broken"

    sky = post_json(
        base,
        "/api/skymap",
        {
            "nu_input": "1e9",
            "include_figure": False,
            "include_exports": False,
        },
    )
    assert "meta" in sky and "compute_seconds" in sky["meta"], "skymap meta.compute_seconds missing"
    assert "figure" not in sky and "exports" not in sky, "skymap include flags broken"

    expect_http_400(base, "/api/lightcurve", {"t_min": 10, "t_max": 1, "include_figure": False, "include_exports": False})

    print("backend smoke ok")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (URLError, HTTPError, AssertionError) as exc:
        print(f"backend smoke failed: {exc}", file=sys.stderr)
        raise SystemExit(1)
