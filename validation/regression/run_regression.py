#!/usr/bin/env python3
"""Regression test runner for VegasAfterglow power-law scaling verification."""

import argparse
import json
import sys
from datetime import datetime
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet, Wind

sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from validation.regression.utils import STANDARD_PARAMS, fit_powerlaw

from validation.colors import (_bold, _green, _red, _yellow, _cyan, _dim, _bold_green, _bold_red,
                               _header, _subheader, _bar)

TIME_RANGES = {
    "coasting":        {"ISM": (1e-1, 1),    "wind": (1e-2, 1e-1)},
    "BM":              {"ISM": (5e2, 5e3),    "wind": (5e3, 5e4)},
    "deep_newtonian":  {"ISM": (1e12, 1e13),  "wind": (1e14, 1e15)},
}

VIZ_TIME_GRID = {"ISM": (1e-2, 1e13, 400), "wind": (1e-3, 1e15, 400)}
VIZ_PARAMS = {"ISM": {"Gamma0": 300, "n_ism": 0.1}, "wind": {"Gamma0": 50, "A_star": 1e-1}}

SHOCK_SCALINGS = {
    "coasting":       {"ISM": {"u": F(0), "r": F(1), "B": F(0), "N_p": F(3)},
                       "wind": {"u": F(0), "r": F(1), "B": F(-1), "N_p": F(1)}},
    "BM":             {"ISM": {"u": F(-3,8), "r": F(1,4), "B": F(-3,8), "N_p": F(3,4)},
                       "wind": {"u": F(-1,4), "r": F(1,2), "B": F(-3,4), "N_p": F(1,2)}},
    "deep_newtonian": {"ISM": {"u": F(-3,5), "r": F(2,5), "B": F(-3,5), "N_p": F(6,5)},
                       "wind": {"u": F(-1,3), "r": F(2,3), "B": F(-1), "N_p": F(2,3)}},
}

FREQ_SCALINGS = {
    "coasting":       {"ISM": {"nu_m": F(0), "nu_c": F(-2), "nu_M": F(0)},
                       "wind": {"nu_m": F(-1), "nu_c": F(-1), "nu_M": F(0)}},
    "BM":             {"ISM": {"nu_m": F(-3,2), "nu_c": F(-1,2), "nu_M": F(-3,8)},
                       "wind": {"nu_m": F(-3,2), "nu_c": F(1,2), "nu_M": F(-1,4)}},
    "deep_newtonian": {"ISM": {"nu_m": F(-3,5), "nu_c": F(-1,5), "nu_M": F(0)},
                       "wind": {"nu_m": F(-1), "nu_c": F(1), "nu_M": F(0)}},
}

SLOPE_TOLERANCE = 0.1
SPECTRAL_TOLERANCE = 0.15

# --- Reverse Shock Constants ---

RVS_TIME_RANGES_THIN = {
    "crossing":        {"ISM": (1e-2, 1e-1),     "wind": (5e-1, 5)},
    "BM":              {"ISM": (1e5, 1e6),      "wind": (1e6, 1e7)},
    "deep_newtonian":  {"ISM": (1e12, 1e13),    "wind": (1e14, 1e15)},
}

RVS_TIME_RANGES_THICK = {
    "crossing":        {"ISM": (1e2, 1e3),     "wind": (1e0, 1e1)},
    "BM":              {"ISM": (1e7, 1e8),      "wind": (1e7, 1e8)},
    "deep_newtonian":  {"ISM": (1e12, 1e13),    "wind": (1e14, 1e15)},
}

# Thin shell (small duration τ) — expected power-law indices
RVS_SHOCK_SCALINGS_THIN = {
    "crossing":       {"ISM": {"u": F(3, 2), "r": 1, "B": 0, "N_p": F(3, 2)},
                       "wind": {"u": F(1, 2), "r": 1, "B": -1, "N_p": F(1, 2)}},
    "BM":             {"ISM": {"u": F(-1, 4), "r": F(1, 4), "B": None, "N_p": 0},
                       "wind": {"u": None, "r": F(1, 2), "B": None, "N_p": 0}},
    "deep_newtonian": {"ISM": {"u": F(-2, 5), "r": F(2, 5), "B": None, "N_p": 0},
                       "wind": {"u": None, "r": F(2, 3), "B": None, "N_p": 0}},
}

# Thick shell (large duration τ) — expected power-law indices
RVS_SHOCK_SCALINGS_THICK = {
    "crossing":       {"ISM": {"u": F(1, 4), "r": F(1, 2), "B": F(-1, 4), "N_p": 1},
                       "wind": {"u": 0, "r": 1, "B": -1, "N_p": 1}},
    "BM":             {"ISM": {"u": F(-1, 4), "r": F(1, 4), "B": None, "N_p": 0},
                       "wind": {"u": None, "r": F(1, 2), "B": None, "N_p": 0}},
    "deep_newtonian": {"ISM": {"u": F(-2, 5), "r": F(2, 5), "B": None, "N_p": 0},
                       "wind": {"u": None, "r": F(2, 3), "B": None, "N_p": 0}},
}

RVS_FREQ_SCALINGS_THIN = {
    "crossing":       {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
    "BM":             {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
    "deep_newtonian": {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
}

RVS_FREQ_SCALINGS_THICK = {
    "crossing":       {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
    "BM":             {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
    "deep_newtonian": {"ISM": {"nu_m": None, "nu_c": None, "nu_M": None},
                       "wind": {"nu_m": None, "nu_c": None, "nu_M": None}},
}

RVS_VIZ_PARAMS_THIN = {
    "ISM":  {"Gamma0": 300, "n_ism": 1, "duration": .01},
    "wind": {"Gamma0": 50,  "A_star": 1e-1, "duration": .01},
}
RVS_VIZ_PARAMS_THICK = {
    "ISM":  {"Gamma0": 300, "n_ism": 0.01, "duration": 1e5},
    "wind": {"Gamma0": 50,  "A_star": 1e-2, "duration": 1e4},
}

RVS_VIZ_TIME_GRID = {"ISM": (1e-2, 1e14, 400), "wind": (1e-3, 1e15, 400)}

SPECTRAL_REGIMES = {
    "I":   {"order": ("nu_a", "nu_m", "nu_c"),
            "segments": [("below_nu_a", 2.0), ("nu_a_to_nu_m", F(1,3)), ("nu_m_to_nu_c", "-(p-1)/2"), ("above_nu_c", "-p/2")]},
    "II":  {"order": ("nu_m", "nu_a", "nu_c"),
            "segments": [("below_nu_m", 2.0), ("nu_m_to_nu_a", F(5,2)), ("nu_a_to_nu_c", "-(p-1)/2"), ("above_nu_c", "-p/2")]},
    "III": {"order": ("nu_a", "nu_c", "nu_m"),
            "segments": [("below_nu_a", 2.0), ("nu_a_to_nu_c", F(1,3)), ("nu_c_to_nu_m", F(-1,2)), ("above_nu_m", "-p/2")]},
    "IV":  {"order": ("nu_c", "nu_a", "nu_m"),
            "segments": [("below_nu_c", 2), ("nu_c_to_nu_a", F(2)), ("nu_a_to_nu_m", F(-1,2)), ("above_nu_m", "-p/2")]},
    "V":   {"order": ("nu_c", "nu_m", "nu_a"),
            "segments": [("below_nu_c", 2.0), ("nu_c_to_nu_m", 2.0), ("nu_m_to_nu_a", F(5,2)), ("above_nu_a", "-p/2")]},
}

REGIME_TEST_CONFIGS = {
    "I":   {"medium": "ISM", "t": 5e3,  "n_ism": 0.1, "eps_B": 1e-2},
    "II":  {"medium": "ISM", "t": 5e4,  "n_ism": 1e6, "eps_B": 1e-5, "eps_e": 5e-2},
    "III": {"medium": "ISM", "t": 1e5,  "n_ism": 30,  "eps_B": 1e-1, "eps_e": 1e-1, "xi_e": 5e-3},
    "IV":  {"medium": "ISM", "t": 1e5,  "eps_B": 0.3, "eps_e": 0.3,  "n_ism": 1e6,  "xi_e": 5e-3},
    "V":   {"medium": "ISM", "t": 1e4,  "n_ism": 1e8, "eps_e": 0.01, "eps_B": 3e-1, "xi_e": 1},
}

# --- Helpers ---

_MEDIUM_KEY = {"ISM": "ISM", "wind": "Wind"}

def _fwd(details, *keys):
    """Extract forward-shock arrays sorted by t_obs. Returns [t, key1, key2, ...]."""
    t = np.array(details.fwd.t_obs)[0, 0, :]
    idx = np.argsort(t)
    return [t[idx]] + [np.array(getattr(details.fwd, k))[0, 0, :][idx] for k in keys]

def _rvs(details, *keys):
    """Extract reverse-shock arrays sorted by t_obs. Returns [t, key1, key2, ...]."""
    t = np.array(details.rvs.t_obs)[0, 0, :]
    idx = np.argsort(t)
    return [t[idx]] + [np.array(getattr(details.rvs, k))[0, 0, :][idx] for k in keys]

def _eval_beta(expr, p):
    if isinstance(expr, str):
        return eval(expr.replace("p", str(p)))
    return float(expr)

def _beta_expr(b):
    if isinstance(b, F):
        return f"{b.numerator}/{b.denominator}" if b.denominator != 1 else str(b.numerator)
    return str(b) if isinstance(b, str) else (str(int(b)) if b == int(b) else str(b))

def _check(measured, expected, tol):
    return bool(not np.isnan(measured) and abs(measured - expected) < tol)

def _segment_range(seg_name, nu_a, nu_m, nu_c):
    freq_map = {"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}
    if seg_name.startswith("below_"):
        nu = freq_map.get(seg_name[6:])
        return (nu / 1000, nu / 30) if nu else (None, None)
    if seg_name.startswith("above_"):
        nu = freq_map.get(seg_name[6:])
        return (nu * 10, nu * 100) if nu else (None, None)
    parts = seg_name.split("_to_")
    if len(parts) == 2 and (lo := freq_map.get(parts[0])) and (hi := freq_map.get(parts[1])):
        return (lo * 5, hi / 20)
    return None, None

def _detect_regime(nu_a, nu_m, nu_c):
    order = tuple(f[0] for f in sorted({"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}.items(), key=lambda x: x[1]))
    return next((name for name, info in SPECTRAL_REGIMES.items() if info["order"] == order), None)

def _measure_segment_slopes(model, t_test, regime_name, nu_a, nu_m, nu_c, p):
    """Measure spectral slopes for all segments of a regime. Returns list of result dicts."""
    results = []
    for seg_name, expected_beta in SPECTRAL_REGIMES[regime_name]["segments"]:
        nu_low, nu_high = _segment_range(seg_name, nu_a, nu_m, nu_c)
        if nu_low is None or nu_high is None or nu_high <= nu_low * 1.5:
            continue
        expected_val = _eval_beta(expected_beta, p)
        nu_grid = np.logspace(np.log10(nu_low * 1.2), np.log10(nu_high * 0.8), 20)
        flux = model.flux_density(np.full_like(nu_grid, t_test), nu_grid)
        measured = fit_powerlaw(nu_grid, flux.total)
        results.append({
            "segment": seg_name, "expected": expected_val, "expected_expr": _beta_expr(expected_beta),
            "measured": float(measured) if not np.isnan(measured) else None,
            "passed": _check(measured, expected_val, SPECTRAL_TOLERANCE),
            "tolerance": SPECTRAL_TOLERANCE,
            "nu_low": float(nu_low), "nu_high": float(nu_high),
        })
    return results


def create_model(medium_type, **overrides):
    defaults = {k: STANDARD_PARAMS[k] for k in ("Gamma0", "A_star", "eps_e", "eps_B", "xi_e", "n_ism")}
    p = {**defaults, **{k: v for k, v in overrides.items() if v is not None}}
    jet = TophatJet(theta_c=STANDARD_PARAMS["theta_c"], E_iso=STANDARD_PARAMS["E_iso"], Gamma0=p["Gamma0"])
    medium = ISM(n_ism=p["n_ism"]) if medium_type == "ISM" else Wind(A_star=p["A_star"])
    observer = Observer(lumi_dist=STANDARD_PARAMS["lumi_dist"], z=STANDARD_PARAMS["z"], theta_obs=0.0)
    radiation = Radiation(eps_e=p["eps_e"], eps_B=p["eps_B"], p=STANDARD_PARAMS["p"], xi_e=p["xi_e"])
    return Model(jet, medium, observer, radiation, resolutions=(0.3, 2, 15))


def create_model_with_rvs(medium_type, duration=1, **overrides):
    """Create a model with both forward and reverse shock radiation."""
    defaults = {k: STANDARD_PARAMS[k] for k in ("Gamma0", "A_star", "eps_e", "eps_B", "xi_e", "n_ism")}
    p = {**defaults, **{k: v for k, v in overrides.items() if v is not None}}
    jet = TophatJet(theta_c=STANDARD_PARAMS["theta_c"], E_iso=STANDARD_PARAMS["E_iso"],
                    Gamma0=p["Gamma0"], duration=duration)
    medium = ISM(n_ism=p["n_ism"]) if medium_type == "ISM" else Wind(A_star=p["A_star"])
    observer = Observer(lumi_dist=STANDARD_PARAMS["lumi_dist"], z=STANDARD_PARAMS["z"], theta_obs=0.0)
    fwd_rad = Radiation(eps_e=p["eps_e"], eps_B=p["eps_B"], p=STANDARD_PARAMS["p"], xi_e=p["xi_e"])
    rvs_rad = Radiation(eps_e=p.get("eps_e_r", p["eps_e"]),
                        eps_B=p.get("eps_B_r", p["eps_B"]),
                        p=p.get("p_r", STANDARD_PARAMS["p"]),
                        xi_e=p.get("xi_e_r", p["xi_e"]))
    return Model(jet, medium, observer, fwd_rad, rvs_rad, resolutions=(0.3, 2, 15))


class RegressionRunner:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(__file__).parent / output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {"timestamp": datetime.now().isoformat(), "tests": [], "categories": {},
                        "summary": {"by_model": {"ISM": {"pass": 0, "fail": 0}, "Wind": {"pass": 0, "fail": 0}}, "by_category": {}}}
        self.p = STANDARD_PARAMS["p"]
        self._models = {}

    def _get_viz_model(self, medium_type):
        key = f"viz_{medium_type}"
        if key not in self._models:
            self._models[key] = create_model(medium_type, **VIZ_PARAMS.get(medium_type, {}))
        return self._models[key]

    def _record(self, result, category, medium):
        self.results["tests"].append(result)
        self.results["categories"].setdefault(category, []).append(result)
        stat = "pass" if result.get("passed", False) else "fail"
        self.results["summary"]["by_model"][_MEDIUM_KEY.get(medium, medium)][stat] += 1
        self.results["summary"]["by_category"].setdefault(category, {"pass": 0, "fail": 0})[stat] += 1

    def _status_line(self, measured, expected, passed):
        tag = _bold_green("PASS") if passed else _bold_red("FAIL")
        if np.isnan(measured):
            return f"      {tag}: measurement failed"
        return f"      {tag}: measured={measured:.3f}, expected={expected:.3f}"

    # --- Core test methods ---

    def run_shock_dynamics(self):
        print(_subheader("1. Forward Shock Dynamics (All 3 Phases)"))
        self.results.setdefault("shock_grid", {})
        for medium in ["ISM", "wind"]:
            model, mk = self._get_viz_model(medium), _MEDIUM_KEY[medium]
            print(f"\n{_bold(f'--- {mk} ---')}")
            self.results["shock_grid"].setdefault(mk, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {_cyan(phase)}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], SHOCK_SCALINGS[phase][medium], {}
                for qty in ["u", "r", "B", "N_p"]:
                    if (exp_val := expected.get(qty)) is None:
                        continue
                    result = self._run_shock_quantity(model, t_range, qty, exp_val, f"{mk} {phase}: {qty}", phase)
                    self._record(result, "shock_dynamics", mk)
                    phase_results[qty] = result
                self.results["shock_grid"][mk][phase] = phase_results

    def _run_shock_quantity(self, model, t_range, qty, expected, name, phase="BM",
                            extractor=_fwd, expand_range=False, min_valid=5, type_prefix="shock"):
        """Run a single shock dynamics test for forward or reverse shock."""
        print(f"    Testing: {_bold(qty)}")
        attr_map = {"u": None, "r": "r", "B": "B_comv", "N_p": "N_p"}
        try:
            t_lo = t_range[0] / 10 if expand_range else t_range[0]
            t_hi = t_range[1] * 10 if expand_range else t_range[1]
            details = model.details(t_lo, t_hi)
            arrays = extractor(details, "Gamma_th") if qty == "u" else extractor(details, attr_map[qty])
            t = arrays[0]
            mask = (t >= t_range[0]) & (t <= t_range[1])
            t = t[mask]
            if qty == "u":
                Gamma = arrays[1][mask]
                y = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                valid = (y > 0) & np.isfinite(y)
                if np.sum(valid) < min_valid:
                    raise ValueError(f"Not enough valid points ({np.sum(valid)})")
                measured = fit_powerlaw(t[valid], y[valid])
            else:
                measured = fit_powerlaw(t, arrays[1][mask])
        except Exception as e:
            print(f"      {_bold_red('ERROR')}: {e}")
            measured = np.nan
        tol = SLOPE_TOLERANCE
        passed = _check(measured, expected, tol)
        result = {"name": name, "quantity": qty, "phase": phase, "type": f"{type_prefix}_{qty}", "expected": float(expected),
                  "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tol, "passed": passed}
        print(self._status_line(measured, float(expected), passed))
        return result

    def run_characteristic_frequencies(self):
        print(_subheader("2. Forward Shock Frequencies (All 3 Phases)"))
        self.results.setdefault("freq_grid", {})
        for medium in ["ISM", "wind"]:
            model, mk = self._get_viz_model(medium), _MEDIUM_KEY[medium]
            print(f"\n{_bold(f'--- {mk} ---')}")
            self.results["freq_grid"].setdefault(mk, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {_cyan(phase)}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], FREQ_SCALINGS[phase][medium], {}
                for freq_name in ["nu_m", "nu_c", "nu_M"]:
                    if (exp_val := expected.get(freq_name)) is None:
                        continue
                    result = self._run_nu_scaling(model, t_range, freq_name, exp_val, f"{mk} {phase}: {freq_name}", phase)
                    self._record(result, "frequencies", mk)
                    phase_results[freq_name] = result
                self.results["freq_grid"][mk][phase] = phase_results

    def _run_nu_scaling(self, model, t_range, freq_name, expected, name, phase="BM",
                        extractor=_fwd, expand_range=False, min_valid=5, type_name="frequency"):
        """Run a single frequency scaling test for forward or reverse shock."""
        print(f"    Testing: {_bold(freq_name)}")
        try:
            t_lo = t_range[0] / 10 if expand_range else t_range[0]
            t_hi = t_range[1] * 10 if expand_range else t_range[1]
            details = model.details(t_lo, t_hi)
            t, Doppler, nu = extractor(details, "Doppler", freq_name)
            nu = nu * Doppler
            mask = (t >= t_range[0]) & (t <= t_range[1])
            t, nu = t[mask], nu[mask]
            valid = (nu > 0) & np.isfinite(nu)
            if np.sum(valid) < min_valid:
                raise ValueError(f"Not enough valid points ({np.sum(valid)})")
            measured = fit_powerlaw(t[valid], nu[valid])
        except Exception as e:
            print(f"      {_bold_red('ERROR')}: {e}")
            measured = np.nan
        tol = SLOPE_TOLERANCE
        passed = _check(measured, expected, tol)
        result = {"name": name, "quantity": freq_name, "phase": phase, "type": type_name, "expected": float(expected),
                  "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tol, "passed": passed}
        print(self._status_line(measured, float(expected), passed))
        return result

    def run_spectrum_shapes(self):
        print(_subheader("3. Forward Shock Spectrum Shapes (5 Regimes)"))
        self.results.setdefault("spectrum_grid", {})
        for regime_name, config in REGIME_TEST_CONFIGS.items():
            medium, t_test = config["medium"], config["t"]
            mk = _MEDIUM_KEY.get(medium, medium)
            print(f"\n  Regime {_bold(regime_name)}")
            model = create_model(medium, **{k: config[k] for k in ("eps_e", "eps_B", "xi_e", "n_ism", "A_star") if k in config})
            try:
                details = model.details(t_test * 0.9, t_test * 1.1)
                idx = np.argmin(np.abs(np.array(details.fwd.t_obs)[0, 0, :] - t_test))
                D = np.array(details.fwd.Doppler)[0, 0, idx]
                nu_a = np.array(details.fwd.nu_a)[0, 0, idx] * D
                nu_m = np.array(details.fwd.nu_m)[0, 0, idx] * D
                nu_c = np.array(details.fwd.nu_c)[0, 0, idx] * D
                print(f"    ν_a = {nu_a:.2e}, ν_m = {nu_m:.2e}, ν_c = {nu_c:.2e}")

                actual = _detect_regime(nu_a, nu_m, nu_c) or regime_name
                if actual != regime_name:
                    print(f"    {_yellow('NOTE')}: Expected regime {regime_name}, found regime {actual}")

                seg_results = _measure_segment_slopes(model, t_test, regime_name, nu_a, nu_m, nu_c, self.p)
                for sr in seg_results:
                    result = {"name": f"Regime {regime_name}: {sr['segment']}", "regime": regime_name, "segment": sr["segment"],
                              "type": "spectrum_shape", **{k: sr[k] for k in ("expected", "expected_expr", "measured", "tolerance", "passed")}}
                    self._record(result, "spectrum_shapes", mk)
                    if sr['measured'] is not None:
                        tag = _bold_green("PASS") if sr['passed'] else _bold_red("FAIL")
                        print(f"    {sr['segment']}: β = {sr['measured']:.2f} (expected {sr['expected']:.2f}) [{tag}]")
                    else:
                        print(f"    {sr['segment']}: {_bold_red('measurement failed')}")

                self.results["spectrum_grid"].setdefault(mk, {})[regime_name] = {
                    "t": t_test, "nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c, "actual_regime": actual,
                    "segments": [{"name": f"Regime {regime_name}: {s['segment']}", **s} for s in seg_results]}
            except Exception as e:
                print(f"    {_bold_red('ERROR')}: {e}")

    # --- Reverse Shock Tests ---

    def _get_rvs_viz_model(self, medium_type, regime):
        """Get or create reverse shock model for visualization. regime='thin' or 'thick'."""
        key = f"rvs_viz_{regime}_{medium_type}"
        if key not in self._models:
            params = RVS_VIZ_PARAMS_THIN if regime == "thin" else RVS_VIZ_PARAMS_THICK
            self._models[key] = create_model_with_rvs(medium_type, **params.get(medium_type, {}))
        return self._models[key]

    def _run_rvs_shock_quantity(self, model, t_range, qty, expected, name, phase="BM"):
        """Run a single reverse-shock dynamics test."""
        return self._run_shock_quantity(
            model, t_range, qty, expected, name, phase,
            extractor=_rvs, expand_range=True, min_valid=3, type_prefix="rvs_shock")

    def _run_rvs_nu_scaling(self, model, t_range, freq_name, expected, name, phase="BM"):
        """Run a single reverse-shock frequency scaling test."""
        return self._run_nu_scaling(
            model, t_range, freq_name, expected, name, phase,
            extractor=_rvs, expand_range=True, min_valid=3, type_name="rvs_frequency")

    def run_rvs_shock_dynamics(self):
        """Test reverse shock dynamics across thin/thick shell regimes and all phases."""
        for i_regime, (regime, scalings) in enumerate([("thin", RVS_SHOCK_SCALINGS_THIN),
                                                        ("thick", RVS_SHOCK_SCALINGS_THICK)]):
            section_num = 4 + i_regime
            print(_subheader(f"{section_num}. Reverse Shock Dynamics — {regime.title()} Shell (All Phases)"))
            grid_key = f"rvs_shock_grid_{regime}"
            self.results.setdefault(grid_key, {})
            for medium in ["ISM", "wind"]:
                model, mk = self._get_rvs_viz_model(medium, regime), _MEDIUM_KEY[medium]
                print(f"\n{_bold(f'--- {mk} ---')}")
                self.results[grid_key].setdefault(mk, {})
                for phase in ["crossing", "BM", "deep_newtonian"]:
                    expected = scalings[phase][medium]
                    # Skip phases with no expected values (TBD)
                    if all(v is None for v in expected.values()):
                        print(f"\n  Phase: {_cyan(phase)} {_dim('(skipped — scalings TBD)')}")
                        continue
                    print(f"\n  Phase: {_cyan(phase)}")
                    rvs_time_ranges = RVS_TIME_RANGES_THIN if regime == "thin" else RVS_TIME_RANGES_THICK
                    t_range = rvs_time_ranges[phase][medium]
                    phase_results = {}
                    for qty in ["u", "r", "B", "N_p"]:
                        exp_val = expected.get(qty)
                        if exp_val is None:
                            continue
                        result = self._run_rvs_shock_quantity(
                            model, t_range, qty, exp_val,
                            f"Reverse Shock {regime} {mk} {phase}: {qty}", phase)
                        self._record(result, f"rvs_shock_dynamics_{regime}", mk)
                        phase_results[qty] = result
                    self.results[grid_key][mk][phase] = phase_results

    def run_rvs_characteristic_frequencies(self):
        """Test reverse shock characteristic frequency scalings across thin/thick shell regimes."""
        for i_regime, (regime, scalings) in enumerate([("thin", RVS_FREQ_SCALINGS_THIN),
                                                        ("thick", RVS_FREQ_SCALINGS_THICK)]):
            section_num = 6 + i_regime
            print(_subheader(f"{section_num}. Reverse Shock Frequencies — {regime.title()} Shell (All Phases)"))
            grid_key = f"rvs_freq_grid_{regime}"
            self.results.setdefault(grid_key, {})
            for medium in ["ISM", "wind"]:
                model, mk = self._get_rvs_viz_model(medium, regime), _MEDIUM_KEY[medium]
                print(f"\n{_bold(f'--- {mk} ---')}")
                self.results[grid_key].setdefault(mk, {})
                for phase in ["crossing", "BM", "deep_newtonian"]:
                    expected = scalings[phase][medium]
                    if all(v is None for v in expected.values()):
                        print(f"\n  Phase: {_cyan(phase)} {_dim('(skipped — scalings TBD)')}")
                        continue
                    print(f"\n  Phase: {_cyan(phase)}")
                    rvs_time_ranges = RVS_TIME_RANGES_THIN if regime == "thin" else RVS_TIME_RANGES_THICK
                    t_range = rvs_time_ranges[phase][medium]
                    phase_results = {}
                    for freq_name in ["nu_m", "nu_c", "nu_M"]:
                        exp_val = expected.get(freq_name)
                        if exp_val is None:
                            continue
                        result = self._run_rvs_nu_scaling(
                            model, t_range, freq_name, exp_val,
                            f"Reverse Shock {regime} {mk} {phase}: {freq_name}", phase)
                        self._record(result, f"rvs_frequencies_{regime}", mk)
                        phase_results[freq_name] = result
                    self.results[grid_key][mk][phase] = phase_results

    # --- Visualization data collection ---

    @staticmethod
    def _make_phase_mask(t, u, t_range, phase):
        """Get mask for a phase time range, returns None if insufficient points."""
        mask = (t >= t_range[0]) & (t <= t_range[1])
        if phase == "deep_newtonian":
            mask &= (u < 0.1)
        return mask if np.sum(mask) >= 5 else None

    @staticmethod
    def _collect_shock_freq_data(t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M,
                                  time_ranges, medium, shock_scalings, freq_scalings):
        """Collect shock dynamics and frequency viz data for a set of phases.

        Returns (shock_data_dict, freq_data_dict).
        """
        shock_phases = {}
        for phase in time_ranges:
            tr = time_ranges[phase][medium]
            mask = RegressionRunner._make_phase_mask(t, u, tr, phase)
            if mask is None:
                continue
            tp = t[mask]
            fits = {}
            for name, arr, key in [("r", r, "r"), ("u", u, "u"), ("B", B, "B"), ("N_p", N_p, "N_p")]:
                vals = arr[mask]
                valid = (vals > 0) & np.isfinite(vals)
                if np.sum(valid) > 5:
                    exp = shock_scalings.get(phase, {}).get(medium, {}).get(key)
                    fits[name] = {"measured": float(fit_powerlaw(tp[valid], vals[valid])),
                                  "expected": float(exp) if exp is not None else None}
            shock_phases[phase] = {"t_range": [float(tp.min()), float(tp.max())], "fits": fits}

        freq_phases = {}
        for phase in time_ranges:
            tr = time_ranges[phase][medium]
            mask = RegressionRunner._make_phase_mask(t, u, tr, phase)
            if mask is None:
                continue
            tp = t[mask]
            fits = {}
            for fname, farr in [("nu_m", nu_m), ("nu_c", nu_c), ("nu_a", nu_a), ("nu_M", nu_M)]:
                valid = (farr[mask] > 0) & np.isfinite(farr[mask])
                if np.sum(valid) > 5:
                    exp = freq_scalings.get(phase, {}).get(medium, {}).get(fname)
                    fits[fname] = {"measured": float(fit_powerlaw(tp[valid], farr[mask][valid])),
                                   "expected": float(exp) if exp is not None else None}
            freq_phases[phase] = {"t_range": [float(tp.min()), float(tp.max())], "fits": fits}

        shock_data = {"t": t.tolist(), "u": u.tolist(), "Gamma": Gamma.tolist(),
                      "r": r.tolist(), "B": B.tolist(), "N_p": N_p.tolist(), "phases": shock_phases}
        freq_data = {"t": t.tolist(), "Doppler": Doppler.tolist(), "nu_m": nu_m.tolist(),
                     "nu_c": nu_c.tolist(), "nu_a": nu_a.tolist(), "nu_M": nu_M.tolist(), "phases": freq_phases}
        return shock_data, freq_data

    @staticmethod
    def _extract_arrays(details, extractor):
        """Extract and process shock/freq arrays from model details."""
        t, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M = extractor(
            details, "Gamma_th", "r", "B_comv", "N_p", "Doppler", "nu_m", "nu_c", "nu_a", "nu_M")
        u = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
        nu_m, nu_c, nu_a, nu_M = nu_m * Doppler, nu_c * Doppler, nu_a * Doppler, nu_M * Doppler
        return t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M

    def _collect_viz(self):
        """Collect all visualization data in one pass."""
        print(_subheader("Collecting Visualization Data"))
        p = self.p

        # Forward shock
        for medium in ["ISM", "wind"]:
            mk = _MEDIUM_KEY[medium]
            viz_model = self._get_viz_model(medium)
            t_min, t_max, _ = VIZ_TIME_GRID[medium]
            print(f"{_bold(f'--- {mk} ---')}")
            try:
                details = viz_model.details(t_min, t_max)
                t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M = self._extract_arrays(details, _fwd)
                shock_data, freq_data = self._collect_shock_freq_data(
                    t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M,
                    TIME_RANGES, medium, SHOCK_SCALINGS, FREQ_SCALINGS)
                self.results.setdefault("viz_shock_dynamics", {})[mk] = shock_data
                self.results.setdefault("viz_frequencies", {})[mk] = freq_data
                print(f"  Shock/freq: {_dim(f'{len(t)} points')}")
            except Exception as e:
                print(f"  Shock/freq {_bold_red('ERROR')}: {e}")

        # Spectrum shapes (regime-specific)
        viz_shapes = {}
        for regime_name, config in REGIME_TEST_CONFIGS.items():
            medium, t_test = config["medium"], config["t"]
            model = create_model(medium, **{k: config[k] for k in ("eps_e", "eps_B", "xi_e", "n_ism", "A_star") if k in config})
            try:
                details = model.details(t_test * 0.9, t_test * 1.1)
                idx = np.argmin(np.abs(np.array(details.fwd.t_obs)[0, 0, :] - t_test))
                D = np.array(details.fwd.Doppler)[0, 0, idx]
                nu_a, nu_m, nu_c = [np.array(getattr(details.fwd, f))[0, 0, idx] * D for f in ("nu_a", "nu_m", "nu_c")]
                nu_g = np.logspace(6, 20, 200)
                flux = model.flux_density(np.full_like(nu_g, t_test), nu_g)
                actual = _detect_regime(nu_a, nu_m, nu_c) or regime_name
                viz_shapes[f"regime_{regime_name}"] = {
                    "medium": _MEDIUM_KEY[medium], "t": t_test, "nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c,
                    "actual_regime": actual, "nu": nu_g.tolist(), "flux": flux.total.tolist(),
                    "expected_segments": _measure_segment_slopes(model, t_test, regime_name, nu_a, nu_m, nu_c, p)}
                print(f"  Regime {regime_name}: {_dim('collected')}")
            except Exception as e:
                print(f"  Regime {regime_name}: {_bold_red('ERROR')} - {e}")
        self.results["viz_spectrum_shapes"] = viz_shapes

        # Reverse shock
        rvs_config = {
            "thin":  (RVS_SHOCK_SCALINGS_THIN, RVS_FREQ_SCALINGS_THIN, RVS_TIME_RANGES_THIN),
            "thick": (RVS_SHOCK_SCALINGS_THICK, RVS_FREQ_SCALINGS_THICK, RVS_TIME_RANGES_THICK),
        }
        for regime, (shock_scalings, freq_scalings, time_ranges) in rvs_config.items():
            for medium in ["ISM", "wind"]:
                mk = _MEDIUM_KEY[medium]
                rvs_model = self._get_rvs_viz_model(medium, regime)
                t_min, t_max, _ = RVS_VIZ_TIME_GRID[medium]
                print(f"{_bold(f'--- Reverse Shock {regime.title()} {mk} ---')}")
                try:
                    details = rvs_model.details(t_min, t_max)
                    t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M = self._extract_arrays(details, _rvs)
                    shock_data, freq_data = self._collect_shock_freq_data(
                        t, u, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a, nu_M,
                        time_ranges, medium, shock_scalings, freq_scalings)
                    self.results.setdefault(f"viz_rvs_shock_dynamics_{regime}", {})[mk] = shock_data
                    self.results.setdefault(f"viz_rvs_frequencies_{regime}", {})[mk] = freq_data
                    print(f"  Reverse Shock {regime} shock/freq: {_dim(f'{len(t)} points')}")
                except Exception as e:
                    print(f"  Reverse Shock {regime} shock/freq {_bold_red('ERROR')}: {e}")

    # --- Run all + summary ---

    def run_all(self, viz=True):
        print(_header("COMPREHENSIVE REGRESSION TEST SUITE"))
        # Forward shock
        self.run_shock_dynamics()
        self.run_characteristic_frequencies()
        self.run_spectrum_shapes()
        # Reverse shock
        self.run_rvs_shock_dynamics()
        self.run_rvs_characteristic_frequencies()
        if viz:
            self._collect_viz()
        self._print_summary()
        return self.results

    def _print_summary(self):
        print(_header("SUMMARY"))
        print(f"\n{_bold('By Model:')}")
        for model, stats in self.results["summary"]["by_model"].items():
            n_pass, n_fail = stats["pass"], stats["fail"]
            total = n_pass + n_fail
            if total > 0:
                pct = n_pass / total * 100
                color = _bold_green if pct == 100 else (_bold_red if pct < 50 else _yellow)
                fail_str = _red(f"{n_fail:2d} fail") if n_fail else _dim(f"{n_fail:2d} fail")
                print(f"  {model:6s}: {_green(f'{n_pass:2d} pass')}, {fail_str} ({color(f'{pct:.0f}%')})")
        print(f"\n{_bold('By Category:')}")
        for cat, stats in self.results["summary"]["by_category"].items():
            n_pass, n_fail = stats["pass"], stats["fail"]
            total = n_pass + n_fail
            if total > 0:
                pct = n_pass / total * 100
                color = _bold_green if pct == 100 else (_bold_red if pct < 50 else _yellow)
                fail_str = _red(f"{n_fail:2d} fail") if n_fail else _dim(f"{n_fail:2d} fail")
                print(f"  {cat:20s}: {_green(f'{n_pass:2d} pass')}, {fail_str} ({color(f'{pct:.0f}%')})")
        n_tests = len(self.results["tests"])
        n_pass = sum(1 for t in self.results["tests"] if t.get("passed", False))
        pct = n_pass / n_tests * 100
        overall_color = _bold_green if pct == 100 else _bold_red
        print(f"\n{_bold('Overall')}: {overall_color(f'{n_pass}/{n_tests} passed ({pct:.0f}%)')}\n{_bar()}")

    def save_results(self, filename="regression_results.json"):
        path = self.output_dir / filename
        path.write_text(json.dumps(self.results, indent=2, default=str))
        print(f"\nResults saved to {path}")
        return path


def main():
    parser = argparse.ArgumentParser(description="Run VegasAfterglow regression tests")
    parser.add_argument("--no-viz", action="store_true", help="Skip visualization data collection")
    parser.add_argument("--output", default="results", help="Output directory for results")
    args = parser.parse_args()
    runner = RegressionRunner(output_dir=args.output)
    runner.run_all(viz=not args.no_viz)
    runner.save_results()


if __name__ == "__main__":
    main()
