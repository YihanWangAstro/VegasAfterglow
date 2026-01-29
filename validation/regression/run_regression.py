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

TIME_RANGES = {
    "coasting":        {"ISM": (1e-1, 1),    "wind": (1e-2, 1e-1)},
    "BM":              {"ISM": (5e2, 5e3),    "wind": (5e3, 5e4)},
    "deep_newtonian":  {"ISM": (1e12, 1e13),  "wind": (1e14, 1e15)},
}

JET_BREAK_TIME = 1.8e5
VIZ_TIME_GRID = {"ISM": (1e-2, 1e13, 400), "wind": (1e-3, 1e15, 400)}
VIZ_FREQ_GRID = (1e8, 1e20, 100)
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
    "coasting":       {"ISM": {"nu_m": F(0), "nu_c": F(-2)},    "wind": {"nu_m": F(-1), "nu_c": F(-1)}},
    "BM":             {"ISM": {"nu_m": F(-3,2), "nu_c": F(-1,2)}, "wind": {"nu_m": F(-3,2), "nu_c": F(1,2)}},
    "deep_newtonian": {"ISM": {"nu_m": F(-3,5), "nu_c": F(-1,5)}, "wind": {"nu_m": F(-1), "nu_c": F(1)}},
}

SLOPE_TOLERANCE = 0.1
SPECTRAL_TOLERANCE = 0.15

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
        return (lo * 10, hi / 10)
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


class RegressionRunner:
    def __init__(self, output_dir="results"):
        self.output_dir = Path(__file__).parent / output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results = {"timestamp": datetime.now().isoformat(), "tests": [], "scaling_results": [], "categories": {},
                        "summary": {"by_model": {"ISM": {"pass": 0, "fail": 0}, "Wind": {"pass": 0, "fail": 0}}, "by_category": {}}}
        self.p = STANDARD_PARAMS["p"]
        self._models = {}

    def _get_model(self, medium_type):
        return self._models.setdefault(medium_type, create_model(medium_type))

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
        if np.isnan(measured):
            return f"      {'PASS' if passed else 'FAIL'}: measurement failed"
        return f"      {'PASS' if passed else 'FAIL'}: measured={measured:.3f}, expected={expected:.3f}"

    # --- Core test methods ---

    def run_shock_dynamics(self):
        print("\n=== 1. Shock Dynamics (All 3 Phases) ===\n")
        self.results.setdefault("shock_grid", {})
        for medium in ["ISM", "wind"]:
            model, mk = self._get_viz_model(medium), _MEDIUM_KEY[medium]
            print(f"\n--- {mk} ---")
            self.results["shock_grid"].setdefault(mk, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {phase}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], SHOCK_SCALINGS[phase][medium], {}
                for qty in ["u", "r", "B", "N_p"]:
                    if (exp_val := expected.get(qty)) is None:
                        continue
                    result = self._run_shock_quantity(model, t_range, qty, exp_val, f"{mk} {phase}: {qty}", phase)
                    self._record(result, "shock_dynamics", mk)
                    phase_results[qty] = result
                self.results["shock_grid"][mk][phase] = phase_results

    def _run_shock_quantity(self, model, t_range, qty, expected, name, phase="BM"):
        print(f"    Testing: {qty}")
        attr_map = {"u": None, "r": "r", "B": "B_comv", "N_p": "N_p"}
        try:
            details = model.details(t_range[0], t_range[1])
            arrays = _fwd(details, "Gamma") if qty == "u" else _fwd(details, attr_map[qty])
            t = arrays[0]
            mask = (t >= t_range[0]) & (t <= t_range[1])
            t = t[mask]
            if qty == "u":
                Gamma = arrays[1][mask]
                y = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                valid = (y > 0) & np.isfinite(y)
                if np.sum(valid) < 5:
                    raise ValueError(f"Not enough valid points ({np.sum(valid)})")
                measured = fit_powerlaw(t[valid], y[valid])
            else:
                measured = fit_powerlaw(t, arrays[1][mask])
        except Exception as e:
            print(f"      ERROR: {e}")
            measured = np.nan
        tol = SLOPE_TOLERANCE
        passed = _check(measured, expected, tol)
        result = {"name": name, "quantity": qty, "phase": phase, "type": f"shock_{qty}", "expected": float(expected),
                  "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tol, "passed": passed}
        print(self._status_line(measured, float(expected), passed))
        return result

    def run_characteristic_frequencies(self):
        print("\n=== 2. Characteristic Frequencies (All 3 Phases) ===\n")
        self.results.setdefault("freq_grid", {})
        for medium in ["ISM", "wind"]:
            model, mk = self._get_viz_model(medium), _MEDIUM_KEY[medium]
            print(f"\n--- {mk} ---")
            self.results["freq_grid"].setdefault(mk, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {phase}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], FREQ_SCALINGS[phase][medium], {}
                for freq_name in ["nu_m", "nu_c"]:
                    if (exp_val := expected.get(freq_name)) is None:
                        continue
                    result = self._run_nu_scaling(model, t_range, freq_name, exp_val, f"{mk} {phase}: {freq_name}", phase)
                    self._record(result, "frequencies", mk)
                    phase_results[freq_name] = result
                self.results["freq_grid"][mk][phase] = phase_results

    def _run_nu_scaling(self, model, t_range, freq_name, expected, name, phase="BM"):
        print(f"    Testing: {freq_name}")
        try:
            details = model.details(t_range[0], t_range[1])
            t, Doppler, nu = _fwd(details, "Doppler", freq_name)
            nu = nu * Doppler
            mask = (t >= t_range[0]) & (t <= t_range[1])
            t, nu = t[mask], nu[mask]
            valid = (nu > 0) & np.isfinite(nu)
            if np.sum(valid) < 5:
                raise ValueError(f"Not enough valid points ({np.sum(valid)})")
            measured = fit_powerlaw(t[valid], nu[valid])
        except Exception as e:
            print(f"      ERROR: {e}")
            measured = np.nan
        tol = SLOPE_TOLERANCE
        passed = _check(measured, expected, tol)
        result = {"name": name, "quantity": freq_name, "phase": phase, "type": "frequency", "expected": float(expected),
                  "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tol, "passed": passed}
        print(self._status_line(measured, float(expected), passed))
        return result

    def run_flux_regimes(self):
        print("\n=== 3. Flux Scalings (Regime I) ===\n")
        p = self.p

        def _run(model, medium, nu_or_t, fixed, is_temporal, expected, name, tol):
            print(f"  Testing: {name}")
            grid = np.logspace(*[np.log10(x) for x in nu_or_t], 30)
            try:
                if is_temporal:
                    flux = model.flux_density(grid, np.full_like(grid, fixed))
                else:
                    flux = model.flux_density(np.full_like(grid, fixed), grid)
                measured = fit_powerlaw(grid, flux.total)
            except Exception as e:
                print(f"    ERROR: {e}")
                measured = np.nan
            passed = _check(measured, expected, tol * abs(expected) if is_temporal else tol)
            result = {"name": name, "type": "temporal" if is_temporal else "spectral", "expected": expected,
                      "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tol, "passed": passed}
            self.results["tests"].append(result)
            self.results["scaling_results"].append(result)
            self._record(result, "flux_regimes", medium)
            print(f"    {'PASS' if passed else 'FAIL'}: measured={measured:.3f}, expected={expected:.3f}")

        ism = self._get_model("ISM")
        print("--- ISM ---")
        _run(ism, "ISM", (1e4, 5e4),  5e13, True,  -(3*p-3)/4,  "ISM: α (ν_m < ν < ν_c)", 0.2)
        _run(ism, "ISM", (1e3, 5e4),  1e17, True,  -(3*p-2)/4,  "ISM: α (ν > ν_c)",        0.2)
        _run(ism, "ISM", (1e13, 5e14), 1e4, False, -(p-1)/2,    "ISM: β (ν_m < ν < ν_c)",  0.3)
        _run(ism, "ISM", (1e16, 1e18), 1e4, False, -p/2,        "ISM: β (ν > ν_c)",         0.3)

        print("--- Wind ---")
        wind = self._get_model("wind")
        _run(wind, "Wind", (1e3, 1e5), 1e14, True, -(3*p-1)/4, "Wind: α (ν_m < ν < ν_c)", 0.2)

    def run_closure_relations(self):
        print("\n=== 4. Closure Relations ===\n")
        ism = self._get_model("ISM")
        try:
            t_arr, nu_arr = np.logspace(4, 4.7, 20), np.logspace(13, 14, 15)
            alpha = fit_powerlaw(t_arr, ism.flux_density(t_arr, np.full_like(t_arr, 5e13)).total)
            beta = fit_powerlaw(nu_arr, ism.flux_density(np.full_like(nu_arr, 2e4), nu_arr).total)
            expected_alpha = 3 * beta / 2
            passed = bool(abs(alpha - expected_alpha) < 0.35)
        except Exception as e:
            print(f"  ERROR: {e}")
            alpha, beta, expected_alpha, passed = np.nan, np.nan, np.nan, False
        result = {"name": "ISM: α = 3β/2", "type": "closure", "alpha_measured": float(alpha) if not np.isnan(alpha) else None,
                  "beta_measured": float(beta) if not np.isnan(beta) else None, "alpha_expected": float(expected_alpha) if not np.isnan(expected_alpha) else None, "passed": passed}
        self._record(result, "closure_relations", "ISM")
        print(f"  {'PASS' if passed else 'FAIL'}: α={alpha:.3f}, β={beta:.3f}, expected α={expected_alpha:.3f}")

    def run_spectrum_shapes(self):
        print("\n=== 5. Spectrum Shapes (5 Regimes) ===\n")
        self.results.setdefault("spectrum_grid", {})
        for regime_name, config in REGIME_TEST_CONFIGS.items():
            medium, t_test = config["medium"], config["t"]
            mk = _MEDIUM_KEY.get(medium, medium)
            print(f"\n  Regime {regime_name}")
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
                    print(f"    NOTE: Expected regime {regime_name}, found regime {actual}")

                seg_results = _measure_segment_slopes(model, t_test, regime_name, nu_a, nu_m, nu_c, self.p)
                for sr in seg_results:
                    result = {"name": f"Regime {regime_name}: {sr['segment']}", "regime": regime_name, "segment": sr["segment"],
                              "type": "spectrum_shape", **{k: sr[k] for k in ("expected", "expected_expr", "measured", "tolerance", "passed")}}
                    self._record(result, "spectrum_shapes", mk)
                    print(f"    {sr['segment']}: β = {sr['measured']:.2f} (expected {sr['expected']:.2f}) [{'PASS' if sr['passed'] else 'FAIL'}]"
                          if sr['measured'] is not None else f"    {sr['segment']}: measurement failed")

                self.results["spectrum_grid"].setdefault(mk, {})[regime_name] = {
                    "t": t_test, "nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c, "actual_regime": actual,
                    "segments": [{"name": f"Regime {regime_name}: {s['segment']}", **s} for s in seg_results]}
            except Exception as e:
                print(f"    ERROR: {e}")

    # --- Visualization data collection ---

    def _phase_mask(self, t, u, phase, medium):
        """Get mask for a phase, returns None if insufficient points."""
        tr = TIME_RANGES[phase][medium]
        mask = (t >= tr[0]) & (t <= tr[1])
        if phase == "deep_newtonian":
            mask &= (u < 0.1)
        return mask if np.sum(mask) >= 5 else None

    def _collect_viz(self):
        """Collect all visualization data in one pass."""
        print("\n=== Collecting Visualization Data ===\n")
        nu_grid = np.logspace(np.log10(VIZ_FREQ_GRID[0]), np.log10(VIZ_FREQ_GRID[1]), VIZ_FREQ_GRID[2])
        p = self.p

        for medium in ["ISM", "wind"]:
            mk = _MEDIUM_KEY[medium]
            model = self._get_model(medium)
            viz_model = self._get_viz_model(medium)
            t_min, t_max, _ = VIZ_TIME_GRID[medium]
            print(f"--- {mk} ---")

            # Shock dynamics & frequencies from viz_model
            try:
                details = viz_model.details(t_min, t_max)
                t, Gamma, r, B, N_p, Doppler, nu_m, nu_c, nu_a = _fwd(
                    details, "Gamma", "r", "B_comv", "N_p", "Doppler", "nu_m", "nu_c", "nu_a")
                u = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                nu_m, nu_c, nu_a = nu_m * Doppler, nu_c * Doppler, nu_a * Doppler

                # Shock dynamics
                shock_phases = {}
                for phase in TIME_RANGES:
                    if (mask := self._phase_mask(t, u, phase, medium)) is None:
                        continue
                    tp = t[mask]
                    fits = {}
                    for name, arr, key in [("r", r, "r"), ("u", u, "u"), ("B", B, "B"), ("N_p", N_p, "N_p")]:
                        vals = arr[mask]
                        valid = (vals > 0) & np.isfinite(vals)
                        if np.sum(valid) > 5:
                            fits[name] = {"measured": float(fit_powerlaw(tp[valid], vals[valid])),
                                          "expected": SHOCK_SCALINGS.get(phase, {}).get(medium, {}).get(key)}
                    shock_phases[phase] = {"t_range": [float(tp.min()), float(tp.max())], "fits": fits}
                self.results.setdefault("viz_shock_dynamics", {})[mk] = {
                    "t": t.tolist(), "u": u.tolist(), "Gamma": Gamma.tolist(),
                    "r": r.tolist(), "B": B.tolist(), "N_p": N_p.tolist(), "phases": shock_phases}

                # Frequencies
                freq_phases = {}
                for phase in TIME_RANGES:
                    if (mask := self._phase_mask(t, u, phase, medium)) is None:
                        continue
                    tp = t[mask]
                    fits = {}
                    for fname, farr in [("nu_m", nu_m), ("nu_c", nu_c), ("nu_a", nu_a)]:
                        valid = (farr[mask] > 0) & np.isfinite(farr[mask])
                        if np.sum(valid) > 5:
                            fits[fname] = {"measured": float(fit_powerlaw(tp[valid], farr[mask][valid])),
                                           "expected": FREQ_SCALINGS.get(phase, {}).get(medium, {}).get(fname)}
                    freq_phases[phase] = {"t_range": [float(tp.min()), float(tp.max())], "fits": fits}
                self.results.setdefault("viz_frequencies", {})[mk] = {
                    "t": t.tolist(), "Doppler": Doppler.tolist(), "nu_m": nu_m.tolist(),
                    "nu_c": nu_c.tolist(), "nu_a": nu_a.tolist(), "phases": freq_phases}
                print(f"  Shock/freq: {len(t)} points")
            except Exception as e:
                print(f"  Shock/freq ERROR: {e}")

            # Spectra from standard model
            spectra = []
            for t_val in [1e2, 1e3, 1e4, 1e5, 1e6]:
                try:
                    flux = model.flux_density(np.full_like(nu_grid, t_val), nu_grid)
                    spectra.append({"t": t_val, "nu": nu_grid.tolist(), "flux": flux.total.tolist()})
                except Exception:
                    pass
            self.results.setdefault("viz_spectra", {})[mk] = {
                "spectra": spectra, "expected_betas": {"below_nu_a": 2.0, "nu_a_to_nu_m": 1/3,
                                                       "nu_m_to_nu_c": -(p-1)/2, "above_nu_c": -p/2}}

            # Lightcurves
            alpha_expected = {"ISM": {"slow_cooling": -(3*p-3)/4, "above_cooling": -(3*p-2)/4, "post_break": -(3*p-3)/4-0.75},
                              "wind": {"slow_cooling": -(3*p-1)/4, "above_cooling": -(3*p-2)/4}}
            lightcurves = []
            for nu in [1e10, 1e12, 1e14, 1e15, 1e17]:
                try:
                    t_pre = np.logspace(2, np.log10(JET_BREAK_TIME * 0.8), 50)
                    flux_pre = model.flux_density(t_pre, np.full_like(t_pre, nu))
                    lc = {"nu": nu, "pre_break": {"t": t_pre.tolist(), "flux": flux_pre.total.tolist(),
                                                   "alpha_measured": float(fit_powerlaw(t_pre, flux_pre.total))}}
                    if medium == "ISM":
                        t_post = np.logspace(np.log10(JET_BREAK_TIME * 1.5), 7, 30)
                        flux_post = model.flux_density(t_post, np.full_like(t_post, nu))
                        lc["post_break"] = {"t": t_post.tolist(), "flux": flux_post.total.tolist(),
                                            "alpha_measured": float(fit_powerlaw(t_post, flux_post.total))}
                    lightcurves.append(lc)
                except Exception:
                    pass
            self.results.setdefault("viz_lightcurves", {})[mk] = {
                "lightcurves": lightcurves, "jet_break_time": JET_BREAK_TIME, "expected_alphas": alpha_expected[medium]}
            print(f"  Spectra: {len(spectra)}, LCs: {len(lightcurves)}")

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
                print(f"  Regime {regime_name}: collected")
            except Exception as e:
                print(f"  Regime {regime_name}: ERROR - {e}")
        self.results["viz_spectrum_shapes"] = viz_shapes

    # --- Run all + summary ---

    def run_all(self, viz=True):
        print("\n" + "="*60 + "\n       COMPREHENSIVE REGRESSION TEST SUITE\n" + "="*60)
        self.run_shock_dynamics()
        self.run_characteristic_frequencies()
        self.run_flux_regimes()
        self.run_closure_relations()
        self.run_spectrum_shapes()
        if viz:
            self._collect_viz()
        self._print_summary()
        return self.results

    def _print_summary(self):
        print("\n" + "="*60 + "\n                    SUMMARY\n" + "="*60)
        print("\nBy Model:")
        for model, stats in self.results["summary"]["by_model"].items():
            total = stats["pass"] + stats["fail"]
            if total > 0:
                print(f"  {model:6s}: {stats['pass']:2d} pass, {stats['fail']:2d} fail ({stats['pass']/total*100:.0f}%)")
        print("\nBy Category:")
        for cat, stats in self.results["summary"]["by_category"].items():
            total = stats["pass"] + stats["fail"]
            if total > 0:
                print(f"  {cat:20s}: {stats['pass']:2d} pass, {stats['fail']:2d} fail ({stats['pass']/total*100:.0f}%)")
        n_tests = len(self.results["tests"])
        n_pass = sum(1 for t in self.results["tests"] if t.get("passed", False))
        print(f"\nOverall: {n_pass}/{n_tests} passed ({n_pass/n_tests*100:.0f}%)\n" + "="*60)

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
