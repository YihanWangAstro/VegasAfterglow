#!/usr/bin/env python3
"""Comprehensive regression test runner for VegasAfterglow power-law scaling verification."""

import argparse
import json
import sys
from datetime import datetime
from fractions import Fraction as F
from pathlib import Path

import numpy as np
from VegasAfterglow import ISM, Model, Observer, Radiation, TophatJet, Wind

# Allow running from project root or tests/regression directory
sys.path.insert(0, str(Path(__file__).parent.parent.parent))
from validation.regression.utils import STANDARD_PARAMS, fit_powerlaw

TIME_RANGES = {
    "coasting": {"ISM": (1e-1, 1), "wind": (1e-2, 1e-1)},
    "BM": {"ISM": (5e2, 5e3), "wind": (5e3, 5e4)},
    "deep_newtonian": {"ISM": (1e12, 1e13), "wind": (1e14, 1e15)},
}

JET_BREAK_TIME = 1.8e5
VIZ_TIME_GRID = {"ISM": (1e-2, 1e13, 400), "wind": (1e-3, 1e15, 400)}
VIZ_FREQ_GRID = (1e8, 1e20, 100)
VIZ_PARAMS = {"ISM": {"Gamma0": 300, "n_ism": 0.1}, "wind": {"Gamma0": 50, "A_star": 1e-1}}

SHOCK_SCALINGS = {
    "coasting": {"ISM": {"u": F(0), "r": F(1), "B": F(0), "N_p": F(3)}, "wind": {"u": F(0), "r": F(1), "B": F(-1), "N_p": F(1)}},
    "BM": {"ISM": {"u": F(-3, 8), "r": F(1, 4), "B": F(-3, 8), "N_p": F(3, 4)}, "wind": {"u": F(-1, 4), "r": F(1, 2), "B": F(-3, 4), "N_p": F(1, 2)}},
    "deep_newtonian": {"ISM": {"u": F(-3, 5), "r": F(2, 5), "B": F(-3, 5), "N_p": F(6, 5)}, "wind": {"u": F(-1, 3), "r": F(2, 3), "B": F(-1), "N_p": F(2, 3)}},
}

FREQ_SCALINGS = {
    "coasting": {"ISM": {"nu_m": F(0), "nu_c": F(-2)}, "wind": {"nu_m": F(-1), "nu_c": F(-1)}},
    "BM": {"ISM": {"nu_m": F(-3, 2), "nu_c": F(-1, 2)}, "wind": {"nu_m": F(-3, 2), "nu_c": F(1, 2)}},
    "deep_newtonian": {"ISM": {"nu_m": F(-3, 5), "nu_c": F(-1, 5)}, "wind": {"nu_m": F(-1), "nu_c": F(1)}},
}

SLOPE_TOLERANCE = 0.1
TOLERANCES = {"u": SLOPE_TOLERANCE, "Gamma": SLOPE_TOLERANCE, "r": SLOPE_TOLERANCE, "B": SLOPE_TOLERANCE,
              "N_p": SLOPE_TOLERANCE, "nu_m": SLOPE_TOLERANCE, "nu_c": SLOPE_TOLERANCE}

# Spectral regimes based on frequency ordering (ν_a, ν_m, ν_c)
# Each regime has specific spectral indices (β where F_ν ∝ ν^β) in different frequency ranges
SPECTRAL_REGIMES = {
    "I": {  # ν_a < ν_m < ν_c (slow cooling, optically thin)
        "order": ("nu_a", "nu_m", "nu_c"),
        "segments": [
            ("below_nu_a", 2.0),           # self-absorbed
            ("nu_a_to_nu_m", F(1, 3)),     # optically thin, rising
            ("nu_m_to_nu_c", "-(p-1)/2"),  # power-law electrons
            ("above_nu_c", "-p/2"),        # cooled electrons
        ],
    },
    "II": {  # ν_m < ν_a < ν_c (slow cooling, self-absorption above ν_m)
        "order": ("nu_m", "nu_a", "nu_c"),
        "segments": [
            ("below_nu_m", 2.0),
            ("nu_m_to_nu_a", F(5, 2)),     # self-absorbed power-law
            ("nu_a_to_nu_c", "-(p-1)/2"),
            ("above_nu_c", "-p/2"),
        ],
    },
    "III": {  # ν_a < ν_c < ν_m (fast cooling)
        "order": ("nu_a", "nu_c", "nu_m"),
        "segments": [
            ("below_nu_a", 2.0),
            ("nu_a_to_nu_c", F(1, 3)),
            ("nu_c_to_nu_m", F(-1, 2)),    # cooled electrons
            ("above_nu_m", "-p/2"),
        ],
    },
    "IV": {  # ν_c < ν_a < ν_m (fast cooling, self-absorption modified)
        "order": ("nu_c", "nu_a", "nu_m"),
        "segments": [
            ("below_nu_c", 2),
            ("nu_c_to_nu_a", F(2)),    # self-absorbed cooling
            ("nu_a_to_nu_m", F(-1, 2)),
            ("above_nu_m", "-p/2"),
        ],
    },
    "V": {  # ν_c < ν_m < ν_a (fast cooling, strong self-absorption)
        "order": ("nu_c", "nu_m", "nu_a"),
        "segments": [
            ("below_nu_c", 2.0),
            ("nu_c_to_nu_m", 2.0),    # cooled electrons
            ("nu_m_to_nu_a", F(5, 2)),     # self-absorbed power-law
            ("above_nu_a", "-p/2"),
        ],
    },
}

SPECTRAL_TOLERANCE = 0.15

REGIME_TEST_CONFIGS = {
    "I": {  # ν_a < ν_m < ν_c (slow cooling, standard)
        "medium": "ISM",
        "t": 5e3,
        "description": "ISM slow cooling (ν_a < ν_m < ν_c)",
        "n_ism": 0.1,
        "eps_B": 1e-2,
        # Standard parameters work well for this regime
    },
    "II": {  # ν_m < ν_a < ν_c (slow cooling, self-absorption above ν_m)
        "medium": "ISM",
        "t": 5e4,
        "description": "ISM late time (ν_m < ν_a < ν_c)",
        "n_ism": 1e6,
        "eps_B": 1e-5,
        "eps_e": 5e-2,
        # At late times, ν_m drops below ν_a
    },
    "III": {  # ν_a < ν_c < ν_m (fast cooling)
        "medium": "ISM",
        "t": 1e5,
        "description": "ISM early fast cooling (ν_a < ν_c < ν_m)",
        "n_ism": 30,
        "eps_B": 1e-1,
        "eps_e": 1e-1,
        "xi_e":5e-3
    },
    "IV": {  # ν_c < ν_a < ν_m (fast cooling, self-absorption modified)
        "medium": "ISM",
        "t": 1e5,
        "description": "ISM fast cooling with absorption (ν_c < ν_a < ν_m)",
        "eps_B": 0.3,
        "eps_e": 0.3,
        "n_ism": 1e6,
        "xi_e": 5e-3,
    },
    "V": {  # ν_c < ν_m < ν_a (fast cooling, strong self-absorption)
        "medium": "ISM",
        "t": 1e4,
        "description": "ISM fast cooling with strong absorption (ν_c < ν_m < ν_a)",
        "n_ism": 1e8,
        "eps_e": 0.01,
        "eps_B": 3e-1,
        "xi_e": 1,
    },
}


def create_model(medium_type, Gamma0=None, A_star=None, eps_e=None, eps_B=None, xi_e=None, n_ism=None):
    """Create a model with customizable parameters."""
    Gamma0 = Gamma0 if Gamma0 is not None else STANDARD_PARAMS["Gamma0"]
    A_star = A_star if A_star is not None else STANDARD_PARAMS["A_star"]
    eps_e = eps_e if eps_e is not None else STANDARD_PARAMS["eps_e"]
    eps_B = eps_B if eps_B is not None else STANDARD_PARAMS["eps_B"]
    xi_e = xi_e if xi_e is not None else STANDARD_PARAMS["xi_e"]
    n_ism = n_ism if n_ism is not None else STANDARD_PARAMS["n_ism"]

    jet = TophatJet(theta_c=STANDARD_PARAMS["theta_c"], E_iso=STANDARD_PARAMS["E_iso"], Gamma0=Gamma0)
    medium = ISM(n_ism=n_ism) if medium_type == "ISM" else Wind(A_star=A_star)
    observer = Observer(lumi_dist=STANDARD_PARAMS["lumi_dist"], z=STANDARD_PARAMS["z"], theta_obs=0.0)
    radiation = Radiation(eps_e=eps_e, eps_B=eps_B, p=STANDARD_PARAMS["p"], xi_e=xi_e)
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
        if medium_type not in self._models:
            self._models[medium_type] = create_model(medium_type)
        return self._models[medium_type]

    def _get_viz_model(self, medium_type):
        key = f"viz_{medium_type}"
        if key not in self._models:
            params = VIZ_PARAMS.get(medium_type, {})
            self._models[key] = create_model(medium_type, Gamma0=params.get("Gamma0"), A_star=params.get("A_star"))
        return self._models[key]

    def _record_result(self, result, category, medium):
        self.results["tests"].append(result)
        self.results["categories"].setdefault(category, []).append(result)
        model_key = "ISM" if medium == "ISM" else "Wind"
        stat_key = "pass" if result.get("passed", False) else "fail"
        self.results["summary"]["by_model"][model_key][stat_key] += 1
        self.results["summary"]["by_category"].setdefault(category, {"pass": 0, "fail": 0})[stat_key] += 1

    def run_shock_dynamics(self):
        print("\n=== 1. Shock Dynamics (All 3 Phases) ===\n")
        if "shock_grid" not in self.results:
            self.results["shock_grid"] = {}
        for medium in ["ISM", "wind"]:
            model = self._get_viz_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"\n--- {medium_key} ---")
            self.results["shock_grid"].setdefault(medium_key, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {phase}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], SHOCK_SCALINGS[phase][medium], {}
                for qty in ["u", "r", "B", "N_p"]:
                    if (exp_val := expected.get(qty)) is None:
                        continue
                    result = self._run_shock_quantity(model, t_range, qty, exp_val, f"{medium_key} {phase}: {qty}", phase)
                    self._record_result(result, "shock_dynamics", medium_key)
                    phase_results[qty] = result
                self.results["shock_grid"][medium_key][phase] = phase_results

    def _run_shock_quantity(self, model, t_range, qty, expected, name, phase="BM"):
        print(f"    Testing: {qty}")
        try:
            details = model.details(t_range[0], t_range[1])
            t = np.array(details.fwd.t_obs)[0, 0, :]
            sort_idx = np.argsort(t)
            t = t[sort_idx]
            time_mask = (t >= t_range[0]) & (t <= t_range[1])
            t = t[time_mask]
            if qty == "u":
                Gamma = np.array(details.fwd.Gamma)[0, 0, :][sort_idx][time_mask]
                y = Gamma * np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                valid = (y > 0) & np.isfinite(y)
                if np.sum(valid) < 5:
                    raise ValueError(f"Not enough valid points ({np.sum(valid)})")
                measured = fit_powerlaw(t[valid], y[valid])
            elif qty == "r":
                measured = fit_powerlaw(t, np.array(details.fwd.r)[0, 0, :][sort_idx][time_mask])
            elif qty == "B":
                measured = fit_powerlaw(t, np.array(details.fwd.B_comv)[0, 0, :][sort_idx][time_mask])
            elif qty == "N_p":
                measured = fit_powerlaw(t, np.array(details.fwd.N_p)[0, 0, :][sort_idx][time_mask])
            else:
                raise ValueError(f"Unknown quantity: {qty}")
        except Exception as e:
            print(f"      ERROR: {e}")
            measured = np.nan
        tolerance = TOLERANCES.get(qty, 0.25)
        passed = bool(not np.isnan(measured) and abs(measured - expected) < tolerance)
        result = {"name": name, "quantity": qty, "phase": phase, "type": f"shock_{qty}", "expected": float(expected),
                  "measured": float(measured) if not np.isnan(measured) else None, "tolerance": tolerance, "passed": passed}
        status = "PASS" if passed else "FAIL"
        print(f"      {status}: measured={measured:.3f}, expected={expected:.3f}" if not np.isnan(measured) else f"      {status}: measurement failed")
        return result

    def run_characteristic_frequencies(self):
        print("\n=== 2. Characteristic Frequencies (All 3 Phases) ===\n")
        if "freq_grid" not in self.results:
            self.results["freq_grid"] = {}
        for medium in ["ISM", "wind"]:
            model = self._get_viz_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"\n--- {medium_key} ---")
            self.results["freq_grid"].setdefault(medium_key, {})
            for phase in ["coasting", "BM", "deep_newtonian"]:
                print(f"\n  Phase: {phase}")
                t_range, expected, phase_results = TIME_RANGES[phase][medium], FREQ_SCALINGS[phase][medium], {}
                for freq_name in ["nu_m", "nu_c"]:
                    if (exp_val := expected.get(freq_name)) is None:
                        continue
                    result = self._run_nu_scaling(model, t_range, freq_name, exp_val, f"{medium_key} {phase}: {freq_name}", phase)
                    self._record_result(result, "frequencies", medium_key)
                    phase_results[freq_name] = result
                self.results["freq_grid"][medium_key][phase] = phase_results

    def _run_nu_scaling(self, model, t_range, freq_name, expected, name, phase="BM"):
        print(f"    Testing: {freq_name}")
        try:
            details = model.details(t_range[0], t_range[1])
            t = np.array(details.fwd.t_obs)[0, 0, :]
            Doppler = np.array(details.fwd.Doppler)[0, 0, :]
            nu = np.array(getattr(details.fwd, freq_name))[0, 0, :] * Doppler
            sort_idx = np.argsort(t)
            t, nu = t[sort_idx], nu[sort_idx]
            time_mask = (t >= t_range[0]) & (t <= t_range[1])
            t, nu = t[time_mask], nu[time_mask]
            valid = (nu > 0) & np.isfinite(nu)
            if np.sum(valid) < 5:
                raise ValueError(f"Not enough valid points ({np.sum(valid)})")
            nu_alpha = fit_powerlaw(t[valid], nu[valid])
        except Exception as e:
            print(f"      ERROR: {e}")
            nu_alpha = np.nan
        tolerance = TOLERANCES.get(freq_name, 0.25)
        passed = bool(not np.isnan(nu_alpha) and abs(nu_alpha - expected) < tolerance)
        result = {"name": name, "quantity": freq_name, "phase": phase, "type": "frequency", "expected": float(expected),
                  "measured": float(nu_alpha) if not np.isnan(nu_alpha) else None, "tolerance": tolerance, "passed": passed}
        status = "PASS" if passed else "FAIL"
        print(f"      {status}: measured={nu_alpha:.3f}, expected={expected:.3f}" if not np.isnan(nu_alpha) else f"      {status}: measurement failed")
        return result

    def run_flux_regimes(self):
        print("\n=== 3. Flux Scalings (Regime I) ===\n")
        ism_model = self._get_model("ISM")
        print("--- ISM ---")

        def run_temporal(nu, t_range, expected, name, medium):
            print(f"  Testing: {name}")
            t = np.logspace(np.log10(t_range[0]), np.log10(t_range[1]), 30)
            try:
                flux = ism_model.flux_density(t, np.full_like(t, nu)) if medium == "ISM" else self._get_model("wind").flux_density(t, np.full_like(t, nu))
                measured = fit_powerlaw(t, flux.total)
            except Exception as e:
                print(f"    ERROR: {e}")
                measured = np.nan
            passed = bool(not np.isnan(measured) and abs(measured - expected) < 0.2 * abs(expected))
            result = {"name": name, "type": "temporal", "expected": expected, "measured": float(measured) if not np.isnan(measured) else None, "tolerance": 0.2, "passed": passed}
            self.results["tests"].append(result)
            self.results["scaling_results"].append(result)
            self._record_result(result, "flux_regimes", medium)
            print(f"    {'PASS' if passed else 'FAIL'}: measured={measured:.3f}, expected={expected:.3f}")

        def run_spectral(t_fixed, nu_range, expected, name, medium):
            print(f"  Testing: {name}")
            nu = np.logspace(np.log10(nu_range[0]), np.log10(nu_range[1]), 30)
            try:
                flux = ism_model.flux_density(np.full_like(nu, t_fixed), nu)
                measured = fit_powerlaw(nu, flux.total)
            except Exception as e:
                print(f"    ERROR: {e}")
                measured = np.nan
            passed = bool(not np.isnan(measured) and abs(measured - expected) < 0.3)
            result = {"name": name, "type": "spectral", "expected": expected, "measured": float(measured) if not np.isnan(measured) else None, "tolerance": 0.3, "passed": passed}
            self.results["tests"].append(result)
            self.results["scaling_results"].append(result)
            self._record_result(result, "flux_regimes", medium)
            print(f"    {'PASS' if passed else 'FAIL'}: measured={measured:.3f}, expected={expected:.3f}")

        run_temporal(5e13, (1e4, 5e4), -(3 * self.p - 3) / 4, "ISM: α (ν_m < ν < ν_c)", "ISM")
        run_temporal(1e17, (1e3, 5e4), -(3 * self.p - 2) / 4, "ISM: α (ν > ν_c)", "ISM")
        run_spectral(1e4, (1e13, 5e14), -(self.p - 1) / 2, "ISM: β (ν_m < ν < ν_c)", "ISM")
        run_spectral(1e4, (1e16, 1e18), -self.p / 2, "ISM: β (ν > ν_c)", "ISM")

        print("--- Wind ---")
        wind_model = self._get_model("wind")
        t = np.logspace(3, 5, 30)
        try:
            flux = wind_model.flux_density(t, np.full_like(t, 1e14))
            measured = fit_powerlaw(t, flux.total)
        except Exception:
            measured = np.nan
        expected = -(3 * self.p - 1) / 4
        passed = bool(not np.isnan(measured) and abs(measured - expected) < 0.2 * abs(expected))
        result = {"name": "Wind: α (ν_m < ν < ν_c)", "type": "temporal", "expected": expected, "measured": float(measured) if not np.isnan(measured) else None, "tolerance": 0.2, "passed": passed}
        self.results["tests"].append(result)
        self._record_result(result, "flux_regimes", "Wind")
        print(f"  Testing: Wind: α (ν_m < ν < ν_c)")
        print(f"    {'PASS' if passed else 'FAIL'}: measured={measured:.3f}, expected={expected:.3f}")

    def run_closure_relations(self):
        print("\n=== 4. Closure Relations ===\n")
        ism_model = self._get_model("ISM")
        t_array, nu_array = np.logspace(4, 4.7, 20), np.logspace(13, 14, 15)
        try:
            flux_t = ism_model.flux_density(t_array, np.full_like(t_array, 5e13))
            alpha = fit_powerlaw(t_array, flux_t.total)
            flux_nu = ism_model.flux_density(np.full_like(nu_array, 2e4), nu_array)
            beta = fit_powerlaw(nu_array, flux_nu.total)
            expected_alpha = 3 * beta / 2
            passed = bool(abs(alpha - expected_alpha) < 0.35)
        except Exception as e:
            print(f"  ERROR: {e}")
            alpha, beta, expected_alpha, passed = np.nan, np.nan, np.nan, False
        result = {"name": "ISM: α = 3β/2", "type": "closure", "alpha_measured": float(alpha) if not np.isnan(alpha) else None,
                  "beta_measured": float(beta) if not np.isnan(beta) else None, "alpha_expected": float(expected_alpha) if not np.isnan(expected_alpha) else None, "passed": passed}
        self._record_result(result, "closure_relations", "ISM")
        print(f"  {'PASS' if passed else 'FAIL'}: α={alpha:.3f}, β={beta:.3f}, expected α={expected_alpha:.3f}")

    def _create_regime_model(self, config):
        """Create a model with custom parameters for a specific regime test."""
        medium = config.get("medium", "ISM")
        return create_model(
            medium,
            eps_e=config.get("eps_e"),
            eps_B=config.get("eps_B"),
            xi_e=config.get("xi_e"),
            n_ism=config.get("n_ism"),
            A_star=config.get("A_star"),
        )

    def run_spectrum_shapes(self):
        """Test spectral shapes for 5 regimes based on ν_a, ν_m, ν_c ordering."""
        print("\n=== 5. Spectrum Shapes (5 Regimes) ===\n")
        if "spectrum_grid" not in self.results:
            self.results["spectrum_grid"] = {}

        for regime_name, config in REGIME_TEST_CONFIGS.items():
            medium = config.get("medium", "ISM")
            t_test = config["t"]
            description = config.get("description", f"Regime {regime_name}")

            print(f"\n  Regime {regime_name}: {description}")

            # Create model with custom parameters for this regime
            model = self._create_regime_model(config)
            medium_key = "ISM" if medium == "ISM" else "Wind"

            try:
                # Get characteristic frequencies at this time
                details = model.details(t_test * 0.9, t_test * 1.1)
                idx = np.argmin(np.abs(np.array(details.fwd.t_obs)[0, 0, :] - t_test))
                Doppler = np.array(details.fwd.Doppler)[0, 0, idx]
                nu_a = np.array(details.fwd.nu_a)[0, 0, idx] * Doppler
                nu_m = np.array(details.fwd.nu_m)[0, 0, idx] * Doppler
                nu_c = np.array(details.fwd.nu_c)[0, 0, idx] * Doppler

                print(f"    ν_a = {nu_a:.2e}, ν_m = {nu_m:.2e}, ν_c = {nu_c:.2e}")

                # Determine actual regime from frequency ordering
                freqs = {"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}
                sorted_freqs = sorted(freqs.items(), key=lambda x: x[1])
                actual_order = tuple(f[0] for f in sorted_freqs)
                actual_regime = None
                for r_name, r_info in SPECTRAL_REGIMES.items():
                    if r_info["order"] == actual_order:
                        actual_regime = r_name
                        break

                if actual_regime is None:
                    print(f"    WARNING: Unexpected frequency order {actual_order}")
                    actual_regime = regime_name  # Fall back to expected

                if actual_regime != regime_name:
                    print(f"    NOTE: Expected regime {regime_name}, found regime {actual_regime}")

                # Test spectral slopes using expected regime's segments
                # (even if actual regime differs, we test what was intended)
                regime_info = SPECTRAL_REGIMES[regime_name]
                segment_results = []

                for seg_name, expected_beta in regime_info["segments"]:
                    # Determine frequency range for this segment
                    nu_low, nu_high = self._get_segment_range(seg_name, nu_a, nu_m, nu_c)
                    if nu_low is None or nu_high is None:
                        continue
                    if nu_high <= nu_low * 1.5:
                        print(f"    {seg_name}: skipped (range too narrow)")
                        continue

                    # Compute expected value and store original expression
                    if isinstance(expected_beta, F):
                        expected_val = float(expected_beta)
                        expected_expr = f"{expected_beta.numerator}/{expected_beta.denominator}" if expected_beta.denominator != 1 else str(expected_beta.numerator)
                    elif isinstance(expected_beta, str):
                        expected_val = eval(expected_beta.replace("p", str(self.p)))
                        expected_expr = expected_beta  # e.g., "-(p-1)/2"
                    else:
                        expected_val = float(expected_beta)
                        expected_expr = str(int(expected_beta)) if expected_beta == int(expected_beta) else str(expected_beta)

                    # Measure spectral slope
                    nu_grid = np.logspace(np.log10(nu_low * 1.2), np.log10(nu_high * 0.8), 20)
                    flux = model.flux_density(np.full_like(nu_grid, t_test), nu_grid)
                    measured = fit_powerlaw(nu_grid, flux.total)

                    passed = bool(not np.isnan(measured) and abs(measured - expected_val) < SPECTRAL_TOLERANCE)
                    status = "PASS" if passed else "FAIL"

                    print(f"    {seg_name}: β = {measured:.2f} (expected {expected_val:.2f}) [{status}]")

                    result = {
                        "name": f"Regime {regime_name}: {seg_name}",
                        "regime": regime_name,
                        "segment": seg_name,
                        "type": "spectrum_shape",
                        "expected": expected_val,
                        "expected_expr": expected_expr,
                        "measured": float(measured) if not np.isnan(measured) else None,
                        "tolerance": SPECTRAL_TOLERANCE,
                        "passed": passed,
                    }
                    self._record_result(result, "spectrum_shapes", medium_key)
                    segment_results.append(result)

                self.results["spectrum_grid"].setdefault(medium_key, {})[regime_name] = {
                    "t": t_test,
                    "nu_a": nu_a,
                    "nu_m": nu_m,
                    "nu_c": nu_c,
                    "actual_regime": actual_regime,
                    "segments": segment_results,
                }

            except Exception as e:
                print(f"    ERROR: {e}")

    def _get_segment_range(self, seg_name, nu_a, nu_m, nu_c):
        """Get frequency range for a spectral segment."""
        freq_map = {"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}
        parts = seg_name.split("_to_")

        if seg_name.startswith("below_"):
            nu_ref = freq_map.get(seg_name.replace("below_", ""))
            return (nu_ref / 1000, nu_ref / 30) if nu_ref else (None, None)
        elif seg_name.startswith("above_"):
            nu_ref = freq_map.get(seg_name.replace("above_", ""))
            return (nu_ref * 10, nu_ref * 100) if nu_ref else (None, None)
        elif len(parts) == 2:
            nu_low = freq_map.get(parts[0])
            nu_high = freq_map.get(parts[1])
            return (nu_low * 10, nu_high / 10) if nu_low and nu_high else (None, None)
        return None, None

    def collect_shock_dynamics_viz(self):
        print("\n=== Collecting Shock Dynamics Visualization Data ===\n")
        viz_data = {}
        for medium in ["ISM", "wind"]:
            model = self._get_viz_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"--- {medium_key} ---")
            t_min, t_max, _ = VIZ_TIME_GRID[medium]
            try:
                details = model.details(t_min, t_max)
                t = np.array(details.fwd.t_obs)[0, 0, :]
                Gamma = np.array(details.fwd.Gamma)[0, 0, :]
                r = np.array(details.fwd.r)[0, 0, :]
                B = np.array(details.fwd.B_comv)[0, 0, :]
                N_p = np.array(details.fwd.N_p)[0, 0, :]
                beta = np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                u = Gamma * beta
                sort_idx = np.argsort(t)
                t, Gamma, r, B, u, N_p = t[sort_idx], Gamma[sort_idx], r[sort_idx], B[sort_idx], u[sort_idx], N_p[sort_idx]
                phases = {}
                for phase, t_range in TIME_RANGES.items():
                    phase_range = t_range[medium]
                    mask = (t >= phase_range[0]) & (t <= phase_range[1])
                    if phase == "deep_newtonian":
                        u_mask = u < 0.1
                        if not np.any(u_mask) or np.sum(mask & u_mask) < 5:
                            continue
                        mask = mask & u_mask
                    elif np.sum(mask) < 5:
                        continue
                    t_phase = t[mask]
                    fits = {"r": {"measured": float(fit_powerlaw(t_phase, r[mask])), "expected": SHOCK_SCALINGS.get(phase, {}).get(medium, {}).get("r")}}
                    u_phase = u[mask]
                    if np.sum((u_phase > 0) & np.isfinite(u_phase)) > 5:
                        fits["u"] = {"measured": float(fit_powerlaw(t_phase[(u_phase > 0) & np.isfinite(u_phase)], u_phase[(u_phase > 0) & np.isfinite(u_phase)])), "expected": SHOCK_SCALINGS.get(phase, {}).get(medium, {}).get("u")}
                    fits["B"] = {"measured": float(fit_powerlaw(t_phase, B[mask])), "expected": SHOCK_SCALINGS.get(phase, {}).get(medium, {}).get("B")}
                    N_p_phase = N_p[mask]
                    if np.sum((N_p_phase > 0) & np.isfinite(N_p_phase)) > 5:
                        fits["N_p"] = {"measured": float(fit_powerlaw(t_phase[(N_p_phase > 0) & np.isfinite(N_p_phase)], N_p_phase[(N_p_phase > 0) & np.isfinite(N_p_phase)])), "expected": SHOCK_SCALINGS.get(phase, {}).get(medium, {}).get("N_p")}
                    phases[phase] = {"t_range": [float(t_phase.min()), float(t_phase.max())], "fits": fits}
                viz_data[medium_key] = {"t": t.tolist(), "u": u.tolist(), "Gamma": Gamma.tolist(), "r": r.tolist(), "B": B.tolist(), "N_p": N_p.tolist(), "phases": phases}
                print(f"  Collected {len(t)} time points")
            except Exception as e:
                print(f"  ERROR: {e}")
        self.results["viz_shock_dynamics"] = viz_data
        return viz_data

    def collect_frequencies_viz(self):
        print("\n=== Collecting Characteristic Frequencies Visualization Data ===\n")
        viz_data = {}
        for medium in ["ISM", "wind"]:
            model = self._get_viz_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"--- {medium_key} ---")
            t_min, t_max, _ = VIZ_TIME_GRID[medium]
            try:
                details = model.details(t_min, t_max)
                t = np.array(details.fwd.t_obs)[0, 0, :]
                Doppler = np.array(details.fwd.Doppler)[0, 0, :]
                Gamma = np.array(details.fwd.Gamma)[0, 0, :]
                beta = np.sqrt(1.0 - 1.0 / (Gamma * Gamma))
                u = Gamma * beta
                nu_m = np.array(details.fwd.nu_m)[0, 0, :] * Doppler
                nu_c = np.array(details.fwd.nu_c)[0, 0, :] * Doppler
                nu_a = np.array(details.fwd.nu_a)[0, 0, :] * Doppler
                try:
                    nu_M = np.array(details.fwd.nu_M)[0, 0, :] * Doppler
                except (AttributeError, KeyError):
                    try:
                        gamma_max = np.array(details.fwd.gamma_max)[0, 0, :]
                        B_comv = np.array(details.fwd.B_comv)[0, 0, :]
                        nu_M = 4.2e6 * B_comv * gamma_max**2 * Doppler
                    except (AttributeError, KeyError):
                        nu_M = np.array([])
                sort_idx = np.argsort(t)
                t, Doppler, u = t[sort_idx], Doppler[sort_idx], u[sort_idx]
                nu_m, nu_c, nu_a = nu_m[sort_idx], nu_c[sort_idx], nu_a[sort_idx]
                if len(nu_M) > 0:
                    nu_M = nu_M[sort_idx]
                phases = {}
                for phase, phase_t_range in TIME_RANGES.items():
                    t_range = phase_t_range[medium]
                    mask = (t >= t_range[0]) & (t <= t_range[1])
                    if phase == "deep_newtonian":
                        u_mask = u < 0.1
                        if not np.any(u_mask) or np.sum(mask & u_mask) < 5:
                            continue
                        mask = mask & u_mask
                    elif np.sum(mask) < 5:
                        continue
                    t_phase = t[mask]
                    fits = {}
                    for fname, farr in [("nu_m", nu_m), ("nu_c", nu_c), ("nu_a", nu_a)]:
                        valid = (farr[mask] > 0) & np.isfinite(farr[mask])
                        if np.sum(valid) > 5:
                            fits[fname] = {"measured": float(fit_powerlaw(t_phase[valid], farr[mask][valid])), "expected": FREQ_SCALINGS.get(phase, {}).get(medium, {}).get(fname)}
                    if len(nu_M) > 0:
                        nu_M_masked = nu_M[mask]
                        valid = (nu_M_masked > 0) & np.isfinite(nu_M_masked)
                        if np.sum(valid) > 5:
                            fits["nu_M"] = {"measured": float(fit_powerlaw(t_phase[valid], nu_M_masked[valid])), "expected": None}
                    phases[phase] = {"t_range": [float(t_phase.min()), float(t_phase.max())], "fits": fits}
                viz_data[medium_key] = {"t": t.tolist(), "Doppler": Doppler.tolist(), "nu_m": nu_m.tolist(), "nu_c": nu_c.tolist(),
                                        "nu_a": nu_a.tolist(), "nu_M": nu_M.tolist() if len(nu_M) > 0 else [], "phases": phases}
                print(f"  Collected {len(t)} time points, {len(phases)} phases")
            except Exception as e:
                print(f"  ERROR: {e}")
        self.results["viz_frequencies"] = viz_data
        return viz_data

    def collect_spectrum_viz(self):
        print("\n=== Collecting Spectrum Visualization Data ===\n")
        viz_data = {}
        sample_times = [1e2, 1e3, 1e4, 1e5, 1e6]
        nu_min, nu_max, n_nu = VIZ_FREQ_GRID
        nu_grid = np.logspace(np.log10(nu_min), np.log10(nu_max), n_nu)
        expected_betas = {"below_nu_a": 2.0, "nu_a_to_nu_m": 1/3, "nu_m_to_nu_c": -(self.p-1)/2, "above_nu_c": -self.p/2}
        for medium in ["ISM", "wind"]:
            model = self._get_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"--- {medium_key} ---")
            spectra = []
            for t in sample_times:
                try:
                    flux = model.flux_density(np.full_like(nu_grid, t), nu_grid)
                    spectra.append({"t": t, "nu": nu_grid.tolist(), "flux": flux.total.tolist()})
                    print(f"  t = {t:.0e} s: collected")
                except Exception as e:
                    print(f"  t = {t:.0e} s: ERROR - {e}")
            viz_data[medium_key] = {"spectra": spectra, "expected_betas": expected_betas}
        self.results["viz_spectra"] = viz_data
        return viz_data

    def collect_lightcurve_viz(self):
        print("\n=== Collecting Light Curve Visualization Data ===\n")
        viz_data = {}
        sample_freqs = [1e10, 1e12, 1e14, 1e15, 1e17]
        expected_alphas_ism = {"slow_cooling": -(3*self.p - 3) / 4, "above_cooling": -(3*self.p - 2) / 4, "post_break": -(3*self.p - 3) / 4 - 0.75}
        expected_alphas_wind = {"slow_cooling": -(3*self.p - 1) / 4, "above_cooling": -(3*self.p - 2) / 4}
        for medium in ["ISM", "wind"]:
            model = self._get_model(medium)
            medium_key = "ISM" if medium == "ISM" else "Wind"
            print(f"--- {medium_key} ---")
            lightcurves = []
            for nu in sample_freqs:
                try:
                    t_pre = np.logspace(2, np.log10(JET_BREAK_TIME * 0.8), 50)
                    flux_pre = model.flux_density(t_pre, np.full_like(t_pre, nu))
                    alpha_pre = fit_powerlaw(t_pre, flux_pre.total)
                    lc_data = {"nu": nu, "pre_break": {"t": t_pre.tolist(), "flux": flux_pre.total.tolist(), "alpha_measured": float(alpha_pre)}}
                    if medium == "ISM":
                        t_post = np.logspace(np.log10(JET_BREAK_TIME * 1.5), 7, 30)
                        flux_post = model.flux_density(t_post, np.full_like(t_post, nu))
                        lc_data["post_break"] = {"t": t_post.tolist(), "flux": flux_post.total.tolist(), "alpha_measured": float(fit_powerlaw(t_post, flux_post.total))}
                    lightcurves.append(lc_data)
                    print(f"  ν = {nu:.0e} Hz: α_pre = {alpha_pre:.2f}")
                except Exception as e:
                    print(f"  ν = {nu:.0e} Hz: ERROR - {e}")
            viz_data[medium_key] = {"lightcurves": lightcurves, "jet_break_time": JET_BREAK_TIME,
                                    "expected_alphas": expected_alphas_ism if medium == "ISM" else expected_alphas_wind}
        self.results["viz_lightcurves"] = viz_data
        return viz_data

    def collect_spectrum_shapes_viz(self):
        """Collect visualization data for spectrum shape tests."""
        print("\n=== Collecting Spectrum Shapes Visualization Data ===\n")
        viz_data = {}

        for regime_name, config in REGIME_TEST_CONFIGS.items():
            medium = config.get("medium", "ISM")
            t_test = config["t"]
            medium_key = "ISM" if medium == "ISM" else "Wind"

            print(f"  Regime {regime_name} ({medium_key}, t={t_test:.0e}s)")

            # Create model with custom parameters for this regime
            model = self._create_regime_model(config)

            try:
                # Get characteristic frequencies
                details = model.details(t_test * 0.9, t_test * 1.1)
                idx = np.argmin(np.abs(np.array(details.fwd.t_obs)[0, 0, :] - t_test))
                Doppler = np.array(details.fwd.Doppler)[0, 0, idx]
                nu_a = np.array(details.fwd.nu_a)[0, 0, idx] * Doppler
                nu_m = np.array(details.fwd.nu_m)[0, 0, idx] * Doppler
                nu_c = np.array(details.fwd.nu_c)[0, 0, idx] * Doppler

                # Compute spectrum with fixed frequency range
                nu_grid = np.logspace(6, 20, 200)  # 1e6 to 1e20 Hz
                flux = model.flux_density(np.full_like(nu_grid, t_test), nu_grid)

                # Determine actual regime
                freqs = {"nu_a": nu_a, "nu_m": nu_m, "nu_c": nu_c}
                sorted_freqs = sorted(freqs.items(), key=lambda x: x[1])
                actual_order = tuple(f[0] for f in sorted_freqs)
                actual_regime = regime_name
                for r_name, r_info in SPECTRAL_REGIMES.items():
                    if r_info["order"] == actual_order:
                        actual_regime = r_name
                        break

                # Compute measured slopes for expected regime's segments
                regime_info = SPECTRAL_REGIMES.get(regime_name, {})
                segment_results = []
                for seg_name, expected_beta in regime_info.get("segments", []):
                    nu_low, nu_high = self._get_segment_range(seg_name, nu_a, nu_m, nu_c)
                    if nu_low is None or nu_high is None or nu_high <= nu_low * 1.5:
                        continue

                    # Compute expected value
                    if isinstance(expected_beta, str):
                        expected_val = eval(expected_beta.replace("p", str(self.p)))
                    else:
                        expected_val = float(expected_beta)

                    # Measure slope in this segment
                    seg_nu = np.logspace(np.log10(nu_low * 1.2), np.log10(nu_high * 0.8), 20)
                    seg_flux = model.flux_density(np.full_like(seg_nu, t_test), seg_nu)
                    measured = fit_powerlaw(seg_nu, seg_flux.total)

                    # Store original expression for display
                    if isinstance(expected_beta, F):
                        expected_expr = f"{expected_beta.numerator}/{expected_beta.denominator}" if expected_beta.denominator != 1 else str(expected_beta.numerator)
                    elif isinstance(expected_beta, str):
                        expected_expr = expected_beta  # e.g., "-(p-1)/2"
                    else:
                        expected_expr = str(expected_beta)

                    segment_results.append({
                        "segment": seg_name,
                        "expected": expected_val,
                        "expected_expr": expected_expr,
                        "measured": float(measured) if not np.isnan(measured) else None,
                        "passed": bool(not np.isnan(measured) and abs(measured - expected_val) < SPECTRAL_TOLERANCE),
                        "nu_low": float(nu_low),
                        "nu_high": float(nu_high),
                    })

                viz_data[f"regime_{regime_name}"] = {
                    "medium": medium_key,
                    "t": t_test,
                    "nu_a": nu_a,
                    "nu_m": nu_m,
                    "nu_c": nu_c,
                    "actual_regime": actual_regime,
                    "nu": nu_grid.tolist(),
                    "flux": flux.total.tolist(),
                    "expected_segments": segment_results,
                }
                print(f"    Collected (actual regime: {actual_regime})")

            except Exception as e:
                print(f"    ERROR: {e}")

        self.results["viz_spectrum_shapes"] = viz_data
        return viz_data

    def run_all(self, quick=False, viz=True):
        print("\n" + "="*60 + "\n       COMPREHENSIVE REGRESSION TEST SUITE\n" + "="*60)
        self.run_shock_dynamics()
        self.run_characteristic_frequencies()
        self.run_flux_regimes()
        if not quick:
            self.run_closure_relations()
            self.run_spectrum_shapes()
        if viz:
            self.collect_shock_dynamics_viz()
            self.collect_frequencies_viz()
            self.collect_spectrum_viz()
            self.collect_lightcurve_viz()
            self.collect_spectrum_shapes_viz()
        self._print_summary()
        return self.results

    def _print_summary(self):
        print("\n" + "="*60 + "\n                    SUMMARY\n" + "="*60)
        print("\nBy Model:")
        for model, stats in self.results["summary"]["by_model"].items():
            total = stats["pass"] + stats["fail"]
            if total > 0:
                print(f"  {model:6s}: {stats['pass']:2d} pass, {stats['fail']:2d} fail ({stats['pass'] / total * 100:.0f}%)")
        print("\nBy Category:")
        for cat, stats in self.results["summary"]["by_category"].items():
            total = stats["pass"] + stats["fail"]
            if total > 0:
                print(f"  {cat:20s}: {stats['pass']:2d} pass, {stats['fail']:2d} fail ({stats['pass'] / total * 100:.0f}%)")
        n_tests = len(self.results["tests"])
        n_pass = sum(1 for t in self.results["tests"] if t.get("passed", False))
        print(f"\nOverall: {n_pass}/{n_tests} passed ({n_pass/n_tests*100:.0f}%)\n" + "="*60)

    def save_results(self, filename="regression_results.json"):
        output_path = self.output_dir / filename
        with open(output_path, "w") as f:
            json.dump(self.results, f, indent=2, default=str)
        print(f"\nResults saved to {output_path}")
        return output_path


def main():
    parser = argparse.ArgumentParser(description="Run VegasAfterglow regression tests")
    parser.add_argument("--quick", action="store_true", help="Run quick subset of tests")
    parser.add_argument("--no-viz", action="store_true", help="Skip visualization data collection")
    parser.add_argument("--output", default="results", help="Output directory for results")
    args = parser.parse_args()
    runner = RegressionRunner(output_dir=args.output)
    runner.run_all(quick=args.quick, viz=not args.no_viz)
    runner.save_results()


if __name__ == "__main__":
    main()
