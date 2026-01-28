#!/usr/bin/env python3
"""VegasAfterglow benchmark suite for systematic performance testing and convergence analysis."""

import argparse
import gc
import json
import multiprocessing as mp
import os
import subprocess
import sys
import time
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from configs import (create_jet, create_medium, create_observer, create_radiation, get_jet_config, get_medium_config,
                     get_radiation_config, get_all_jet_names, get_all_medium_names, get_all_radiation_names,
                     get_radiation_names_no_ssc)


@dataclass
class ComponentTiming:
    flux_single_ms: float = 0.0
    total_ms: float = 0.0


@dataclass
class ResolutionPoint:
    resolution: Tuple[float, float, float]
    flux_value: float
    time_ms: float


@dataclass
class DimensionConvergence:
    dimension: str
    fiducial_resolution: Tuple[float, float, float]
    values: List[float]
    times_ms: List[float]
    errors_by_band: Dict[str, List[float]] = field(default_factory=dict)
    mean_errors_by_band: Dict[str, List[float]] = field(default_factory=dict)
    times_by_band: Dict[str, List[float]] = field(default_factory=dict)
    t_array: List[float] = field(default_factory=list)
    flux_by_band: Dict[str, List[List[float]]] = field(default_factory=dict)
    relative_errors: List[float] = field(default_factory=list)


@dataclass
class ConfigResult:
    jet_type: str
    medium: str
    radiation: str
    theta_obs: float
    theta_obs_ratio: float
    spreading: bool
    timing: ComponentTiming
    lightcurve_t: List[float] = field(default_factory=list)
    lightcurve_flux: List[float] = field(default_factory=list)
    spectrum_nu: List[float] = field(default_factory=list)
    spectrum_flux: List[float] = field(default_factory=list)
    phi_convergence: Optional[DimensionConvergence] = None
    theta_convergence: Optional[DimensionConvergence] = None
    t_convergence: Optional[DimensionConvergence] = None
    resolution_costs: List[ResolutionPoint] = field(default_factory=list)


@dataclass
class BenchmarkSession:
    timestamp: str
    commit: str
    platform: str
    python_version: str
    vegasafterglow_version: str = "unknown"
    compiler: str = "unknown"
    compile_flags: str = "unknown"
    configs: List[ConfigResult] = field(default_factory=list)


FIDUCIAL_RESOLUTION = (0.3, 0.3, 10)
PHI_VALUES = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
THETA_VALUES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
T_VALUES = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25]
CONVERGENCE_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18}
VIEWING_ANGLE_RATIOS = [0, 2, 4]


def get_git_commit() -> str:
    try:
        result = subprocess.run(["git", "rev-parse", "--short", "HEAD"], capture_output=True, text=True, check=True)
        return result.stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        return "unknown"


def get_platform_info() -> str:
    import platform
    return f"{platform.system()} {platform.machine()}"


def get_vegasafterglow_version() -> str:
    try:
        import VegasAfterglow
        return getattr(VegasAfterglow, "__version__", "unknown")
    except ImportError:
        return "not installed"


def get_compiler_info() -> Tuple[str, str]:
    compiler, flags = "unknown", "unknown"
    compile_commands = Path(__file__).parent.parent.parent / "compile_commands.json"
    if compile_commands.exists():
        try:
            with open(compile_commands, "r") as f:
                data = json.load(f)
            if data:
                cmd = data[0].get("command", "")
                parts = cmd.split()
                if parts:
                    compiler = Path(parts[0]).name
                    flag_set = set()
                    for p in parts[1:]:
                        if p.startswith("-O") or p.startswith("-march") or p.startswith("-std"):
                            flag_set.add(p)
                        elif p in ["-fopenmp", "-ffast-math", "-DNDEBUG"]:
                            flag_set.add(p)
                    if flag_set:
                        flags = " ".join(sorted(flag_set))
        except (json.JSONDecodeError, KeyError, IndexError):
            pass
    if compiler != "unknown":
        try:
            result = subprocess.run([compiler, "--version"], capture_output=True, text=True, check=True)
            compiler = result.stdout.split("\n")[0].strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass
    return compiler, flags


def time_function_min(func: Callable, iterations: int = 1) -> float:
    """Time a function and return minimum time in ms (filters OS scheduling noise)."""
    gc.collect()
    gc.disable()
    times = []
    try:
        for _ in range(iterations):
            start = time.perf_counter_ns()
            func()
            times.append((time.perf_counter_ns() - start) / 1e6)
    finally:
        gc.enable()
    return float(np.min(times))


def get_jet_theta_c(jet_name: str) -> float:
    return get_jet_config(jet_name).params.get("theta_c", 0.1)


def _benchmark_worker(args: Tuple) -> Optional[Dict]:
    jet_name, medium_name, radiation_name, theta_obs, spreading, run_convergence, iterations, reference_resolution = args
    try:
        benchmark = ComprehensiveBenchmark(iterations=iterations, reference_resolution=reference_resolution)
        result = benchmark.benchmark_config(jet_name, medium_name, radiation_name, theta_obs, spreading,
                                            run_convergence, verbose=False)
        return asdict(result)
    except Exception as e:
        print(f"Worker error for {jet_name}/{medium_name}: {e}")
        return None


def _init_worker():
    os.environ["OMP_NUM_THREADS"] = "1"


class ComprehensiveBenchmark:
    def __init__(self, iterations: int = 1, reference_resolution: Tuple[float, float, float] = FIDUCIAL_RESOLUTION):
        self.iterations = iterations
        self.reference_resolution = reference_resolution
        self.results: List[ConfigResult] = []
        self.t_ref = 1e5
        self.nu_ref = 4.84e14

    def _create_model(self, jet_name: str, medium_name: str, radiation_name: str, theta_obs: float,
                      resolution: Tuple[float, float, float], spreading: bool = False):
        import VegasAfterglow as va
        jet = create_jet(get_jet_config(jet_name), spreading=spreading)
        medium = create_medium(get_medium_config(medium_name))
        observer = create_observer(theta_obs)
        radiation = create_radiation(get_radiation_config(radiation_name))
        return va.Model(jet, medium, observer, radiation, resolutions=resolution)

    def _compute_timing(self, jet_name: str, medium_name: str, radiation_name: str, theta_obs: float,
                        resolution: Tuple[float, float, float], spreading: bool = False) -> ComponentTiming:
        model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, resolution, spreading)
        model.details(t_min=1e2, t_max=1e7)
        t = np.logspace(2, 7, 30)
        nu = np.full_like(t, self.nu_ref)
        model.flux_density(t, nu)  # warmup
        flux_time = time_function_min(lambda: model.flux_density(t, nu), self.iterations)
        return ComponentTiming(flux_single_ms=flux_time, total_ms=flux_time)

    def _compute_lightcurve(self, model, n_points: int = 50) -> Tuple[List[float], List[float]]:
        t = np.logspace(2, 7, n_points)
        nu = np.full_like(t, self.nu_ref)
        flux = model.flux_density(t, nu)
        return t.tolist(), flux.total.tolist()

    def _compute_spectrum(self, model, n_points: int = 50) -> Tuple[List[float], List[float]]:
        nu = np.logspace(9, 20, n_points)
        t = np.full_like(nu, self.t_ref)
        flux = model.flux_density(t, nu)
        return nu.tolist(), flux.total.tolist()

    def _run_dimension_convergence(self, jet_name: str, medium_name: str, radiation_name: str, theta_obs: float,
                                   dimension: str, values: List[float], spreading: bool = False) -> DimensionConvergence:
        fiducial = list(FIDUCIAL_RESOLUTION)
        dim_idx = {"phi": 0, "theta": 1, "t": 2}[dimension]

        times_ms, times_by_band = [], {band: [] for band in CONVERGENCE_BANDS}
        errors_by_band = {band: [] for band in CONVERGENCE_BANDS}
        mean_errors_by_band = {band: [] for band in CONVERGENCE_BANDS}
        flux_by_band = {band: [] for band in CONVERGENCE_BANDS}

        n_lc_points = 150
        t_lc = np.logspace(2, 7, n_lc_points)

        ref_res = fiducial.copy()
        ref_res[dim_idx] = values[-1] * 2

        try:
            ref_model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(ref_res), spreading)
            ref_model.details(t_min=1e2, t_max=1e7)
            ref_lcs = {}
            for band_name, nu in CONVERGENCE_BANDS.items():
                nu_arr = np.full_like(t_lc, nu)
                ref_lcs[band_name] = ref_model.flux_density(t_lc, nu_arr).total.copy()
        except Exception as e:
            print(f"      Error computing reference: {e}")
            return DimensionConvergence(
                dimension=dimension, fiducial_resolution=tuple(fiducial), values=values,
                times_ms=[np.nan] * len(values),
                errors_by_band={band: [np.nan] * len(values) for band in CONVERGENCE_BANDS},
                mean_errors_by_band={band: [np.nan] * len(values) for band in CONVERGENCE_BANDS},
                times_by_band={band: [np.nan] * len(values) for band in CONVERGENCE_BANDS},
                t_array=t_lc.tolist(), flux_by_band={band: [] for band in CONVERGENCE_BANDS},
            )

        for val in values:
            res = fiducial.copy()
            res[dim_idx] = val

            try:
                model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(res), spreading)
                model.details(t_min=1e2, t_max=1e7)
                for band_name, nu in CONVERGENCE_BANDS.items():
                    model.flux_density(t_lc, np.full_like(t_lc, nu))  # warmup

                band_times = []
                n_timing_iters = self.iterations if theta_obs == 0 else 1
                gc.collect()
                gc.disable()
                try:
                    for band_name, nu in CONVERGENCE_BANDS.items():
                        nu_arr = np.full_like(t_lc, nu)
                        iter_times = []
                        for _ in range(n_timing_iters):
                            start = time.perf_counter_ns()
                            flux = model.flux_density(t_lc, nu_arr)
                            iter_times.append((time.perf_counter_ns() - start) / 1e6)
                        elapsed = np.min(iter_times) if n_timing_iters > 1 else iter_times[0]
                        times_by_band[band_name].append(elapsed)
                        band_times.append(elapsed)
                        flux_by_band[band_name].append(flux.total.tolist())

                        ref_flux = ref_lcs[band_name]
                        valid = (ref_flux > 0) & np.isfinite(ref_flux) & np.isfinite(flux.total)
                        if np.any(valid):
                            rel_err = np.abs(flux.total[valid] - ref_flux[valid]) / ref_flux[valid]
                            max_err, mean_err = float(np.max(rel_err)), float(np.mean(rel_err))
                        else:
                            max_err, mean_err = np.nan, np.nan
                        errors_by_band[band_name].append(max_err)
                        mean_errors_by_band[band_name].append(mean_err)
                finally:
                    gc.enable()
                times_ms.append(np.mean(band_times))
            except Exception as e:
                print(f"      Error at {dimension}={val}: {e}")
                times_ms.append(np.nan)
                for band_name in CONVERGENCE_BANDS:
                    times_by_band[band_name].append(np.nan)
                    errors_by_band[band_name].append(np.nan)
                    mean_errors_by_band[band_name].append(np.nan)
                    flux_by_band[band_name].append([])

        all_errors = np.array([errors_by_band[b] for b in CONVERGENCE_BANDS])
        relative_errors = np.nanmean(all_errors, axis=0).tolist()

        return DimensionConvergence(
            dimension=dimension, fiducial_resolution=tuple(fiducial), values=values, times_ms=times_ms,
            errors_by_band=errors_by_band, mean_errors_by_band=mean_errors_by_band, times_by_band=times_by_band,
            t_array=t_lc.tolist(), flux_by_band=flux_by_band, relative_errors=relative_errors,
        )

    def benchmark_config(self, jet_name: str, medium_name: str, radiation_name: str, theta_obs: float,
                         spreading: bool = False, run_convergence: bool = True, verbose: bool = True) -> ConfigResult:
        theta_c = get_jet_theta_c(jet_name)
        theta_obs_ratio = theta_obs / theta_c if theta_c > 0 else 0.0
        config_str = f"{jet_name}/{medium_name}/{radiation_name}/θ_v/θ_c={theta_obs_ratio:.1f}"
        if spreading:
            config_str += "/spread"

        if verbose:
            print(f"\n  Benchmarking: {config_str}")
            print("    Computing timing...")

        timing = self._compute_timing(jet_name, medium_name, radiation_name, theta_obs, self.reference_resolution, spreading)
        model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, self.reference_resolution, spreading)

        if verbose:
            print("    Computing light curve...")
        lc_t, lc_flux = self._compute_lightcurve(model)
        if verbose:
            print("    Computing spectrum...")
        spec_nu, spec_flux = self._compute_spectrum(model)

        result = ConfigResult(
            jet_type=jet_name, medium=medium_name, radiation=radiation_name, theta_obs=theta_obs,
            theta_obs_ratio=theta_obs_ratio, spreading=spreading, timing=timing,
            lightcurve_t=lc_t, lightcurve_flux=lc_flux, spectrum_nu=spec_nu, spectrum_flux=spec_flux,
        )

        if run_convergence:
            if verbose:
                print("    Running phi convergence...")
            result.phi_convergence = self._run_dimension_convergence(jet_name, medium_name, radiation_name, theta_obs,
                                                                     "phi", PHI_VALUES, spreading)
            if verbose:
                print("    Running theta convergence...")
            result.theta_convergence = self._run_dimension_convergence(jet_name, medium_name, radiation_name, theta_obs,
                                                                       "theta", THETA_VALUES, spreading)
            if verbose:
                print("    Running t convergence...")
            result.t_convergence = self._run_dimension_convergence(jet_name, medium_name, radiation_name, theta_obs,
                                                                   "t", T_VALUES, spreading)
            if verbose:
                print("    Computing resolution costs...")
            for scale in [0.5, 0.75, 1.0, 1.5, 2.0]:
                res = tuple(r * scale for r in self.reference_resolution)
                try:
                    m = self._create_model(jet_name, medium_name, radiation_name, theta_obs, res, spreading)
                    elapsed = time_function_min(lambda: m.flux_density(np.array([self.t_ref]), np.array([self.nu_ref])),
                                                self.iterations)
                    flux = m.flux_density(np.array([self.t_ref]), np.array([self.nu_ref]))
                    result.resolution_costs.append(ResolutionPoint(resolution=res, flux_value=float(flux.total[0]),
                                                                   time_ms=elapsed))
                except Exception as e:
                    print(f"      Error at scale {scale}: {e}")

        if verbose:
            print(f"    Single-freq LC time: {timing.flux_single_ms:.2f} ms")
        self.results.append(result)
        return result

    def run_quick_suite(self) -> List[ConfigResult]:
        print("\n" + "=" * 70 + "\nQUICK BENCHMARK SUITE\n" + "=" * 70)
        for jet, medium, rad, theta in [("tophat", "ISM", "synchrotron_only", 0.0),
                                        ("gaussian", "ISM", "synchrotron_only", 0.0)]:
            self.benchmark_config(jet, medium, rad, theta, spreading=False, run_convergence=False)
        return self.results

    def run_standard_suite(self, jet_types: Optional[List[str]] = None, medium_types: Optional[List[str]] = None,
                           radiation_types: Optional[List[str]] = None, viewing_angles: Optional[List[float]] = None,
                           spreading: bool = False, run_convergence: bool = True) -> List[ConfigResult]:
        print("\n" + "=" * 70 + "\nSTANDARD BENCHMARK SUITE\n" + "=" * 70)
        jet_types = jet_types or get_all_jet_names()
        medium_types = medium_types or ["ISM", "wind"]
        radiation_types = radiation_types or ["synchrotron_only"]
        viewing_angles = viewing_angles or [0.0]

        total = len(jet_types) * len(medium_types) * len(radiation_types) * len(viewing_angles)
        print(f"Running {total} configurations...")

        count = 0
        for jet_name in jet_types:
            for medium_name in medium_types:
                for radiation_name in radiation_types:
                    for theta_obs in viewing_angles:
                        count += 1
                        print(f"\n[{count}/{total}]", end="")
                        self.benchmark_config(jet_name, medium_name, radiation_name, theta_obs,
                                              spreading=spreading, run_convergence=run_convergence)
        return self.results

    def run_full_suite(self) -> List[ConfigResult]:
        print("\n" + "=" * 70 + "\nFULL BENCHMARK SUITE - ALL COMBINATIONS\n" + "=" * 70)
        jet_types = get_all_jet_names()
        medium_types = ["ISM", "wind"]
        radiation_types = get_radiation_names_no_ssc()

        total = len(jet_types) * len(medium_types) * len(radiation_types) * len(VIEWING_ANGLE_RATIOS)
        print(f"Running {total} configurations (θ_v/θ_c = {VIEWING_ANGLE_RATIOS})...")

        count = 0
        for jet_name in jet_types:
            theta_c = get_jet_theta_c(jet_name)
            for medium_name in medium_types:
                for radiation_name in radiation_types:
                    for ratio in VIEWING_ANGLE_RATIOS:
                        count += 1
                        print(f"\n[{count}/{total}]", end="")
                        self.benchmark_config(jet_name, medium_name, radiation_name, theta_c * ratio,
                                              spreading=False, run_convergence=True)
        return self.results

    def run_convergence_only(self, jet_types: Optional[List[str]] = None,
                             medium_types: Optional[List[str]] = None) -> List[ConfigResult]:
        print("\n" + "=" * 70 + "\nCONVERGENCE TEST SUITE\n" + "=" * 70)
        jet_types = jet_types or ["tophat", "gaussian", "powerlaw"]
        medium_types = medium_types or ["ISM", "wind"]
        for jet_name in jet_types:
            for medium_name in medium_types:
                self.benchmark_config(jet_name, medium_name, "synchrotron_only", 0.0, spreading=False, run_convergence=True)
        return self.results

    def run_parallel(self, jet_types: Optional[List[str]] = None, medium_types: Optional[List[str]] = None,
                     radiation_types: Optional[List[str]] = None, viewing_angles: Optional[List[float]] = None,
                     viewing_ratios: Optional[List[float]] = None, spreading: bool = False,
                     run_convergence: bool = True, n_workers: Optional[int] = None) -> List[ConfigResult]:
        print("\n" + "=" * 70 + "\nPARALLEL BENCHMARK SUITE\n" + "=" * 70)
        jet_types = jet_types or get_all_jet_names()
        medium_types = medium_types or ["ISM", "wind"]
        radiation_types = radiation_types or ["synchrotron_only"]
        n_workers = n_workers or max(1, mp.cpu_count() - 1)

        configs_to_run = []
        for jet_name in jet_types:
            theta_c = get_jet_theta_c(jet_name)
            if viewing_ratios is not None:
                jet_viewing_angles = [theta_c * ratio for ratio in viewing_ratios]
            elif viewing_angles is not None:
                jet_viewing_angles = viewing_angles
            else:
                jet_viewing_angles = [0.0]

            for medium_name in medium_types:
                for radiation_name in radiation_types:
                    for theta_obs in jet_viewing_angles:
                        configs_to_run.append((jet_name, medium_name, radiation_name, theta_obs, spreading,
                                               run_convergence, self.iterations, self.reference_resolution))

        total = len(configs_to_run)
        if viewing_ratios is not None:
            print(f"Running {total} configurations with {n_workers} workers (θ_v/θ_c = {viewing_ratios})...")
        else:
            print(f"Running {total} configurations with {n_workers} workers...")

        completed, results_dicts = 0, []
        with mp.Pool(n_workers, initializer=_init_worker) as pool:
            for result_dict in pool.imap_unordered(_benchmark_worker, configs_to_run):
                completed += 1
                if result_dict is not None:
                    results_dicts.append(result_dict)
                    jet = result_dict.get("jet_type", "?")
                    medium = result_dict.get("medium", "?")
                    lc_time = result_dict.get("timing", {}).get("flux_single_ms", 0)
                    print(f"  [{completed}/{total}] {jet}/{medium}: single-freq LC {lc_time:.1f} ms")
                else:
                    print(f"  [{completed}/{total}] FAILED")

        print(f"\nCompleted {len(results_dicts)}/{total} configurations")

        for rd in results_dicts:
            timing = ComponentTiming(**rd.pop("timing"))
            phi_conv = DimensionConvergence(**rd.pop("phi_convergence")) if rd.get("phi_convergence") else None
            theta_conv = DimensionConvergence(**rd.pop("theta_convergence")) if rd.get("theta_convergence") else None
            t_conv = DimensionConvergence(**rd.pop("t_convergence")) if rd.get("t_convergence") else None
            rd.pop("phi_convergence", None)
            rd.pop("theta_convergence", None)
            rd.pop("t_convergence", None)
            res_costs = [ResolutionPoint(**rc) for rc in rd.pop("resolution_costs", [])]
            self.results.append(ConfigResult(timing=timing, phi_convergence=phi_conv, theta_convergence=theta_conv,
                                             t_convergence=t_conv, resolution_costs=res_costs, **rd))
        return self.results

    def save_results(self, filepath: str):
        compiler, compile_flags = get_compiler_info()
        session = BenchmarkSession(
            timestamp=datetime.now().isoformat(), commit=get_git_commit(), platform=get_platform_info(),
            python_version=sys.version.split()[0], vegasafterglow_version=get_vegasafterglow_version(),
            compiler=compiler, compile_flags=compile_flags, configs=self.results,
        )
        history_path = Path(filepath)
        if history_path.exists():
            with open(history_path, "r") as f:
                history = json.load(f)
        else:
            history = {"sessions": []}
        history["sessions"].append(asdict(session))
        history_path.parent.mkdir(parents=True, exist_ok=True)
        with open(history_path, "w") as f:
            json.dump(history, f, indent=2)
        print(f"\nResults saved to {filepath}")

    def print_summary(self):
        if not self.results:
            print("No results to summarize")
            return

        print("\n" + "=" * 70 + "\nBENCHMARK SUMMARY\n" + "=" * 70)
        by_jet = {}
        for r in self.results:
            by_jet.setdefault(r.jet_type, []).append(r)

        print("\nTiming by Jet Type (ms):")
        print("-" * 50)
        print(f"{'Jet Type':<20} {'SingleLC':<15}")
        print("-" * 50)

        for jet_type, results in by_jet.items():
            avg_flux = np.mean([r.timing.flux_single_ms for r in results])
            print(f"{jet_type:<20} {avg_flux:<15.2f}")

        print("\nConvergence Summary (max relative error):")
        print("-" * 80)
        for r in self.results:
            if r.phi_convergence:
                phi_err = max([e for e in r.phi_convergence.relative_errors if not np.isnan(e)], default=0)
                theta_err = max([e for e in r.theta_convergence.relative_errors if not np.isnan(e)], default=0)
                t_err = max([e for e in r.t_convergence.relative_errors if not np.isnan(e)], default=0)
                print(f"{r.jet_type}/{r.medium}: phi={phi_err:.2e}, theta={theta_err:.2e}, t={t_err:.2e}")


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Comprehensive Benchmark Suite",
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quick", action="store_true", help="Run minimal benchmark")
    parser.add_argument("--full", action="store_true", help="Run full benchmark (all combinations)")
    parser.add_argument("--convergence", action="store_true", help="Run convergence tests only")
    parser.add_argument("-j", "--parallel", type=int, default=0, metavar="N", help="Run with N parallel workers")
    parser.add_argument("--jet", type=str, nargs="+", help=f"Jet type(s): {get_all_jet_names()}")
    parser.add_argument("--medium", type=str, nargs="+", help=f"Medium type(s): {get_all_medium_names()}")
    parser.add_argument("--radiation", type=str, nargs="+", help=f"Radiation config(s): {get_all_radiation_names()}")
    parser.add_argument("--viewing-angle", type=float, nargs="+", help="Viewing angle(s) in radians")
    parser.add_argument("--iterations", type=int, default=1, help="Number of timing iterations")
    parser.add_argument("--output", type=str, default="results/benchmark_history.json", help="Output file path")
    args = parser.parse_args()

    print("=" * 70 + "\nVegasAfterglow Comprehensive Benchmark Suite\n" + "=" * 70)
    print(f"Commit: {get_git_commit()}")
    print(f"Platform: {get_platform_info()}")
    print(f"Python: {sys.version.split()[0]}")
    if args.parallel > 0:
        print(f"Parallel workers: {args.parallel}")

    benchmark = ComprehensiveBenchmark(iterations=args.iterations)

    if args.quick:
        benchmark.run_quick_suite()
    elif args.parallel > 0:
        if args.full:
            benchmark.run_parallel(jet_types=get_all_jet_names(), medium_types=["ISM", "wind"],
                                   radiation_types=get_radiation_names_no_ssc(), viewing_ratios=VIEWING_ANGLE_RATIOS,
                                   spreading=False, run_convergence=True, n_workers=args.parallel)
        else:
            benchmark.run_parallel(jet_types=args.jet, medium_types=args.medium, radiation_types=args.radiation,
                                   viewing_angles=args.viewing_angle, spreading=False, run_convergence=True,
                                   n_workers=args.parallel)
    elif args.full:
        benchmark.run_full_suite()
    elif args.convergence:
        benchmark.run_convergence_only(jet_types=args.jet, medium_types=args.medium)
    else:
        benchmark.run_standard_suite(jet_types=args.jet, medium_types=args.medium, radiation_types=args.radiation,
                                     viewing_angles=args.viewing_angle, run_convergence=True)

    benchmark.print_summary()
    benchmark.save_results(str(Path(__file__).parent / args.output))
    print("\nTo generate visualizations, run:")
    print("  python tests/visualization/dashboard.py --full")


if __name__ == "__main__":
    main()
