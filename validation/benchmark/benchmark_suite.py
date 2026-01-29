#!/usr/bin/env python3
"""VegasAfterglow benchmark suite for performance testing and convergence analysis."""

import argparse
import json
import multiprocessing as mp
import os
import subprocess
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime
from itertools import product
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(Path(__file__).parent))
from configs import (create_jet, create_medium, create_observer, create_radiation, get_jet_config, get_medium_config,
                     get_radiation_config, get_all_jet_names, get_radiation_names_no_ssc)

sys.path.insert(0, str(Path(__file__).parent.parent))
from validation.colors import _bold, _green, _red, _cyan, _dim, _bold_red, _header

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ComponentTiming:
    flux_single_ms: float = 0.0
    total_ms: float = 0.0
    stage_breakdown: Dict[str, float] = field(default_factory=dict)

@dataclass
class DimensionConvergence:
    dimension: str
    values: List[float]
    times_ms: List[float]
    errors_by_band: Dict[str, List[float]] = field(default_factory=dict)
    mean_errors_by_band: Dict[str, List[float]] = field(default_factory=dict)
    times_by_band: Dict[str, List[float]] = field(default_factory=dict)
    t_array: List[float] = field(default_factory=list)
    flux_by_band: Dict[str, List[List[float]]] = field(default_factory=dict)

@dataclass
class ConfigResult:
    jet_type: str
    medium: str
    radiation: str
    theta_obs: float
    theta_obs_ratio: float
    spreading: bool
    timing: ComponentTiming
    phi_convergence: Optional[DimensionConvergence] = None
    theta_convergence: Optional[DimensionConvergence] = None
    t_convergence: Optional[DimensionConvergence] = None

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

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

FIDUCIAL_RESOLUTION = (0.3, 0.3, 10)
PHI_VALUES = [0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5]
THETA_VALUES = [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
T_VALUES = [5, 7.5, 10, 12.5, 15, 17.5, 20, 22.5, 25]
CONVERGENCE_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18}
VIEWING_ANGLE_RATIOS = [0, 2, 4]
DIM_INDEX = {"phi": 0, "theta": 1, "t": 2}

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def get_git_commit() -> str:
    try:
        return subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                              capture_output=True, text=True, check=True).stdout.strip()
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
    from validation.visualization.common import get_runtime_build_info
    info = get_runtime_build_info()
    return info.get("Compiler", "unknown"), info.get("Flags", "unknown")

def get_jet_theta_c(jet_name: str) -> float:
    return get_jet_config(jet_name).params.get("theta_c", 0.1)

# ---------------------------------------------------------------------------
# Parallel helpers
# ---------------------------------------------------------------------------

def _init_worker():
    os.environ["OMP_NUM_THREADS"] = "1"

def _benchmark_worker(args: Tuple) -> Optional[Dict]:
    jet, med, rad, theta, spread, iters, ref_res = args
    try:
        b = ComprehensiveBenchmark(iterations=iters, reference_resolution=ref_res)
        return asdict(b.benchmark_config(jet, med, rad, theta, spread, verbose=False))
    except Exception as e:
        print(f"Worker {_bold_red('error')} for {jet}/{med}: {e}")
        return None

# ---------------------------------------------------------------------------
# Main benchmark class
# ---------------------------------------------------------------------------

class ComprehensiveBenchmark:
    def __init__(self, iterations=1, reference_resolution=FIDUCIAL_RESOLUTION):
        self.iterations = iterations
        self.reference_resolution = reference_resolution
        self.results: List[ConfigResult] = []
        self.t_ref, self.nu_ref = 1e5, 4.84e14

    def _create_model(self, jet_name, medium_name, radiation_name, theta_obs, resolution, spreading=False):
        import VegasAfterglow as va
        return va.Model(
            create_jet(get_jet_config(jet_name), spreading=spreading),
            create_medium(get_medium_config(medium_name)),
            create_observer(theta_obs),
            create_radiation(get_radiation_config(radiation_name)),
            resolutions=resolution,
        )

    def _compute_timing(self, jet_name, medium_name, radiation_name, theta_obs, resolution, spreading=False):
        import VegasAfterglow as va
        t = np.logspace(2, 7, 30)
        nu = np.full_like(t, self.nu_ref)
        all_breakdowns = []
        for _ in range(self.iterations):
            model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, resolution, spreading)
            va.Model.profile_reset()
            model.flux_density(t, nu)
            all_breakdowns.append(va.Model.profile_data())
        # Average stage times across iterations
        all_stages = set()
        for bd in all_breakdowns:
            all_stages.update(bd.keys())
        stage_breakdown = {}
        for stage in all_stages:
            vals = [bd.get(stage, 0.0) for bd in all_breakdowns]
            stage_breakdown[stage] = np.mean(vals)
        ms = stage_breakdown.get("total", 0.0)
        return ComponentTiming(flux_single_ms=ms, total_ms=ms, stage_breakdown=stage_breakdown)

    def _run_dimension_convergence(self, jet_name, medium_name, radiation_name, theta_obs,
                                   dimension, values, spreading=False):
        import VegasAfterglow as va
        fiducial = list(FIDUCIAL_RESOLUTION)
        dim_idx = DIM_INDEX[dimension]
        t_lc = np.logspace(2, 7, 50)
        bands = list(CONVERGENCE_BANDS)
        n_vals = len(values)

        # Pre-allocate nu arrays for each band (avoid repeated allocation in inner loop)
        nu_arrays = {b: np.full_like(t_lc, nu) for b, nu in CONVERGENCE_BANDS.items()}

        times_ms = []
        times_by_band = {b: [] for b in bands}
        errors_by_band = {b: [] for b in bands}
        mean_errors_by_band = {b: [] for b in bands}
        flux_by_band = {b: [] for b in bands}

        def _nan_result():
            return DimensionConvergence(
                dimension=dimension, values=values, times_ms=[np.nan]*n_vals,
                errors_by_band={b: [np.nan]*n_vals for b in bands},
                mean_errors_by_band={b: [np.nan]*n_vals for b in bands},
                times_by_band={b: [np.nan]*n_vals for b in bands},
                t_array=t_lc.tolist(), flux_by_band={b: [] for b in bands},
            )

        # Reference at 2x highest value
        ref_res = fiducial.copy()
        ref_res[dim_idx] = values[-1] * 2
        try:
            ref_model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(ref_res), spreading)
            ref_lcs = {b: ref_model.flux_density(t_lc, nu_arr).total.copy()
                       for b, nu_arr in nu_arrays.items()}
        except Exception as e:
            print(f"      {_bold_red('Error')} computing reference: {e}")
            return _nan_result()

        for val in values:
            res = fiducial.copy()
            res[dim_idx] = val
            try:
                model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(res), spreading)

                band_times = []
                for band_name in bands:
                    nu_arr = nu_arrays[band_name]
                    va.Model.profile_reset()
                    flux = model.flux_density(t_lc, nu_arr)
                    elapsed = va.Model.profile_data().get("total", 0.0)
                    times_by_band[band_name].append(elapsed)
                    band_times.append(elapsed)
                    flux_by_band[band_name].append(flux.total.tolist())

                    ref = ref_lcs[band_name]
                    valid = (ref > 0) & np.isfinite(ref) & np.isfinite(flux.total)
                    if np.any(valid):
                        rel = np.abs(flux.total[valid] - ref[valid]) / ref[valid]
                        errors_by_band[band_name].append(float(np.max(rel)))
                        mean_errors_by_band[band_name].append(float(np.mean(rel)))
                    else:
                        errors_by_band[band_name].append(np.nan)
                        mean_errors_by_band[band_name].append(np.nan)
                times_ms.append(np.mean(band_times))
            except Exception as e:
                print(f"      {_bold_red('Error')} at {dimension}={val}: {e}")
                times_ms.append(np.nan)
                for b in bands:
                    times_by_band[b].append(np.nan)
                    errors_by_band[b].append(np.nan)
                    mean_errors_by_band[b].append(np.nan)
                    flux_by_band[b].append([])

        return DimensionConvergence(
            dimension=dimension, values=values, times_ms=times_ms,
            errors_by_band=errors_by_band, mean_errors_by_band=mean_errors_by_band, times_by_band=times_by_band,
            t_array=t_lc.tolist(), flux_by_band=flux_by_band,
        )

    def benchmark_config(self, jet_name, medium_name, radiation_name, theta_obs,
                         spreading=False, verbose=True):
        theta_c = get_jet_theta_c(jet_name)
        ratio = theta_obs / theta_c if theta_c > 0 else 0.0
        log = print if verbose else (lambda *_a, **_k: None)

        log(f"\n  Benchmarking: {_bold(f'{jet_name}/{medium_name}/{radiation_name}')}/θ_v/θ_c={ratio:.1f}")
        timing = self._compute_timing(jet_name, medium_name, radiation_name, theta_obs, self.reference_resolution, spreading)

        result = ConfigResult(
            jet_type=jet_name, medium=medium_name, radiation=radiation_name, theta_obs=theta_obs,
            theta_obs_ratio=ratio, spreading=spreading, timing=timing,
        )

        for dim, vals, attr in [("phi", PHI_VALUES, "phi_convergence"),
                                ("theta", THETA_VALUES, "theta_convergence"),
                                ("t", T_VALUES, "t_convergence")]:
            log(f"    {_cyan(dim)} convergence...")
            setattr(result, attr, self._run_dimension_convergence(
                jet_name, medium_name, radiation_name, theta_obs, dim, vals, spreading))

        log(f"    {_green('Done')}: {_dim(f'{timing.flux_single_ms:.2f} ms/LC')}")
        self.results.append(result)
        return result

    # -----------------------------------------------------------------------
    # Suite runners
    # -----------------------------------------------------------------------

    def _run_suite(self, title, configs):
        print(_header(title))
        print(f"Running {_bold(str(len(configs)))} configurations...")
        for i, (jet, med, rad, theta) in enumerate(configs, 1):
            print(f"\n{_dim(f'[{i}/{len(configs)}]')}", end="")
            self.benchmark_config(jet, med, rad, theta)
        return self.results

    def run_full_suite(self):
        jets, meds, rads = get_all_jet_names(), ["ISM", "wind"], get_radiation_names_no_ssc()
        configs = [(j, m, r, get_jet_theta_c(j) * ratio)
                   for j, m, r, ratio in product(jets, meds, rads, VIEWING_ANGLE_RATIOS)]
        return self._run_suite("FULL BENCHMARK SUITE", configs)

    def run_convergence_only(self, jet_types=None, medium_types=None):
        jets = jet_types or ["tophat", "gaussian", "powerlaw"]
        meds = medium_types or ["ISM", "wind"]
        return self._run_suite("CONVERGENCE TEST SUITE",
            [(j, m, "synchrotron_only", 0.0) for j, m in product(jets, meds)])

    def run_parallel(self, jet_types=None, medium_types=None, radiation_types=None,
                     viewing_ratios=None, n_workers=None, timeout=600):
        print(_header("PARALLEL BENCHMARK SUITE"))
        jets = jet_types or get_all_jet_names()
        meds = medium_types or ["ISM", "wind"]
        rads = radiation_types or ["synchrotron_only"]
        ratios = viewing_ratios or [0]
        n_workers = n_workers or max(1, mp.cpu_count() - 1)

        configs = [(j, m, r, get_jet_theta_c(j) * ratio, False, self.iterations, self.reference_resolution)
                   for j, m, r, ratio in product(jets, meds, rads, ratios)]

        print(f"Running {_bold(str(len(configs)))} configurations with {_bold(str(n_workers))} workers "
              f"(timeout {timeout}s per config)...")
        completed, results_dicts, timed_out = 0, [], 0
        with mp.Pool(n_workers, initializer=_init_worker) as pool:
            async_results = [(cfg, pool.apply_async(_benchmark_worker, (cfg,))) for cfg in configs]
            for cfg, ar in async_results:
                completed += 1
                jet, med, rad = cfg[0], cfg[1], cfg[2]
                try:
                    rd = ar.get(timeout=timeout)
                except mp.TimeoutError:
                    timed_out += 1
                    print(f"  {_dim(f'[{completed}/{len(configs)}]')} {_bold_red('TIMEOUT')} "
                          f"{_bold(f'{jet}/{med}/{rad}')}: exceeded {timeout}s")
                    continue
                if rd:
                    results_dicts.append(rd)
                    ms = rd.get('timing', {}).get('flux_single_ms', 0)
                    print(f"  {_dim(f'[{completed}/{len(configs)}]')} {_bold(f'{jet}/{med}')}: "
                          f"single-freq LC {_dim(f'{ms:.1f} ms')}")
                else:
                    print(f"  {_dim(f'[{completed}/{len(configs)}]')} {_bold_red('FAILED')}")
            pool.terminate()

        status = f"Completed {_bold(f'{len(results_dicts)}/{len(configs)}')} configurations"
        if timed_out:
            status += f" ({_bold_red(f'{timed_out} timed out')})"
        print(f"\n{status}")
        for rd in results_dicts:
            timing = ComponentTiming(**rd.pop("timing"))
            convs = {k: DimensionConvergence(**rd.pop(k)) if rd.get(k) else None
                     for k in ["phi_convergence", "theta_convergence", "t_convergence"]}
            for k in convs: rd.pop(k, None)
            self.results.append(ConfigResult(timing=timing, **convs, **rd))
        return self.results

    # -----------------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------------

    def save_results(self, filepath: str):
        compiler, compile_flags = get_compiler_info()
        session = BenchmarkSession(
            timestamp=datetime.now().isoformat(), commit=get_git_commit(), platform=get_platform_info(),
            python_version=sys.version.split()[0], vegasafterglow_version=get_vegasafterglow_version(),
            compiler=compiler, compile_flags=compile_flags, configs=self.results,
        )
        session_dict = asdict(session)
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(session_dict))
        print(f"\nResults saved to {_dim(filepath)}")

    def print_summary(self):
        if not self.results:
            print("No results to summarize"); return
        print(_header("BENCHMARK SUMMARY"))

        by_jet = {}
        for r in self.results:
            by_jet.setdefault(r.jet_type, []).append(r)

        print(f"\n{_bold('Jet Type'):<28} {_bold('SingleLC (ms)'):<15}")
        print(_dim("-" * 35))
        for jet, results in by_jet.items():
            ms = np.mean([r.timing.flux_single_ms for r in results])
            print(f"{_cyan(jet):<28} {ms:<15.2f}")

        print(f"\n{_bold('Convergence')} {_dim('(mean relative error)')}:")
        print(_dim("-" * 80))
        for r in self.results:
            if r.phi_convergence:
                def _mean_err(c):
                    vals = [v for band in c.mean_errors_by_band.values() for v in band if not np.isnan(v)]
                    return np.mean(vals) if vals else 0.0
                label = _bold(f"{r.jet_type}/{r.medium}")
                phi_e = _mean_err(r.phi_convergence)
                theta_e = _mean_err(r.theta_convergence)
                t_e = _mean_err(r.t_convergence)
                fmt = lambda v: _green(f"{v:.2e}") if v < 0.05 else (_red(f"{v:.2e}") if v > 0.1 else f"{v:.2e}")
                print(f"{label}: phi={fmt(phi_e)}, theta={fmt(theta_e)}, t={fmt(t_e)}")


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Benchmark Suite")
    parser.add_argument("--full", action="store_true", help="Full benchmark (all combinations)")
    parser.add_argument("--convergence", action="store_true", help="Convergence tests only")
    parser.add_argument("-j", "--parallel", type=int, default=0, metavar="N", help="Parallel workers")
    parser.add_argument("--jet", type=str, nargs="+", help="Jet type(s)")
    parser.add_argument("--medium", type=str, nargs="+", help="Medium type(s)")
    parser.add_argument("--iterations", type=int, default=10, help="Timing iterations for overview stage breakdown")
    parser.add_argument("--timeout", type=int, default=600, help="Per-config timeout in seconds (default: 600)")
    parser.add_argument("--output", type=str, default="results/benchmark_history.json", help="Output file")
    args = parser.parse_args()

    print(_header("VegasAfterglow Benchmark Suite"))
    print(f"Commit: {_bold(get_git_commit())}\nPlatform: {get_platform_info()}\nPython: {sys.version.split()[0]}")

    b = ComprehensiveBenchmark(iterations=args.iterations)
    if args.parallel > 0:
        print(f"Parallel workers: {_bold(str(args.parallel))}")
        ratios = VIEWING_ANGLE_RATIOS if args.full else [0]
        rads = get_radiation_names_no_ssc() if args.full else ["synchrotron_only"]
        b.run_parallel(jet_types=args.jet, medium_types=args.medium, radiation_types=rads,
                       viewing_ratios=ratios, n_workers=args.parallel, timeout=args.timeout)
    elif args.full:
        b.run_full_suite()
    elif args.convergence:
        b.run_convergence_only(jet_types=args.jet, medium_types=args.medium)
    else:
        b.run_full_suite()  # Default to full suite

    b.print_summary()
    b.save_results(str(Path(__file__).parent / args.output))


if __name__ == "__main__":
    main()
