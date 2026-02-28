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
from configs import (create_jet, create_medium, create_observer, create_radiation, create_rvs_radiation,
                     get_jet_config, get_medium_config, get_radiation_config, get_duration,
                     get_all_jet_names, get_fwd_only_radiation_names, get_rvs_radiation_names)

sys.path.insert(0, str(Path(__file__).parent.parent))
from validation.colors import _bold, _green, _red, _cyan, _dim, _bold_red, _header
from validation.visualization.common import (FIDUCIAL_RESOLUTION, DIM_INDEX, FIDUCIAL_VALUES,
                                             get_git_commit)

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
    ref_flux_by_band: Dict[str, List[float]] = field(default_factory=dict)

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
    cpu: str = "unknown"
    configs: List[ConfigResult] = field(default_factory=list)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

PHI_VALUES = [0.1, 0.15, 0.2, 0.25]
THETA_VALUES = [0.5, 1.0, 1.5, 2.0]
T_VALUES = [5, 10, 15, 20]
REFERENCE_RESOLUTION = (PHI_VALUES[-1] * 1.2, THETA_VALUES[-1] * 1.2, T_VALUES[-1] * 1.2)
CONVERGENCE_T_LC = np.logspace(0, 8, 100)
CONVERGENCE_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18}
SSC_BANDS = {"Radio": 1e9, "Optical": 4.84e14, "X-ray": 1e18, "TeV": 2.4e26}  # TeV ~ 1 TeV
VIEWING_ANGLE_RATIOS = [0, 2, 4]

# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

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

def get_cpu_info() -> str:
    from validation.visualization.common import get_runtime_build_info
    info = get_runtime_build_info()
    return info.get("CPU", "unknown")

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
        rad_config = get_radiation_config(radiation_name)
        rvs_rad = create_rvs_radiation(rad_config)
        duration = get_duration(rad_config)
        return va.Model(
            create_jet(get_jet_config(jet_name), spreading=spreading, duration=duration),
            create_medium(get_medium_config(medium_name)),
            create_observer(theta_obs),
            create_radiation(rad_config),
            rvs_rad=rvs_rad,
            resolutions=resolution,
        )

    def _compute_timing(self, jet_name, medium_name, radiation_name, theta_obs, resolution, spreading=False):
        import VegasAfterglow as va
        t = np.logspace(2, 8, 30)
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

    @staticmethod
    def _extract_flux(flux, has_rvs, has_ssc=False):
        """Extract the single flux array to compare against reference."""
        if has_rvs:
            return np.asarray(flux.rvs.sync)
        if has_ssc:
            return np.asarray(flux.fwd.ssc)
        return np.asarray(flux.total)

    def _run_dimension_convergence(self, jet_name, medium_name, radiation_name, theta_obs,
                                   dimension, values, spreading=False):
        import VegasAfterglow as va
        fiducial = list(FIDUCIAL_RESOLUTION)
        dim_idx = DIM_INDEX[dimension]
        t_lc = CONVERGENCE_T_LC
        n_vals = len(values)
        rad_config = get_radiation_config(radiation_name)
        has_rvs = rad_config.rvs_params is not None
        has_ssc = rad_config.params.get("ssc", False)

        # Use extended bands (including TeV) for SSC configs
        bands_dict = SSC_BANDS if has_ssc else CONVERGENCE_BANDS
        bands = list(bands_dict)

        # Pre-allocate nu arrays for each frequency band
        nu_arrays = {b: np.full_like(t_lc, nu) for b, nu in bands_dict.items()}

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

        # Compute per-dimension reference: 2x the highest sweep value in this dimension,
        # fiducial in all other dimensions (isolates convergence in one dimension)
        ref_res = fiducial.copy()
        ref_res[dim_idx] = max(values) * 1.2
        try:
            ref_model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(ref_res), spreading)
            ref_lcs = {}
            for fb in bands:
                ref_flux = ref_model.flux_density(t_lc, nu_arrays[fb])
                ref_lcs[fb] = self._extract_flux(ref_flux, has_rvs, has_ssc).copy()
        except Exception as e:
            print(f"      {_bold_red('Error')} computing reference: {e}")
            return _nan_result()

        for val in values:
            res = fiducial.copy()
            res[dim_idx] = val
            try:
                model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, tuple(res), spreading)

                band_times = []
                for fb in bands:
                    nu_arr = nu_arrays[fb]
                    va.Model.profile_reset()
                    flux = model.flux_density(t_lc, nu_arr)
                    elapsed = va.Model.profile_data().get("total", 0.0)
                    band_times.append(elapsed)

                    flux_arr = self._extract_flux(flux, has_rvs, has_ssc)
                    times_by_band[fb].append(elapsed)
                    flux_by_band[fb].append(flux_arr.tolist())

                    ref = ref_lcs[fb]
                    valid = (ref > 0) & np.isfinite(ref) & (flux_arr > 0) & np.isfinite(flux_arr)
                    if np.any(valid):
                        # Exclude exponential-decay regions: compute log-log
                        # slope of reference curve and mask |slope| > 4
                        log_t = np.log10(t_lc)
                        log_ref = np.full_like(log_t, np.nan)
                        pos = ref > 0
                        log_ref[pos] = np.log10(ref[pos])
                        slope = np.gradient(log_ref, log_t)
                        valid = valid & (np.abs(slope) <= 4)
                        # Also exclude tail > 6 orders below peak
                        if np.any(valid):
                            ref_peak = np.max(ref[valid])
                            valid = valid & (ref > ref_peak * 1e-6)
                        if np.any(valid):
                            err = np.abs(flux_arr[valid] - ref[valid]) / ref[valid]
                        else:
                            err = np.array([np.nan])
                        errors_by_band[fb].append(float(np.max(err)))
                        mean_errors_by_band[fb].append(float(np.mean(err)))
                    else:
                        errors_by_band[fb].append(np.nan)
                        mean_errors_by_band[fb].append(np.nan)
                times_ms.append(np.mean(band_times))
            except Exception as e:
                print(f"      {_bold_red('Error')} at {dimension}={val}: {e}")
                times_ms.append(np.nan)
                for b in bands:
                    times_by_band[b].append(np.nan)
                    errors_by_band[b].append(np.nan)
                    mean_errors_by_band[b].append(np.nan)
                    flux_by_band[b].append([])

        ref_flux_by_band = {b: arr.tolist() for b, arr in ref_lcs.items()}

        return DimensionConvergence(
            dimension=dimension, values=values, times_ms=times_ms,
            errors_by_band=errors_by_band, mean_errors_by_band=mean_errors_by_band, times_by_band=times_by_band,
            t_array=t_lc.tolist(), flux_by_band=flux_by_band, ref_flux_by_band=ref_flux_by_band,
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
        jets, meds = get_all_jet_names(), ["ISM", "wind"]
        all_rads = get_fwd_only_radiation_names() + get_rvs_radiation_names()
        configs = [(j, m, r, get_jet_theta_c(j) * ratio)
                   for j, m, r, ratio in product(jets, meds, all_rads, VIEWING_ANGLE_RATIOS)]
        return self._run_suite("FULL BENCHMARK SUITE", configs)

    def run_convergence_only(self, jet_types=None, medium_types=None):
        jets = jet_types or ["tophat", "gaussian", "powerlaw"]
        meds = medium_types or ["ISM", "wind"]
        return self._run_suite("CONVERGENCE TEST SUITE",
            [(j, m, "synchrotron", 0.0) for j, m in product(jets, meds)])

    def run_parallel(self, jet_types=None, medium_types=None, radiation_types=None,
                     viewing_ratios=None, n_workers=None, timeout=600):
        print(_header("PARALLEL BENCHMARK SUITE"))
        jets = jet_types or get_all_jet_names()
        meds = medium_types or ["ISM", "wind"]
        rads = radiation_types or ["synchrotron"]
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
        conv_keys = ["phi_convergence", "theta_convergence", "t_convergence"]
        for rd in results_dicts:
            timing = ComponentTiming(**rd.pop("timing"))
            convs = {}
            for k in conv_keys:
                raw = rd.pop(k, None)
                convs[k] = DimensionConvergence(**raw) if raw else None
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
            compiler=compiler, compile_flags=compile_flags, cpu=get_cpu_info(), configs=self.results,
        )
        session_dict = asdict(session)
        path = Path(filepath)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(session_dict))
        print(f"\nResults saved to {_dim(filepath)}")

    @staticmethod
    def _mean_convergence_error(conv):
        """Compute mean error across all bands for a convergence result."""
        vals = [v for band in conv.mean_errors_by_band.values() for v in band if not np.isnan(v)]
        return np.mean(vals) if vals else 0.0

    @staticmethod
    def _format_error(value):
        """Color-code an error value: green if <5%, red if >15%, plain otherwise."""
        color = _green if value < 0.05 else (_red if value > 0.15 else None)
        return color(f"{value:.2e}") if color else f"{value:.2e}"

    def print_summary(self):
        if not self.results:
            print("No results to summarize")
            return
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
                label = _bold(f"{r.jet_type}/{r.medium}")
                phi_e = self._mean_convergence_error(r.phi_convergence)
                theta_e = self._mean_convergence_error(r.theta_convergence)
                t_e = self._mean_convergence_error(r.t_convergence)
                fmt = self._format_error
                print(f"{label}: phi={fmt(phi_e)}, theta={fmt(theta_e)}, t={fmt(t_e)}")


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Benchmark Suite")
    parser.add_argument("--full", action="store_true", help="Full benchmark (all combinations)")
    parser.add_argument("--convergence", action="store_true", help="Convergence tests only")
    parser.add_argument("-j", "--parallel", type=int, default=0, metavar="N", help="Parallel workers")
    parser.add_argument("--jet", type=str, nargs="+", help="Jet type(s)")
    parser.add_argument("--medium", type=str, nargs="+", help="Medium type(s)")
    parser.add_argument("--iterations", type=int, default=5, help="Timing iterations for overview stage breakdown")
    parser.add_argument("--timeout", type=int, default=600, help="Per-config timeout in seconds (default: 600)")
    parser.add_argument("--output", type=str, default="results/benchmark_history.json", help="Output file")
    args = parser.parse_args()

    print(_header("VegasAfterglow Benchmark Suite"))
    print(f"Commit: {_bold(get_git_commit())}\nPlatform: {get_platform_info()}\nPython: {sys.version.split()[0]}")

    b = ComprehensiveBenchmark(iterations=args.iterations)
    if args.parallel > 0:
        print(f"Parallel workers: {_bold(str(args.parallel))}")
        ratios = VIEWING_ANGLE_RATIOS if args.full else [0]
        rads = get_fwd_only_radiation_names() + get_rvs_radiation_names() if args.full else ["synchrotron"]
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
