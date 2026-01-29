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

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class ComponentTiming:
    flux_single_ms: float = 0.0
    total_ms: float = 0.0
    stage_breakdown: Dict[str, float] = field(default_factory=dict)

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
    compiler, flags = "unknown", "unknown"
    cc_path = Path(__file__).parent.parent.parent / "compile_commands.json"
    if not cc_path.exists():
        return compiler, flags
    try:
        data = json.loads(cc_path.read_text())
        if data:
            parts = data[0].get("command", "").split()
            if parts:
                compiler = Path(parts[0]).name
                flag_set = {p for p in parts[1:] if p.startswith(("-O", "-march", "-std"))
                            or p in ("-fopenmp", "-ffast-math", "-DNDEBUG")}
                if flag_set:
                    flags = " ".join(sorted(flag_set))
    except (json.JSONDecodeError, KeyError, IndexError):
        pass
    if compiler != "unknown":
        try:
            compiler = subprocess.run([compiler, "--version"], capture_output=True,
                                      text=True, check=True).stdout.split("\n")[0].strip()
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass
    return compiler, flags

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
        print(f"Worker error for {jet}/{med}: {e}")
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
        model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, resolution, spreading)
        t = np.logspace(2, 7, 30)
        nu = np.full_like(t, self.nu_ref)
        va.Model.profile_reset()
        model.flux_density(t, nu)
        stage_breakdown = va.Model.profile_data()
        ms = stage_breakdown.get("total", 0.0)
        return ComponentTiming(flux_single_ms=ms, total_ms=ms, stage_breakdown=stage_breakdown)

    def _compute_lightcurve(self, model, n=50):
        t = np.logspace(2, 7, n)
        return t.tolist(), model.flux_density(t, np.full_like(t, self.nu_ref)).total.tolist()

    def _compute_spectrum(self, model, n=50):
        nu = np.logspace(9, 20, n)
        return nu.tolist(), model.flux_density(np.full_like(nu, self.t_ref), nu).total.tolist()

    def _run_dimension_convergence(self, jet_name, medium_name, radiation_name, theta_obs,
                                   dimension, values, spreading=False):
        import VegasAfterglow as va
        fiducial = list(FIDUCIAL_RESOLUTION)
        dim_idx = DIM_INDEX[dimension]
        t_lc = np.logspace(2, 7, 150)
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
                dimension=dimension, fiducial_resolution=tuple(fiducial), values=values,
                times_ms=[np.nan]*n_vals,
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
            print(f"      Error computing reference: {e}")
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
                print(f"      Error at {dimension}={val}: {e}")
                times_ms.append(np.nan)
                for b in bands:
                    times_by_band[b].append(np.nan)
                    errors_by_band[b].append(np.nan)
                    mean_errors_by_band[b].append(np.nan)
                    flux_by_band[b].append([])

        return DimensionConvergence(
            dimension=dimension, fiducial_resolution=tuple(fiducial), values=values, times_ms=times_ms,
            errors_by_band=errors_by_band, mean_errors_by_band=mean_errors_by_band, times_by_band=times_by_band,
            t_array=t_lc.tolist(), flux_by_band=flux_by_band,
            relative_errors=np.nanmean([errors_by_band[b] for b in bands], axis=0).tolist(),
        )

    def benchmark_config(self, jet_name, medium_name, radiation_name, theta_obs,
                         spreading=False, verbose=True):
        import VegasAfterglow as va
        theta_c = get_jet_theta_c(jet_name)
        ratio = theta_obs / theta_c if theta_c > 0 else 0.0
        log = print if verbose else (lambda *_a, **_k: None)

        log(f"\n  Benchmarking: {jet_name}/{medium_name}/{radiation_name}/θ_v/θ_c={ratio:.1f}")
        timing = self._compute_timing(jet_name, medium_name, radiation_name, theta_obs, self.reference_resolution, spreading)

        model = self._create_model(jet_name, medium_name, radiation_name, theta_obs, self.reference_resolution, spreading)
        lc_t, lc_flux = self._compute_lightcurve(model)
        spec_nu, spec_flux = self._compute_spectrum(model)

        result = ConfigResult(
            jet_type=jet_name, medium=medium_name, radiation=radiation_name, theta_obs=theta_obs,
            theta_obs_ratio=ratio, spreading=spreading, timing=timing,
            lightcurve_t=lc_t, lightcurve_flux=lc_flux, spectrum_nu=spec_nu, spectrum_flux=spec_flux,
        )

        for dim, vals, attr in [("phi", PHI_VALUES, "phi_convergence"),
                                ("theta", THETA_VALUES, "theta_convergence"),
                                ("t", T_VALUES, "t_convergence")]:
            log(f"    {dim} convergence...")
            setattr(result, attr, self._run_dimension_convergence(
                jet_name, medium_name, radiation_name, theta_obs, dim, vals, spreading))

        for scale in [0.5, 0.75, 1.0, 1.5, 2.0]:
            res = tuple(r * scale for r in self.reference_resolution)
            try:
                m = self._create_model(jet_name, medium_name, radiation_name, theta_obs, res, spreading)
                va.Model.profile_reset()
                flux = m.flux_density(np.array([self.t_ref]), np.array([self.nu_ref]))
                elapsed = va.Model.profile_data().get("total", 0.0)
                result.resolution_costs.append(ResolutionPoint(res, float(flux.total[0]), elapsed))
            except Exception as e:
                print(f"      Error at scale {scale}: {e}")

        log(f"    Done: {timing.flux_single_ms:.2f} ms/LC")
        self.results.append(result)
        return result

    # -----------------------------------------------------------------------
    # Suite runners
    # -----------------------------------------------------------------------

    def _run_suite(self, title, configs):
        print(f"\n{'='*70}\n{title}\n{'='*70}")
        print(f"Running {len(configs)} configurations...")
        for i, (jet, med, rad, theta) in enumerate(configs, 1):
            print(f"\n[{i}/{len(configs)}]", end="")
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
                     viewing_ratios=None, n_workers=None):
        print(f"\n{'='*70}\nPARALLEL BENCHMARK SUITE\n{'='*70}")
        jets = jet_types or get_all_jet_names()
        meds = medium_types or ["ISM", "wind"]
        rads = radiation_types or ["synchrotron_only"]
        ratios = viewing_ratios or [0]
        n_workers = n_workers or max(1, mp.cpu_count() - 1)

        configs = [(j, m, r, get_jet_theta_c(j) * ratio, False, self.iterations, self.reference_resolution)
                   for j, m, r, ratio in product(jets, meds, rads, ratios)]

        print(f"Running {len(configs)} configurations with {n_workers} workers...")
        completed, results_dicts = 0, []
        with mp.Pool(n_workers, initializer=_init_worker) as pool:
            for rd in pool.imap_unordered(_benchmark_worker, configs):
                completed += 1
                if rd:
                    results_dicts.append(rd)
                    print(f"  [{completed}/{len(configs)}] {rd.get('jet_type','?')}/{rd.get('medium','?')}: "
                          f"single-freq LC {rd.get('timing',{}).get('flux_single_ms',0):.1f} ms")
                else:
                    print(f"  [{completed}/{len(configs)}] FAILED")

        print(f"\nCompleted {len(results_dicts)}/{len(configs)} configurations")
        for rd in results_dicts:
            timing = ComponentTiming(**rd.pop("timing"))
            convs = {k: DimensionConvergence(**rd.pop(k)) if rd.get(k) else None
                     for k in ["phi_convergence", "theta_convergence", "t_convergence"]}
            for k in convs: rd.pop(k, None)
            costs = [ResolutionPoint(**rc) for rc in rd.pop("resolution_costs", [])]
            self.results.append(ConfigResult(timing=timing, **convs, resolution_costs=costs, **rd))
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
        path = Path(filepath)
        history = json.loads(path.read_text()) if path.exists() else {"sessions": []}
        history["sessions"].append(asdict(session))
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(history, indent=2))
        print(f"\nResults saved to {filepath}")

    def print_summary(self):
        if not self.results:
            print("No results to summarize"); return
        print(f"\n{'='*70}\nBENCHMARK SUMMARY\n{'='*70}")

        by_jet = {}
        for r in self.results:
            by_jet.setdefault(r.jet_type, []).append(r)

        print(f"\n{'Jet Type':<20} {'SingleLC (ms)':<15}")
        print("-" * 35)
        for jet, results in by_jet.items():
            print(f"{jet:<20} {np.mean([r.timing.flux_single_ms for r in results]):<15.2f}")

        print("\nConvergence (max relative error):")
        print("-" * 80)
        for r in self.results:
            if r.phi_convergence:
                get_max = lambda c: max((e for e in c.relative_errors if not np.isnan(e)), default=0)
                print(f"{r.jet_type}/{r.medium}: phi={get_max(r.phi_convergence):.2e}, "
                      f"theta={get_max(r.theta_convergence):.2e}, t={get_max(r.t_convergence):.2e}")


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Benchmark Suite")
    parser.add_argument("--full", action="store_true", help="Full benchmark (all combinations)")
    parser.add_argument("--convergence", action="store_true", help="Convergence tests only")
    parser.add_argument("-j", "--parallel", type=int, default=0, metavar="N", help="Parallel workers")
    parser.add_argument("--jet", type=str, nargs="+", help="Jet type(s)")
    parser.add_argument("--medium", type=str, nargs="+", help="Medium type(s)")
    parser.add_argument("--iterations", type=int, default=1, help="Timing iterations")
    parser.add_argument("--output", type=str, default="results/benchmark_history.json", help="Output file")
    args = parser.parse_args()

    print(f"{'='*70}\nVegasAfterglow Benchmark Suite\n{'='*70}")
    print(f"Commit: {get_git_commit()}\nPlatform: {get_platform_info()}\nPython: {sys.version.split()[0]}")

    b = ComprehensiveBenchmark(iterations=args.iterations)
    if args.parallel > 0:
        print(f"Parallel workers: {args.parallel}")
        ratios = VIEWING_ANGLE_RATIOS if args.full else [0]
        rads = get_radiation_names_no_ssc() if args.full else ["synchrotron_only"]
        b.run_parallel(jet_types=args.jet, medium_types=args.medium, radiation_types=rads,
                       viewing_ratios=ratios, n_workers=args.parallel)
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
