#!/usr/bin/env python3
"""Unified validation runner for VegasAfterglow CI/CD pipelines."""

import argparse
import json
import math
import os
import subprocess
import sys
from pathlib import Path

DEFAULT_WORKERS = os.cpu_count() or 4

sys.path.insert(0, str(Path(__file__).parent.parent))
from validation.colors import (_bold, _green, _red, _bold_green, _bold_red,
                               _header, _bar)
BENCH_MEAN_THRESH = 0.05   # 5%
BENCH_MAX_THRESH = 0.10    # 10%

VALIDATION_DIR = Path(__file__).parent


def _load_json(path):
    if not path.exists():
        return None, f"Results not found: {path}"
    try:
        return json.loads(path.read_text()), None
    except json.JSONDecodeError as e:
        return None, f"Invalid JSON: {e}"


def check_benchmark_results(json_path):
    data, err = _load_json(json_path)
    if err:
        return False, err
    configs = data.get("configs", [])
    if not configs:
        return False, "No benchmark configurations found"

    failed, passed, acceptable, total = [], 0, 0, 0
    for config in configs:
        for dim in ("phi_convergence", "theta_convergence", "t_convergence"):
            conv = config.get(dim)
            if not conv:
                continue
            total += 1
            # Collect the highest-resolution error across all bands
            mean_by_band = conv.get("mean_errors_by_band", {})
            errs_by_band = conv.get("errors_by_band", {})
            if not mean_by_band and not errs_by_band:
                continue
            # Worst-case mean error across bands (last element = highest resolution)
            mean_err = 0.0
            for band_errs in mean_by_band.values():
                if band_errs:
                    v = band_errs[-1]
                    if not math.isnan(v):
                        mean_err = max(mean_err, v)
            # Worst-case max error across bands
            max_err = 0.0
            for band_errs in errs_by_band.values():
                if band_errs:
                    v = band_errs[-1]
                    if not math.isnan(v):
                        max_err = max(max_err, v)
            name = f"{config.get('jet_type', '?')}/{config.get('medium', '?')}"
            if mean_err >= BENCH_MEAN_THRESH:
                failed.append(f"{name} ({dim}): mean_error={mean_err:.1%} >= {BENCH_MEAN_THRESH:.0%}")
            elif max_err >= BENCH_MAX_THRESH:
                acceptable += 1
            else:
                passed += 1

    if failed:
        msg = f"{_bold_red('Benchmark FAILED')}: {len(failed)} configurations failed\n"
        for f in failed[:5]:
            msg += f"  - {_red(f)}\n"
        if len(failed) > 5:
            msg += f"  ... and {len(failed) - 5} more\n"
        return False, msg
    return True, f"{_bold_green('Benchmark PASSED')}: {_green(f'{passed} pass')}, {acceptable} acceptable out of {total} total"


def check_regression_results(json_path):
    data, err = _load_json(json_path)
    if err:
        return False, err
    tests = data.get("tests", [])
    if not tests:
        return False, "No regression tests found"
    n_pass = sum(1 for t in tests if t.get("passed", False))
    n_fail = len(tests) - n_pass
    if n_fail > 0:
        failed = [t for t in tests if not t.get("passed", False)]
        msg = f"{_bold_red('Regression FAILED')}: {n_fail}/{len(tests)} tests failed\n"
        for t in failed[:5]:
            detail = t.get('message', '')
            if not detail and t.get('measured') is not None:
                detail = f"measured={t['measured']:.4f}, expected={t.get('expected', '?')}, tol={t.get('tolerance', '?')}"
            msg += f"  - {_red(t.get('name', '?'))}: {detail or 'no details'}\n"
        if len(failed) > 5:
            msg += f"  ... and {len(failed) - 5} more\n"
        return False, msg
    return True, f"{_bold_green('Regression PASSED')}: {_green(f'{n_pass}/{len(tests)} tests passed')}"


def _run_suite(name, script, extra_args=()):
    print(_header(f"Running {name}"))
    try:
        result = subprocess.run([sys.executable, str(script), *extra_args], cwd=VALIDATION_DIR)
        return result.returncode == 0
    except Exception as e:
        print(f"{_bold_red('Error')} running {name}: {e}")
        return False


def run_benchmark(parallel=0):
    args = ["--full"] + (["-j", str(parallel)] if parallel > 0 else [])
    path = VALIDATION_DIR / "benchmark" / "results" / "benchmark_history.json"
    return _run_suite("Benchmark Suite", VALIDATION_DIR / "benchmark" / "benchmark_suite.py", args), path


def run_regression():
    path = VALIDATION_DIR / "regression" / "results" / "regression_results.json"
    return _run_suite("Regression Tests", VALIDATION_DIR / "regression" / "run_regression.py"), path


def generate_report(benchmark_path=None, regression_path=None, output_dir=None, parallel=0):
    print(_header("Generating Validation Report"))
    try:
        from validation.visualization.dashboard import ComprehensiveDashboard

        out = Path(output_dir or str(VALIDATION_DIR))
        out.mkdir(parents=True, exist_ok=True)
        dashboard = ComprehensiveDashboard(out)
        dashboard.generate_full_report(
            benchmark_file=str(benchmark_path) if benchmark_path and benchmark_path.exists() else None,
            regression_file=str(regression_path) if regression_path and regression_path.exists() else None,
            n_workers=parallel if parallel > 0 else DEFAULT_WORKERS)
        return True
    except Exception as e:
        print(f"{_bold_red('Error')} generating report: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(description="VegasAfterglow Validation Runner")
    parser.add_argument("--all", action="store_true", help="Run all validation tests")
    parser.add_argument("--benchmark", action="store_true", help="Run benchmark tests")
    parser.add_argument("--regression", action="store_true", help="Run regression tests")
    parser.add_argument("--check-only", action="store_true", help="Only check existing results")
    parser.add_argument("--no-report", action="store_true", help="Skip PDF report generation")
    parser.add_argument("-j", "--parallel", type=int, default=DEFAULT_WORKERS,
                        help=f"Number of parallel workers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--output", type=str, default=None, help="Output directory for reports")
    parser.add_argument("--strict", action="store_true", help="Fail on ACCEPTABLE status")
    args = parser.parse_args()

    if not any([args.all, args.benchmark, args.regression, args.check_only]):
        args.all = True

    bench_path = VALIDATION_DIR / "benchmark" / "results" / "benchmark_history.json"
    regr_path = VALIDATION_DIR / "regression" / "results" / "regression_results.json"
    all_passed, messages = True, []

    print(_header("VegasAfterglow Validation Runner"))

    # Run tests
    if not args.check_only:
        if args.all or args.benchmark:
            ok, bench_path = run_benchmark(parallel=args.parallel)
            if not ok:
                messages.append("Benchmark execution failed")
                all_passed = False
        if args.all or args.regression:
            ok, regr_path = run_regression()
            if not ok:
                messages.append("Regression execution failed")
                all_passed = False

    # Check results
    print(_header("Checking Results"))
    for should_run, checker, path, label in [
        (args.all or args.benchmark or args.check_only, check_benchmark_results, bench_path, "Benchmark"),
        (args.all or args.regression or args.check_only, check_regression_results, regr_path, "Regression"),
    ]:
        if should_run:
            passed, msg = checker(path)
            print(msg)
            if not passed:
                all_passed = False
                messages.append(f"{label} validation failed")

    # Generate report
    if not args.no_report:
        ok = generate_report(
            benchmark_path=bench_path if (args.all or args.benchmark) else None,
            regression_path=regr_path if (args.all or args.regression) else None,
            output_dir=args.output, parallel=args.parallel)
        if not ok:
            messages.append("Report generation failed")

    # Summary
    print(f"\n{_bar()}")
    print(_bold_green("VALIDATION PASSED") if all_passed else _bold_red("VALIDATION FAILED"))
    for msg in messages:
        print(f"  - {_red(msg)}")
    print(_bar())
    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
