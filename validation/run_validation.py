#!/usr/bin/env python3
"""Unified validation runner for CI/CD pipelines.

This script orchestrates benchmark and regression tests, checks pass/fail criteria,
and generates validation reports. Returns exit code 0 on success, 1 on failure.

Usage:
    python -m validation.run_validation --all           # Run everything
    python -m validation.run_validation --benchmark     # Benchmarks only
    python -m validation.run_validation --regression    # Regression only
    python -m validation.run_validation --check-only    # Check existing results
"""

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Tuple

# Default to number of CPUs for parallel execution
DEFAULT_WORKERS = os.cpu_count() or 4

# Thresholds for convergence tests (from benchmark_guide.md)
BENCHMARK_MEAN_ERROR_THRESHOLD = 0.05  # 5%
BENCHMARK_MAX_ERROR_THRESHOLD = 0.10   # 10%

# Thresholds for regression tests (from regression_guide.md)
REGRESSION_SHOCK_TOLERANCE = 0.1
REGRESSION_FREQ_TOLERANCE = 0.1
REGRESSION_SPECTRAL_TOLERANCE = 0.15


def get_validation_dir() -> Path:
    """Get the validation directory path."""
    return Path(__file__).parent


def check_benchmark_results(json_path: Path) -> Tuple[bool, str]:
    """Check if benchmark results pass convergence criteria.

    Criteria:
        PASS: mean_error < 5% AND max_error < 10%
        ACCEPTABLE: mean_error < 5% AND max_error >= 10%
        FAIL: mean_error >= 5%

    Returns:
        Tuple of (passed: bool, message: str)
    """
    if not json_path.exists():
        return False, f"Benchmark results not found: {json_path}"

    try:
        with open(json_path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, f"Invalid JSON in benchmark results: {e}"

    sessions = data.get("sessions", [])
    if not sessions:
        return False, "No benchmark sessions found in results"

    latest = sessions[-1]
    configs = latest.get("configs", latest.get("results", []))

    if not configs:
        return False, "No benchmark configurations found"

    failed_configs = []
    total_configs = 0
    passed_configs = 0
    acceptable_configs = 0

    for config in configs:
        # Check convergence data
        for dim in ["phi_convergence", "theta_convergence", "t_convergence"]:
            conv_data = config.get(dim)
            if not conv_data:
                continue

            total_configs += 1
            mean_errors = conv_data.get("mean_errors", [])
            max_errors = conv_data.get("relative_errors", conv_data.get("max_errors", []))

            if not mean_errors or not max_errors:
                continue

            # Get error at fiducial resolution (last element typically)
            mean_err = mean_errors[-1] if mean_errors else 1.0
            max_err = max_errors[-1] if max_errors else 1.0

            # Handle NaN values
            if mean_err != mean_err:  # NaN check
                mean_err = 1.0
            if max_err != max_err:
                max_err = 1.0

            config_name = f"{config.get('jet_type', 'unknown')}/{config.get('medium', 'unknown')}"

            if mean_err >= BENCHMARK_MEAN_ERROR_THRESHOLD:
                failed_configs.append(f"{config_name} ({dim}): mean_error={mean_err:.1%} >= {BENCHMARK_MEAN_ERROR_THRESHOLD:.0%}")
            elif max_err >= BENCHMARK_MAX_ERROR_THRESHOLD:
                acceptable_configs += 1
            else:
                passed_configs += 1

    if failed_configs:
        msg = f"Benchmark FAILED: {len(failed_configs)} configurations failed\n"
        for fail in failed_configs[:5]:  # Show first 5 failures
            msg += f"  - {fail}\n"
        if len(failed_configs) > 5:
            msg += f"  ... and {len(failed_configs) - 5} more\n"
        return False, msg

    msg = f"Benchmark PASSED: {passed_configs} pass, {acceptable_configs} acceptable out of {total_configs} total"
    return True, msg


def check_regression_results(json_path: Path) -> Tuple[bool, str]:
    """Check if regression results pass theoretical predictions.

    Returns:
        Tuple of (passed: bool, message: str)
    """
    if not json_path.exists():
        return False, f"Regression results not found: {json_path}"

    try:
        with open(json_path) as f:
            data = json.load(f)
    except json.JSONDecodeError as e:
        return False, f"Invalid JSON in regression results: {e}"

    tests = data.get("tests", [])
    if not tests:
        return False, "No regression tests found in results"

    n_total = len(tests)
    n_passed = sum(1 for t in tests if t.get("passed", False))
    n_failed = n_total - n_passed

    if n_failed > 0:
        failed_tests = [t for t in tests if not t.get("passed", False)]
        msg = f"Regression FAILED: {n_failed}/{n_total} tests failed\n"
        for t in failed_tests[:5]:
            msg += f"  - {t.get('name', 'unknown')}: {t.get('message', 'no details')}\n"
        if len(failed_tests) > 5:
            msg += f"  ... and {len(failed_tests) - 5} more\n"
        return False, msg

    msg = f"Regression PASSED: {n_passed}/{n_total} tests passed"
    return True, msg


def run_benchmark(parallel: int = 0, quick: bool = False) -> Tuple[bool, Path]:
    """Run benchmark suite.

    Returns:
        Tuple of (success: bool, results_path: Path)
    """
    validation_dir = get_validation_dir()
    script = validation_dir / "benchmark" / "benchmark_suite.py"
    results_path = validation_dir / "benchmark" / "results" / "benchmark_history.json"

    cmd = [sys.executable, str(script)]
    if quick:
        cmd.append("--quick")
    else:
        cmd.append("--full")  # Always run full benchmark suite
        if parallel > 0:
            cmd.extend(["-j", str(parallel)])

    print(f"\n{'='*60}")
    print("Running Benchmark Suite")
    print(f"{'='*60}\n")

    try:
        result = subprocess.run(cmd, cwd=validation_dir)
        return result.returncode == 0, results_path
    except Exception as e:
        print(f"Error running benchmark: {e}")
        return False, results_path


def run_regression(quick: bool = False) -> Tuple[bool, Path]:
    """Run regression tests.

    Returns:
        Tuple of (success: bool, results_path: Path)
    """
    validation_dir = get_validation_dir()
    script = validation_dir / "regression" / "run_regression.py"
    results_path = validation_dir / "regression" / "results" / "regression_results.json"

    cmd = [sys.executable, str(script)]
    if quick:
        cmd.append("--quick")

    print(f"\n{'='*60}")
    print("Running Regression Tests")
    print(f"{'='*60}\n")

    try:
        result = subprocess.run(cmd, cwd=validation_dir)
        return result.returncode == 0, results_path
    except Exception as e:
        print(f"Error running regression: {e}")
        return False, results_path


def generate_report(benchmark_path: Path = None, regression_path: Path = None,
                    output_dir: str = "output", parallel: int = 0) -> bool:
    """Generate PDF validation report.

    Returns:
        success: bool
    """
    print(f"\n{'='*60}")
    print("Generating Validation Report")
    print(f"{'='*60}\n")

    try:
        # Import dashboard from relative path
        validation_dir = get_validation_dir()
        if str(validation_dir.parent) not in sys.path:
            sys.path.insert(0, str(validation_dir.parent))
        from validation.visualization.dashboard import ComprehensiveDashboard

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        dashboard = ComprehensiveDashboard(output_path)

        # Load data files
        benchmark_file = str(benchmark_path) if benchmark_path and benchmark_path.exists() else None
        regression_file = str(regression_path) if regression_path and regression_path.exists() else None

        # Generate full report with both data sources
        dashboard.generate_full_report(
            benchmark_file=benchmark_file,
            regression_file=regression_file,
            n_workers=parallel if parallel > 0 else DEFAULT_WORKERS
        )
        return True
    except Exception as e:
        print(f"Error generating report: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    parser = argparse.ArgumentParser(
        description="VegasAfterglow Validation Runner",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python -m validation.run_validation --all              # Run all tests (parallel by default)
    python -m validation.run_validation --benchmark        # Benchmark only
    python -m validation.run_validation --regression       # Regression only
    python -m validation.run_validation --check-only       # Check existing results
    python -m validation.run_validation --all --quick      # Quick validation
    python -m validation.run_validation --all -j 2         # Use 2 workers
    python -m validation.run_validation --all --no-report  # Skip report generation
        """
    )

    parser.add_argument("--all", action="store_true", help="Run all validation tests")
    parser.add_argument("--benchmark", action="store_true", help="Run benchmark tests")
    parser.add_argument("--regression", action="store_true", help="Run regression tests")
    parser.add_argument("--check-only", action="store_true", help="Only check existing results (no test execution)")
    parser.add_argument("--no-report", action="store_true", help="Skip PDF report generation")
    parser.add_argument("--quick", action="store_true", help="Run quick/minimal test suite")
    parser.add_argument("-j", "--parallel", type=int, default=DEFAULT_WORKERS,
                        help=f"Number of parallel workers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--output", type=str, default="output", help="Output directory for reports")
    parser.add_argument("--strict", action="store_true", help="Fail on ACCEPTABLE status (require PASS)")

    args = parser.parse_args()

    # Default to --all if no specific test selected
    if not any([args.all, args.benchmark, args.regression, args.check_only]):
        args.all = True

    validation_dir = get_validation_dir()
    benchmark_path = validation_dir / "benchmark" / "results" / "benchmark_history.json"
    regression_path = validation_dir / "regression" / "results" / "regression_results.json"

    all_passed = True
    messages = []

    print("="*60)
    print("VegasAfterglow Validation Runner")
    print("="*60)

    # Run tests (unless check-only mode)
    if not args.check_only:
        if args.all or args.benchmark:
            success, benchmark_path = run_benchmark(parallel=args.parallel, quick=args.quick)
            if not success:
                messages.append("Benchmark execution failed")
                all_passed = False

        if args.all or args.regression:
            success, regression_path = run_regression(quick=args.quick)
            if not success:
                messages.append("Regression execution failed")
                all_passed = False

    # Check results
    print(f"\n{'='*60}")
    print("Checking Results")
    print(f"{'='*60}\n")

    if args.all or args.benchmark or args.check_only:
        passed, msg = check_benchmark_results(benchmark_path)
        print(msg)
        if not passed:
            all_passed = False
            messages.append("Benchmark validation failed")

    if args.all or args.regression or args.check_only:
        passed, msg = check_regression_results(regression_path)
        print(msg)
        if not passed:
            all_passed = False
            messages.append("Regression validation failed")

    # Always generate report (even if tests fail) unless --no-report
    if not args.no_report:
        report_success = generate_report(
            benchmark_path=benchmark_path if (args.all or args.benchmark) else None,
            regression_path=regression_path if (args.all or args.regression) else None,
            output_dir=args.output,
            parallel=args.parallel
        )
        if not report_success:
            messages.append("Report generation failed")
            # Don't fail overall validation for report issues

    # Final summary
    print(f"\n{'='*60}")
    if all_passed:
        print("VALIDATION PASSED")
    else:
        print("VALIDATION FAILED")
        for msg in messages:
            print(f"  - {msg}")
    print("="*60)

    sys.exit(0 if all_passed else 1)


if __name__ == "__main__":
    main()
