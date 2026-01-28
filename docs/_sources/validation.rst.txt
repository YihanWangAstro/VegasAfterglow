Validation & Testing
====================

VegasAfterglow includes a comprehensive validation framework to ensure numerical accuracy and physical correctness. The validation suite is located in the ``validation/`` directory.

Running Validation
------------------

**Prerequisite:** The VegasAfterglow package must be installed before running validation tests. See :doc:`installation` for instructions.

The unified validation runner orchestrates all tests and generates reports:

.. code-block:: bash

   # Run full validation suite (benchmark + regression + PDF report)
   python validation/run_validation.py --all

   # Run only benchmark tests
   python validation/run_validation.py --benchmark

   # Run only regression tests
   python validation/run_validation.py --regression

   # Check existing results without re-running tests
   python validation/run_validation.py --check-only

   # Quick validation (reduced test matrix)
   python validation/run_validation.py --all --quick

   # Control parallelism (default: all CPU cores)
   python validation/run_validation.py --all -j 4

The validation runner returns exit code 0 on success and 1 on failure, making it suitable for CI/CD pipelines.

Benchmark Tests
---------------

Benchmark tests measure computational performance and verify numerical convergence. Located in ``validation/benchmark/``.

**What is tested:**

- **Performance Timing**: Execution time for various jet/medium/radiation configurations
- **Resolution Convergence**: Convergence analysis in three dimensions:

  - **Phi (azimuthal)**: Tests values from 0.1 to 0.5 radians
  - **Theta (polar)**: Tests values from 0.2 to 1.0 radians
  - **Time**: Tests values from 5 to 25 points per decade

- **Configuration Matrix**: All combinations of jet types (tophat, gaussian, powerlaw), media (ISM, wind), and viewing angles

**Pass Criteria:**

+------------+------------------+-------------------+
| Status     | Mean Error       | Max Error         |
+============+==================+===================+
| PASS       | < 5%             | < 10%             |
+------------+------------------+-------------------+
| ACCEPTABLE | < 5%             | >= 10%            |
+------------+------------------+-------------------+
| FAIL       | >= 5%            | any               |
+------------+------------------+-------------------+

Regression Tests
----------------

Regression tests verify that simulation outputs match theoretical predictions from GRB afterglow theory. Located in ``validation/regression/``.

**Shock Dynamics Tests:**

Validates power-law scaling relations Q ~ t^alpha for:

- Lorentz factor (u)
- Radius (r)
- Magnetic field (B)
- Particle number (N_p)

Across three evolutionary phases:

- **Coasting**: Free expansion (Gamma ~ const)
- **Blandford-McKee**: Self-similar relativistic deceleration
- **Sedov-Taylor**: Non-relativistic regime

**Characteristic Frequency Tests:**

Verifies time evolution of:

- Injection frequency (nu_m)
- Cooling frequency (nu_c)
- Self-absorption frequency (nu_a)

**Spectral Shape Tests:**

Checks power-law spectral indices beta = d(log F)/d(log nu) across five regimes:

- **Regime I**: nu_a < nu_m < nu_c (slow cooling)
- **Regime II**: nu_m < nu_a < nu_c
- **Regime III**: nu_a < nu_c < nu_m (fast cooling)
- **Regime IV**: nu_c < nu_a < nu_m
- **Regime V**: nu_c < nu_m < nu_a

**Tolerances:**

- Shock dynamics: 0.1
- Characteristic frequencies: 0.1
- Spectral shapes: 0.15

Validation Reports
------------------

The validation runner generates a comprehensive PDF report (``output/comprehensive_report.pdf``) containing:

1. **Title Page**: Version info, commit hash, platform details
2. **Table of Contents**: Clickable navigation with page numbers
3. **Benchmark Section**:

   - Guide pages explaining how to interpret results
   - Overview plots (on-axis and off-axis configurations)
   - Convergence summary grid with color-coded pass/fail status
   - Detailed convergence plots for each configuration

4. **Regression Section**:

   - Guide pages with theoretical background
   - Summary grid comparing measured vs expected power-law exponents
   - Detailed diagnostic plots for shock dynamics and frequencies

Directory Structure
-------------------

.. code-block:: text

   validation/
   ├── run_validation.py      # Unified CLI runner
   ├── benchmark/
   │   ├── benchmark_suite.py # Benchmark test implementation
   │   ├── configs.py         # Test configurations
   │   └── results/           # JSON output files
   ├── regression/
   │   ├── run_regression.py  # Regression test runner
   │   └── results/           # JSON output files
   ├── visualization/
   │   ├── dashboard.py       # PDF report generator
   │   ├── common.py          # Shared utilities
   │   ├── benchmark_plots.py # Benchmark visualizations
   │   ├── regression_plots.py# Regression visualizations
   │   └── guides/            # Markdown guide documents
   └── output/                # Generated PDF reports

CI/CD Integration
-----------------

The validation framework is designed for CI/CD pipelines:

.. code-block:: yaml

   # Example GitHub Actions workflow
   - name: Run validation
     run: python validation/run_validation.py --all

   - name: Upload report
     uses: actions/upload-artifact@v3
     with:
       name: validation-report
       path: output/comprehensive_report.pdf

The validation runner:

- Returns exit code 0 only if all tests pass
- Generates reports even when tests fail (for debugging)
- Supports parallel execution for faster CI runs
- Produces machine-readable JSON results
