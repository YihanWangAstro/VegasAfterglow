Validation & Testing
====================

VegasAfterglow includes a comprehensive validation framework to ensure numerical accuracy and physical correctness. The validation suite is located in the ``tests/validation/`` directory, as the slow tier of the unified test framework (see ``TESTING.md`` in the repository root for the fast tiers: C++/Python unit tests and physics regression tests that run on every push).

Running Validation
------------------

**Prerequisite:** Install with validation support, which includes per-stage CPU profiling:

.. code-block:: bash

   pip install -e ".[test]" --config-settings=cmake.define.AFTERGLOW_PROFILE=ON

The unified validation runner orchestrates all tests and generates reports:

.. code-block:: bash

   # Run full validation suite (benchmark + regression)
   python tests/validation/run_validation.py --all

   # Run only benchmark tests
   python tests/validation/run_validation.py --benchmark

   # Run only regression tests
   python tests/validation/run_validation.py --regression

   # Check existing results without re-running tests
   python tests/validation/run_validation.py --check-only

   # Fast mode: skip resolution convergence scans (timing + regression only)
   python tests/validation/run_validation.py --all --fast

   # Control parallelism (default: all CPU cores)
   python tests/validation/run_validation.py --all -j 4

``--fast`` is useful for a quick sanity check — it runs performance timing and regression tests but skips the resolution convergence scans, so it completes significantly faster.

The validation runner returns exit code 0 on success and 1 on failure, making it suitable for CI/CD pipelines.

Benchmark Tests
---------------

Benchmark tests measure computational performance and verify numerical convergence. Located in ``tests/validation/benchmark/``.

**What is tested:**

- **Performance Timing**: Execution time for various jet/medium/radiation configurations
- **Resolution Convergence**: Convergence analysis in three dimensions:

  - **Phi (azimuthal)**: Tests values from 0.15 to 0.3 radians
  - **Theta (polar)**: Tests values from 0.5 to 1.25 radians
  - **Time**: Tests values from 10 to 25 points per decade

- **Configuration Matrix**: All combinations of jet types (tophat, gaussian, powerlaw), media (ISM, wind), and viewing angles

**Pass Criteria:**

+------------+------------------+-------------------+
| Status     | Mean Error       | Max Error         |
+============+==================+===================+
| PASS       | < 5%             | < 15%             |
+------------+------------------+-------------------+
| ACCEPTABLE | < 5%             | >= 15%            |
+------------+------------------+-------------------+
| FAIL       | >= 5%            | any               |
+------------+------------------+-------------------+

Regression Tests
----------------

Regression tests verify that simulation outputs match theoretical predictions from GRB afterglow theory. Located in ``tests/validation/regression/``.

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
- Maximum synchrotron frequency (nu_M)

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

Validation results feed the unified HTML test report (``tests/report.py``, or
``make test-report`` from the repository root), which renders:

- **Physics validation**: evolution and spectral figures with fitted phase
  windows and expected-slope guides, plus every check as a
  measured-vs-expected interval chart with tolerance band
- **Performance**: benchmark timing by radiation type, per-stage CPU
  breakdown, and resolution convergence classification

The report is a single self-contained page with no external assets. The
latest release report is published to GitHub Pages at
`reports/latest <https://yihanwangastro.github.io/VegasAfterglow/reports/latest/>`_.
Methodology guides live in ``tests/validation/guides/``.

Directory Structure
-------------------

.. code-block:: text

   tests/validation/
   ├── run_validation.py          # Unified CLI runner
   ├── common.py                  # Shared constants and helpers
   ├── benchmark_svg.py           # README performance SVG generator
   ├── benchmark/
   │   ├── benchmark_suite.py     # Benchmark test implementation
   │   ├── configs.py             # Test configurations
   │   └── results/               # JSON output files
   ├── regression/
   │   ├── run_regression.py      # Regression test runner
   │   └── results/               # JSON output files
   └── guides/                    # Markdown methodology guides

Profiling Per-Stage CPU Cost
----------------------------

When installed with validation support (see above), per-stage CPU profiling is automatically enabled. The unified report's performance section shows a stacked bar chart breaking down CPU time by internal C++ computation stage:

+------------------+------------------------------------------+
| Stage            | Description                              |
+==================+==========================================+
| mesh             | Coordinate grid generation               |
+------------------+------------------------------------------+
| shock_dynamics   | Forward/reverse shock ODE integration    |
+------------------+------------------------------------------+
| EAT_grid         | Equal arrival time surface grid          |
+------------------+------------------------------------------+
| syn_electrons    | Synchrotron electron distribution        |
+------------------+------------------------------------------+
| syn_photons      | Synchrotron photon spectrum              |
+------------------+------------------------------------------+
| cooling          | SSC/IC cooling corrections               |
+------------------+------------------------------------------+
| sync_flux        | Synchrotron flux integration             |
+------------------+------------------------------------------+
| ic_photons       | Inverse Compton photon spectrum           |
+------------------+------------------------------------------+
| ssc_flux         | SSC flux integration                     |
+------------------+------------------------------------------+

Profiling data is also accessible directly from Python when built with profiling:

.. code-block:: python

   import VegasAfterglow as ag

   model = ag.Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad)

   # Reset, run, and retrieve stage timing
   ag.Model.profile_reset()
   result = model.flux_density(t, nu)
   print(ag.Model.profile_data())
   # {'total': 6.4, 'shock_dynamics': 5.1, 'syn_photons': 0.4, ...}

To return to the normal (zero-overhead) build:

.. code-block:: bash

   pip install -e .

CI/CD Integration
-----------------

The validation framework is designed for CI/CD pipelines:

.. code-block:: yaml

   # Example GitHub Actions workflow
   - name: Run validation
     run: python tests/validation/run_validation.py --all

   - name: Build unified report
     run: |
       python tests/report.py py-junit.xml \
         --validation tests/validation/regression/results/regression_results.json \
         --benchmark tests/validation/benchmark/results/benchmark_history.json \
         -o test-report.html

The validation runner:

- Returns exit code 0 only if all tests pass
- Supports parallel execution for faster CI runs
- Produces machine-readable JSON results consumed by the unified report
