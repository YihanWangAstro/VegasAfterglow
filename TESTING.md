# Testing VegasAfterglow

The test infrastructure has three tiers with one shared source of physics
truth. Every tier runs the real compiled code; they differ in depth, runtime,
and when they run.

| Tier | What | Runtime | When it runs |
| --- | --- | --- | --- |
| 1. Unit / API | C++ unit tests (`tests/cpp/`, Boost.Test) + Python API tests (`tests/python/`) | seconds | every push/PR (3 platforms) |
| 2. Physics regression | closure relations, exact invariants, golden baselines, corner sweep (`tests/python/`, pytest markers) | seconds | every push/PR (3 platforms) |
| 3. Full validation | three-phase dynamics scalings, all spectral regimes, reverse shock, resolution convergence, timing (`tests/validation/`) | minutes–hours | release / weekly schedule / manual |

## Running

One entry point runs everything with a per-tier summary:

```bash
make test                  # everything: C++ + pytest + full validation + test-report.html (slow)
make test-quick            # C++ unit tests + pytest only, still writes the report (seconds)
make test-cpp              # C++ only
make test-py               # pytest only
make test-physics          # pytest -m "physics or golden" only
make test-cov              # pytest + coverage summary + htmlcov/index.html
make test-report           # alias of make test (the report is always written)
make test-validation       # the full validation suite only (slow)
make test-demos            # the Makefile-built demo/benchmark executables
```

These wrap `tests/run_all.sh`, which also takes `--build` (rebuild the C++
tests and reinstall the Python module first) and composable flags
(`--quick --cov --no-html`). Individual invocations:

```bash
# everything fast (tiers 1+2)
python -m pytest tests/ -q

# only the quantitative physics checks
python -m pytest tests/ -m physics -q

# only golden-baseline regression
python -m pytest tests/ -m golden -q

# C++ unit tests
cmake -B build -DAFTERGLOW_TESTS=ON && cmake --build build -j && ./build/afterglow_tests

# full validation suite (slow; use --fast to skip convergence scans)
python tests/validation/run_validation.py --all            # add --strict to fail on ACCEPTABLE
```

## Tier 2: the physics test families

- **`test_physics_invariants.py`** — invariants that hold to numerical
  precision by construction: flux ∝ 1/d_L² exactly, total = Σ components,
  the exact redshift transformation, axisymmetric-flag consistency,
  series-vs-grid agreement, band flux vs integrated flux density. A violation
  here is a bug, never "physics tolerance".
- **`test_closure_relations.py`** (`-m physics`) — temporal/spectral
  power-law indices vs standard afterglow theory (SPN98, Granot & Sari 2002):
  α and β in ISM/wind, above/below the cooling break, the α(p) closure, jet
  break steepening, magnetar plateau flattening, thick-shell reverse-shock
  peak at T(1+z), SSC/synchrotron ratio scaling with ε_e/ε_B. Every tolerance
  in the file is calibrated against the code's measured output (smooth
  spectral shapes curve local slopes ~0.04–0.10 away from the asymptotic
  textbook values); the comment on each test records the measured value.
- **`test_shock_scalings.py`** (`-m physics`) — the bridge into tier 3: it
  imports the *same* scaling tables, model factory, and power-law fitter that
  `tests/validation/regression` uses for the release report, and runs the
  Blandford-McKee-phase subset on every push. Expected exponents live in ONE
  place: `tests/validation/regression/run_regression.py`.
- **`test_golden.py`** (`-m golden`) — recomputes three canonical
  configurations and compares all five flux components to committed baselines
  (`tests/python/golden/*.npz`) at rtol=2e-3 (absorbs cross-platform
  fast-math/SIMD differences, catches real drift). This is the guard against
  *silent* drift that stays inside the analytic tolerances of the other
  tests. After an **intended** physics change: inspect the diff, regenerate
  via `python tests/python/golden/regenerate.py`, and commit the new
  baselines together with the change.
- **`test_parameter_corners.py`** (`-m corners`) — every jet type × medium ×
  radiation switch × geometry corner the API exposes, asserting finite,
  positive, correctly-shaped output through every output method.

## C++ physics tests

`tests/cpp/` carries the quantitative unit-level checks: Granot & Sari
asymptotic spectral slopes (slow *and* fast cooling, several p values),
electron-distribution segment slopes, self-absorbed spectral rise,
strong-shock jump limits, Klein-Nishina cross-section values, Compton-Y
limits, synchrotron kernel asymptotics, and the analytic scale relations.
Convention: every assertion is probed against the code first; where the
smooth spectral model deviates from a sharp-asymptote textbook value, the
test pins the measured behavior with a comment saying so.

## Reports and visualization

- **Every `make test` run writes `test-report.html`** (skip with `--no-html`) —
  a single self-contained page covering everything, organized by *what kind of
  correctness* each check establishes: **code correctness** (C++ unit + Python API tests), **physics
  tests** (closure relations, exact invariants, golden baselines), **physics
  validation** (the analytic power-law scaling checks from
  `tests/validation/regression`: evolution and spectral figures with fitted
  phase windows and expected-slope guides, plus every check as a
  measured-vs-expected interval chart with tolerance band), and
  **performance** (benchmark timing, per-stage breakdown, resolution
  convergence). Verdict banner, stat tiles, outcomes-by-suite and slowest-test
  charts, searchable per-file tables with duration bars and expandable
  failures. Light/dark aware, zero external assets. The validation and
  performance sections appear whenever the corresponding results JSONs exist
  (run `python tests/validation/run_validation.py --all` to refresh them; the
  run timestamp is shown in the report).
- **Every CI run** builds the same unified report per OS
  (`test-report-<os>.html`, uploaded as an artifact on the Tests workflow).
  CI runs the C++ suite on Linux/macOS and pytest on all three platforms.
- **Full validation** publishes the unified report to GitHub Pages for all
  three platforms: `reports/latest/` carries the macOS report, with
  `linux/` and `windows/` siblings switchable from the report's nav bar.
  It is refreshed by every release, by the weekly scheduled run (Deploy
  Validation Report workflow, Mondays), and by manual dispatch; each
  release additionally gets an immutable snapshot at `reports/v<version>/`,
  listed in the `reports/` index. Per-push reports live only as CI
  artifacts (90-day expiry). Tier 2's
  `test_shock_scalings.py` shares the validation suite's scaling tables, so
  a table change is exercised on the next push, not at the next release.

## Where a new test belongs

- Testing a formula/function in isolation → `tests/cpp/` (exact or analytic
  tolerance).
- Testing an end-to-end observable against theory → `test_closure_relations.py`
  (calibrate the tolerance empirically first; record the measured value in a
  comment).
- Guarding numerics you can't derive analytically → add a golden config.
- New API surface or physics switch → `test_parameter_corners.py` + an API
  test file.

> **Per-environment rebuilds:** the editable install compiles a separate
> module per conda/venv environment — after any C++ change, rerun
> `pip install -e .` in each environment you test from.
