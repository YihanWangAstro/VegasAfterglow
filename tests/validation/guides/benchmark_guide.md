# Benchmark Report Guide

This document describes the benchmark validation framework for VegasAfterglow, which assesses numerical convergence and computational performance across physical configurations.

---

## 1. Configuration Space

The benchmark suite systematically tests a grid of physical scenarios and numerical settings to characterize code behavior across the parameter space relevant for GRB afterglow modeling.

### 1.1 Physical Parameters

| Parameter | Options | Description |
|-----------|---------|-------------|
| Jet Structure | tophat, gaussian, powerlaw, two_component | Angular energy profile |
| External Medium | ISM, Wind | Density profile: constant (ISM) or $\propto r^{-2}$ (Wind) |
| Radiation | synchrotron, full_ssc, ssc_kn, rvs_sync_thin, rvs_sync_thick | Radiation physics |
| Viewing Angle | $\theta_v/\theta_c$ = 0, 1, 2, 4 | On-axis (0) vs off-axis (>=1) |

Microphysics-only variants (fast cooling, steep/flat electron spectrum) are covered by the regression suite instead.

### 1.2 Numerical Resolution

Resolution controls the discretization density in each computational dimension. The fiducial values represent the recommended default settings. Note that $\theta$ and $t$ grids enforce a minimum total point count regardless of the ppd value, so low ppd settings may not reduce the actual grid size.

| Dimension | Symbol | Unit | Fiducial | Test Range |
|-----------|--------|------|----------|------------|
| Azimuthal angle | $\phi$ | per degree | 0.1 | 0.1 - 0.25 |
| Polar angle | $\theta$ | per degree | 0.25 | 0.25 - 1.0 |
| Observer time | $t$ | per decade | 10 | 5 - 20 |

### 1.3 Frequency Bands

Convergence is evaluated independently at representative frequencies spanning the electromagnetic spectrum.

| Band | Frequency (Hz) |
|------|----------------|
| Radio | 10^9 |
| Optical | 4.84 $\times$ 10^14 |
| X-ray | 10^18 |
| TeV (SSC configurations only) | 2.4 $\times$ 10^26 |

---

## 2. Convergence Criteria

Numerical accuracy is quantified by comparing results at each resolution to a high-resolution reference run. Both mean and maximum relative errors across the light curve are evaluated.

| Status | Criteria |
|--------|----------|
| PASS | mean error < 5% AND max error < 15% |
| ACCEPTABLE | mean error < 5% AND max error >= 15% |
| FAIL | mean error >= 5% |

---

## 3. Convergence Status Grid

The HTML report (`tests/report.py`) renders a compact convergence status grid, one table per radiation mode.

### 3.1 Layout

Rows are jet/medium combinations; columns are viewing angle ratios, subdivided by resolution dimension ($\phi$, $\theta$, $t$). Each cell classifies one configuration and dimension at the fiducial resolution; hovering a cell shows its mean and maximum relative errors.

### 3.2 Status Icons

| Icon | Status |
|------|--------|
| ✓ | PASS |
| △ | ACCEPTABLE |
| ✕ | FAIL |
| · | No data |

---

## 4. Timing Charts

The performance section profiles execution time across configurations, helping identify computational bottlenecks and compare execution times.

### 4.1 Chart Layout

| Chart | Content |
|-------|---------|
| Median light-curve time | Per radiation mode, split into on-axis ($\theta_v$ = 0) and off-axis ($\theta_v/\theta_c$ = 1, 2, 4) |
| Per-configuration breakdown | One chart per radiation mode: jet $\times$ medium $\times$ viewing angle, stacked by stage |

The same per-stage data drives the README performance SVGs (`tests/validation/benchmark_svg.py`).

### 4.2 Timing Metric

Each configuration is timed by computing a 30-point broadband light curve ($t$ = 10^2 to 10^8 s) at the fiducial resolution. The reported time includes dynamics computation and flux evaluation in a single `flux_density` call.

### 4.3 Stage Breakdown

When built with profiling enabled (`pip install -e ".[test]" --config-settings=cmake.define.AFTERGLOW_PROFILE=ON`), timings are broken down by stage. The stages correspond to the internal C++ computation pipeline:

| Stage | Description |
|-------|-------------|
| dynamics | Forward/reverse shock ODE integration |
| EAT_grid | Equal arrival time surface grid |
| syn_electrons | Synchrotron electron distribution |
| syn_photons | Synchrotron photon spectrum |
| cooling | SSC/IC cooling corrections |
| sync_flux | Synchrotron flux integration |
| ic_photons | Inverse Compton photon spectrum |
| ssc_flux | SSC flux integration |

Without profiling, only the total wall time is recorded.

---

## 5. Convergence Curves

Error-vs-resolution panels show how accuracy varies with resolution in each dimension ($\phi$, $\theta$, $t$).

### 5.1 Panel Contents

| Panel set | Content |
|-----------|---------|
| Aggregated | One panel per radiation mode and dimension: median mean error per frequency band across configurations |
| Detail | The individual configuration/dimension checks with the largest fiducial errors |

### 5.2 Plot Features

| Feature | Meaning |
|---------|---------|
| Colored curves | Mean relative error per frequency band |
| Gray curve | Worst maximum error across bands |
| Vertical tick | Fiducial resolution |
| Dashed lines | 5% mean and 15% max error thresholds |

---

## 6. Interpreting Results

### 6.1 Error Convergence

The error curves show how errors decrease as resolution increases. The shape of these curves indicates convergence quality.

| Pattern | Interpretation |
|---------|----------------|
| Monotonic decrease | Stable convergence |
| Plateau | Noise floor or discretization limit |
| Non-monotonic | Potential numerical instability |

### 6.2 Performance Scaling

Computational cost typically scales as:
- Linear with time resolution ($t_{\rm ppd}$)
- Quadratic with angular resolution ($\phi_{\rm ppd} \times \theta_{\rm ppd}$)

Wind medium simulations generally require more computation due to the radially-varying density profile. However, the medium-aware adaptive grid often produces smaller grids for wind (earlier deceleration time), which can offset this cost.
