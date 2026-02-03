# Benchmark Report Guide

This document describes the benchmark validation framework for VegasAfterglow, which assesses numerical convergence and computational performance across physical configurations.

---

## 1. Configuration Space

The benchmark suite systematically tests a grid of physical scenarios and numerical settings to characterize code behavior across the parameter space relevant for GRB afterglow modeling.

### 1.1 Physical Parameters

| Parameter | Options | Description |
|-----------|---------|-------------|
| Jet Structure | tophat, gaussian, powerlaw, two_component | Angular energy profile |
| External Medium | ISM, Wind | Density profile: constant (ISM) or r^(-2) (Wind) |
| Radiation | synchrotron, with_ssc_cooling, fast_cooling, steep/flat_spectrum, rvs_sync_thin, rvs_sync_thick | Radiation physics |
| Viewing Angle | theta_v/theta_c = 0, 2, 4 | On-axis (0) vs off-axis (>1) |

### 1.2 Numerical Resolution

Resolution controls the discretization density in each computational dimension. The fiducial values represent the recommended default settings. Note that theta and t grids enforce a minimum total point count regardless of the ppd value, so low ppd settings may not reduce the actual grid size.

| Dimension | Symbol | Unit | Fiducial | Test Range |
|-----------|--------|------|----------|------------|
| Azimuthal angle | phi | per degree | 0.15 | 0.15 - 0.6 |
| Polar angle | theta | per degree | 0.5 | 0.5 - 2.0 |
| Observer time | t | per decade | 10 | 10 - 30 |

### 1.3 Frequency Bands

Convergence is evaluated independently at three representative frequencies spanning the electromagnetic spectrum.

| Band | Frequency (Hz) |
|------|----------------|
| Radio | 10^(9) |
| Optical | 4.84 x 10^(14) |
| X-ray | 10^(18) |

---

## 2. Convergence Criteria

Numerical accuracy is quantified by comparing results at each resolution to a high-resolution reference run. Both mean and maximum relative errors across the light curve are evaluated.

| Status | Criteria |
|--------|----------|
| PASS | mean error < 5% AND max error < 15% |
| ACCEPTABLE | mean error < 5% AND max error >= 15% |
| FAIL | mean error >= 5% |

---

## 3. Summary Grid

The summary grid provides a compact overview of convergence status for all tested configurations on a single page.

### 3.1 Layout

Each cell represents one model configuration (jet, medium, viewing angle). Cells are arranged to facilitate comparison across jet types and viewing angles.

### 3.2 Cell Contents

| Line | Content |
|------|---------|
| 1 | Model ID |
| 2 | Configuration shorthand (jet/medium/angle_ratio) |
| 3 | Maximum error at fiducial resolution |

### 3.3 Color Coding

| Color | Status |
|-------|--------|
| Green | PASS |
| Blue | ACCEPTABLE |
| Pink | FAIL |
| Gray | No data |

---

## 4. Overview Plots

The overview page provides performance profiling across configurations, helping identify computational bottlenecks and compare execution times.

### 4.1 Panel Layout

| Position | Content |
|----------|---------|
| Top-left | Light curve computation time by jet type |
| Top-right | Stage breakdown (profiling build) or resolution cost scaling |
| Bottom-left | Medium comparison (ISM vs Wind) |
| Bottom-right | Wind/ISM speed ratio |

### 4.2 Timing Metric

Each configuration is timed by computing a 30-point broadband light curve (t = 10^2 to 10^8 s) at the fiducial resolution. The reported time includes dynamics computation and flux evaluation in a single `flux_density` call.

### 4.3 Stage Breakdown

When built with profiling enabled (`pip install -e ".[test]" --config-settings=cmake.define.AFTERGLOW_PROFILE=ON`), the top-right panel shows a stacked bar chart of per-stage CPU cost for each jet/medium combination. The stages correspond to the internal C++ computation pipeline:

| Stage | Description |
|-------|-------------|
| mesh | Coordinate grid generation |
| shock_dynamics | Forward/reverse shock ODE integration |
| EAT_grid | Equal arrival time surface grid |
| syn_electrons | Synchrotron electron distribution |
| syn_photons | Synchrotron photon spectrum |
| cooling | SSC/IC cooling corrections |
| sync_flux | Synchrotron flux integration |
| ic_photons | Inverse Compton photon spectrum |
| ssc_flux | SSC flux integration |

Without profiling, the panel falls back to showing total resolution cost scaling.

---

## 5. Per-Model Convergence Pages

Each configuration receives a detailed convergence analysis page showing how accuracy and performance vary with resolution.

### 5.1 Row Contents

The page displays a 4x3 grid where each column corresponds to one resolution dimension (phi, theta, t) and each row shows a different metric.

| Row | Y-axis | Threshold |
|-----|--------|-----------|
| 1 | Maximum relative error | 15% |
| 2 | Mean relative error | 5% |
| 3 | CPU time (ms) | - |
| 4 | Flux (mJy) | - |

### 5.2 Plot Features

| Feature | Meaning |
|---------|---------|
| Star marker | Fiducial resolution |
| Dashed line | Error threshold |
| Solid curves | Resolution >= fiducial |
| Dotted curves | Resolution < fiducial |

### 5.3 Status Indicator

The page header displays convergence status with color coding matching Section 3.3.

---

## 6. Interpreting Results

### 6.1 Light Curve Convergence

The bottom row plots all tested resolutions together on the same axes. Visual spread between curves indicates resolution dependence.

| Pattern | Interpretation |
|---------|----------------|
| No spread (curves overlap) | Converged |
| Spread in dashed lines only | Not converged below fiducial, acceptable at fiducial |
| Spread in solid lines | Not converged even above fiducial resolution |

### 6.2 Error Convergence

The top two rows show how errors decrease as resolution increases. The shape of these curves indicates convergence quality.

| Pattern | Interpretation |
|---------|----------------|
| Monotonic decrease | Stable convergence |
| Plateau | Noise floor or discretization limit |
| Non-monotonic | Potential numerical instability |

### 6.3 Performance Scaling

Computational cost typically scales as:
- Linear with time resolution (t_ppd)
- Quadratic with angular resolution (phi_ppd x theta_ppd)

Wind medium simulations generally require more computation due to the radially-varying density profile. However, the medium-aware adaptive grid often produces smaller grids for wind (earlier deceleration time), which can offset this cost.
