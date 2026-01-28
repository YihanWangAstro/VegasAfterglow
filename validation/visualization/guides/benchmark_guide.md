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
| Viewing Angle | theta_v/theta_c = 0, 2, 4 | On-axis (0) vs off-axis (>1) |

### 1.2 Numerical Resolution

Resolution controls the discretization density in each computational dimension. The fiducial values represent the recommended default settings.

| Dimension | Symbol | Unit | Fiducial | Test Range |
|-----------|--------|------|----------|------------|
| Azimuthal angle | phi | per degree | 0.3 | 0.1 - 0.5 |
| Polar angle | theta | per degree | 0.3 | 0.2 - 1.0 |
| Observer time | t | per decade | 10 | 5 - 25 |

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
| PASS | mean error < 5% AND max error < 10% |
| ACCEPTABLE | mean error < 5% AND max error >= 10% |
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
| Top-right | Component timing breakdown |
| Bottom-left | Medium comparison (ISM vs Wind) |
| Bottom-right | Wind/ISM speed ratio |

### 4.2 Component Timing Categories

The timing breakdown reveals which computational stages dominate the total execution time.

| Component | Description |
|-----------|-------------|
| Initialization | Model setup and parameter validation |
| Shock/Electron | Blast wave dynamics and particle distribution |
| Flux | Radiation transfer calculation |
| Multi-frequency | Grid evaluation across frequency bands |

---

## 5. Per-Model Convergence Pages

Each configuration receives a detailed convergence analysis page showing how accuracy and performance vary with resolution.

### 5.1 Row Contents

The page displays a 4x3 grid where each column corresponds to one resolution dimension (phi, theta, t) and each row shows a different metric.

| Row | Y-axis | Threshold |
|-----|--------|-----------|
| 1 | Maximum relative error | 10% |
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

Wind medium simulations generally require more computation due to the radially-varying density profile.
