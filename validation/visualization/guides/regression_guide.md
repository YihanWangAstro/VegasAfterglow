# Regression Test Report Guide

This document describes the regression validation framework for VegasAfterglow, which verifies that simulation outputs reproduce the expected power-law scaling relations from GRB afterglow theory.

---

## 1. Test Categories

The regression suite validates physical quantities against analytical predictions derived from standard synchrotron afterglow theory. Tests cover both forward and reverse shocks:

**Forward shock:**

1. **Shock dynamics**: Lorentz factor, radius, magnetic field, particle number
2. **Characteristic frequencies**: nu_m (injection), nu_c (cooling), nu_M (maximum)
3. **Spectral shapes**: Power-law indices in different frequency regimes

**Reverse shock:**

4. **Thin-shell dynamics**: Shock dynamics during and after reverse shock crossing (short engine duration)
5. **Thick-shell dynamics**: Shock dynamics during and after reverse shock crossing (long engine duration)
6. **Thin-shell frequencies**: Characteristic frequency evolution during crossing phase
7. **Thick-shell frequencies**: Characteristic frequency evolution during crossing phase

All regression tests use a model resolution of `(0.3, 2, 15)` to ensure well-resolved evolution for accurate power-law fitting.

---

## 2. Evolutionary Phases

### 2.1 Forward Shock

The blast wave passes through distinct dynamical phases as it decelerates. Each phase exhibits characteristic power-law behavior that serves as a validation target.

| Phase | ISM Time Range | Wind Time Range | Physics |
|-------|----------------|-----------------|---------|
| Coasting | 0.1 - 1 s | 0.01 - 0.1 s | Free expansion, Gamma ~ const |
| Blandford-McKee | 500 - 5000 s | 10^(4) - 10^(5) s | Self-similar deceleration |
| Deep Newtonian | 10^(12) - 10^(13) s | 10^(14) - 10^(15) s | Non-relativistic, u < 0.1 |

### 2.2 Reverse Shock (Thin Shell)

Short engine duration (tau ~ 0.01 s). The reverse shock crosses the ejecta quickly, then evolves as a decaying blast wave.

| Phase | ISM Time Range | Wind Time Range | Physics |
|-------|----------------|-----------------|---------|
| Crossing | 0.01 - 0.1 s | 5 - 50 s | Reverse shock traversing ejecta |
| Blandford-McKee | 10^(5) - 10^(6) s | 10^(6) - 10^(7) s | Post-crossing deceleration |
| Deep Newtonian | 10^(12) - 10^(13) s | 10^(14) - 10^(15) s | Non-relativistic |

### 2.3 Reverse Shock (Thick Shell)

Long engine duration (tau ~ 10^4-10^5 s). The reverse shock crosses while the engine is still active, producing different scaling.

| Phase | ISM Time Range | Wind Time Range | Physics |
|-------|----------------|-----------------|---------|
| Crossing | 5000 - 50000 s | 100 - 1000 s | Reverse shock traversing during engine activity |
| Blandford-McKee | 10^(7) - 10^(8) s | 10^(7) - 10^(8) s | Post-crossing deceleration |
| Deep Newtonian | 10^(12) - 10^(13) s | 10^(14) - 10^(15) s | Non-relativistic |

---

## 3. Expected Scaling Relations

### 3.1 Forward Shock Dynamics

Physical quantities at the shock front scale as power laws with observer time: Q ~ t^(alpha). The exponents depend on the external medium density profile.

**ISM (constant density n)**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Coasting | 0 | 1 | 0 | 3 |
| Blandford-McKee | -3/8 | 1/4 | -3/8 | 3/4 |
| Deep Newtonian | -3/5 | 2/5 | -3/5 | 6/5 |

**Wind (density ~ r^(-2))**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Coasting | 0 | 1 | -1 | 1 |
| Blandford-McKee | -1/4 | 1/2 | -3/4 | 1/2 |
| Deep Newtonian | -1/3 | 2/3 | -1 | 2/3 |

### 3.2 Forward Characteristic Frequencies

The synchrotron spectrum is characterized by break frequencies that evolve with time: nu ~ t^(alpha).

**ISM**

| Phase | nu_m | nu_c | nu_M |
|-------|------|------|------|
| Coasting | 0 | -2 | 0 |
| Blandford-McKee | -3/2 | -1/2 | -3/8 |
| Deep Newtonian | -3/5 | -1/5 | 0 |

**Wind**

| Phase | nu_m | nu_c | nu_M |
|-------|------|------|------|
| Coasting | -1 | -1 | 0 |
| Blandford-McKee | -3/2 | 1/2 | -1/4 |
| Deep Newtonian | -1 | 1 | 0 |

### 3.3 Reverse Shock Dynamics (Thin Shell)

**ISM**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Crossing | 3/2 | 1 | 0 | 3/2 |
| Blandford-McKee | -1/4 | 1/4 | — | 0 |
| Deep Newtonian | -2/5 | 2/5 | — | 0 |

**Wind**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Crossing | 1/2 | 1 | -1 | 1/2 |
| Blandford-McKee | — | 1/2 | — | 0 |
| Deep Newtonian | — | 2/3 | — | 0 |

"—" indicates quantities not tested (insufficiently clean power-law for reliable fitting).

### 3.4 Reverse Shock Dynamics (Thick Shell)

**ISM**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Crossing | 1/4 | 1/2 | -1/4 | 1 |
| Blandford-McKee | -1/4 | 1/4 | — | 0 |
| Deep Newtonian | -2/5 | 2/5 | — | 0 |

**Wind**

| Phase | u | r | B | N_p |
|-------|---|---|---|-----|
| Crossing | 0 | 1 | -1 | 1 |
| Blandford-McKee | — | 1/2 | — | 0 |
| Deep Newtonian | — | 2/3 | — | 0 |

### 3.5 Reverse Shock Characteristic Frequencies

Frequency scaling is validated during the crossing phase. Post-crossing (BM and deep Newtonian) frequency tests are not performed for reverse shocks since the reverse shock material is no longer freshly shocked.

**Thin Shell — Crossing Phase**

| Medium | nu_m | nu_c | nu_M |
|--------|------|------|------|
| ISM | 0 | -2 | 0 |
| Wind | -1 | 1 | 0 |

**Thick Shell — Crossing Phase**

| Medium | nu_m | nu_c | nu_M |
|--------|------|------|------|
| ISM | — | -1 | -1/4 |
| Wind | -1 | 1 | 0 |

---

## 4. Spectral Regimes

The synchrotron spectrum consists of power-law segments joined at break frequencies. The spectral index beta = d(log F)/d(log nu) depends on the ordering of nu_a, nu_m, and nu_c.

### Regime I: nu_a < nu_m < nu_c (Slow Cooling)

The standard slow-cooling spectrum where electrons cool on timescales longer than the dynamical time.

| Frequency Range | beta |
|-----------------|------|
| nu < nu_a | 2 |
| nu_a < nu < nu_m | 1/3 |
| nu_m < nu < nu_c | -(p-1)/2 |
| nu > nu_c | -p/2 |

### Regime II: nu_m < nu_a < nu_c

Self-absorption frequency lies between the injection and cooling frequencies.

| Frequency Range | beta |
|-----------------|------|
| nu < nu_m | 2 |
| nu_m < nu < nu_a | 5/2 |
| nu_a < nu < nu_c | -(p-1)/2 |
| nu > nu_c | -p/2 |

### Regime III: nu_a < nu_c < nu_m (Fast Cooling)

Fast-cooling regime where electrons radiate most of their energy before the next dynamical time.

| Frequency Range | beta |
|-----------------|------|
| nu < nu_a | 2 |
| nu_a < nu < nu_c | 1/3 |
| nu_c < nu < nu_m | -1/2 |
| nu > nu_m | -p/2 |

### Regime IV: nu_c < nu_a < nu_m

Self-absorption occurs above the cooling frequency in fast-cooling conditions.

| Frequency Range | beta |
|-----------------|------|
| nu < nu_c | 2 |
| nu_c < nu < nu_a | 2 |
| nu_a < nu < nu_m | -1/2 |
| nu > nu_m | -p/2 |

### Regime V: nu_c < nu_m < nu_a

Heavily self-absorbed fast-cooling spectrum.

| Frequency Range | beta |
|-----------------|------|
| nu < nu_c | 2 |
| nu_c < nu < nu_m | 2 |
| nu_m < nu < nu_a | 5/2 |
| nu > nu_a | -p/2 |

The electron index p = 2.2 is used, giving -(p-1)/2 = -0.6 and -p/2 = -1.1.

---

## 5. Summary Grid Interpretation

The summary page provides a compact view of all regression tests, allowing rapid identification of any deviations from expected behavior.

**Grid Layout**

- Rows: Physical quantities (u, r, B, N_p for shock; nu_m, nu_c for frequencies; regimes I-V for spectra)
- Columns: Medium (ISM/Wind) subdivided by phase (Coasting, BM, Deep Newtonian)

**Cell Contents**

Each cell displays the comparison between simulation and theory:

- Top value: Measured power-law exponent
- Bottom value: Expected theoretical value

**Color Coding**

| Color | Status | Criterion |
|-------|--------|-----------|
| Green | Pass | |measured - expected| < tolerance |
| Red | Fail | |measured - expected| >= tolerance |
| Gray | N/A | Insufficient data or not applicable |

**Tolerances**

- Shock dynamics: 0.1
- Characteristic frequencies: 0.1
- Spectral shapes: 0.15

---

## 6. Detailed Plot Interpretation

The detailed plots provide diagnostic information for understanding any discrepancies identified in the summary grid.

### Forward Shock Dynamics Plots (2x2 grid per medium)

Each panel shows two subplots:

- Upper: Quantity vs time (log-log), with phase regions color-coded
- Lower: Local power-law exponent d(log Q)/d(log t)

Dashed lines indicate expected scaling. Markers show fitted values. Agreement between markers and dashed lines confirms correct implementation.

### Reverse Shock Dynamics Plots

Same layout as forward shock, with separate pages for thin-shell and thick-shell configurations. The crossing phase is highlighted distinctly from the post-crossing BM and deep Newtonian phases.

### Characteristic Frequency Plots

Same layout as shock dynamics, tracking nu_m, nu_c, and nu_M evolution. The frequency ordering determines which spectral regime applies at each time.

### Spectrum Shape Plots (per regime)

Each regime has a dedicated plot showing:

- Upper: F_nu vs nu with break frequencies marked as vertical lines
- Lower: Spectral index beta = d(log F)/d(log nu)

Colored regions indicate frequency segments between breaks. Flat regions in the lower panel confirm correct power-law behavior within each segment.
