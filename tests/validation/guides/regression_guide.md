# Regression Test Report Guide

This document describes the regression validation framework for VegasAfterglow, which verifies that simulation outputs reproduce the expected power-law scaling relations from GRB afterglow theory.

---

## 1. Test Categories

The regression suite validates physical quantities against analytical predictions derived from standard synchrotron afterglow theory. Tests cover both forward and reverse shocks:

**Forward shock:**

1. **Shock dynamics**: Lorentz factor, radius, magnetic field, particle number
2. **Characteristic frequencies**: $\nu_m$ (injection), $\nu_c$ (cooling), $\nu_M$ (maximum)
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

| Phase | Physics |
|-------|---------|
| Coasting | Free expansion, $\Gamma \approx$ const |
| Blandford-McKee | Self-similar deceleration |
| Sedov-Taylor | Non-relativistic, $u < 0.1$ |

### 2.2 Reverse Shock (Thin Shell)

Short engine duration. The reverse shock crosses the ejecta quickly, then evolves as a decaying blast wave.

| Phase | Physics |
|-------|---------|
| Crossing | Reverse shock traversing ejecta |
| Post-crossing | Post-crossing deceleration |
| Sedov-Taylor | Non-relativistic |

### 2.3 Reverse Shock (Thick Shell)

Long engine duration. The reverse shock crosses while the engine is still active, producing different scaling.

| Phase | Physics |
|-------|---------|
| Crossing | Reverse shock traversing during engine activity |
| Post-crossing | Post-crossing deceleration |
| Sedov-Taylor | Non-relativistic |

---

## 3. Expected Scaling Relations

### 3.1 Forward Shock Dynamics

Physical quantities at the shock front scale as power laws with observer time: $Q \propto t^{\alpha}$. The exponents depend on the external medium density profile.

**ISM (constant density n)**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Coasting | 0 | 1 | 0 | 3 |
| Blandford-McKee | -3/8 | 1/4 | -3/8 | 3/4 |
| Sedov-Taylor | -3/5 | 2/5 | -3/5 | 6/5 |

**Wind (density $\propto r^{-2}$)**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Coasting | 0 | 1 | -1 | 1 |
| Blandford-McKee | -1/4 | 1/2 | -3/4 | 1/2 |
| Sedov-Taylor | -1/3 | 2/3 | -1 | 2/3 |

### 3.2 Forward Characteristic Frequencies

The synchrotron spectrum is characterized by break frequencies that evolve with time: $\nu \propto t^{\alpha}$.

**ISM**

| Phase | $\nu_m$ | $\nu_c$ | $\nu_M$ |
|-------|---------|---------|---------|
| Coasting | 0 | -2 | 0 |
| Blandford-McKee | -3/2 | -1/2 | -3/8 |
| Sedov-Taylor | -3/5 | -1/5 | 0 |

**Wind**

| Phase | $\nu_m$ | $\nu_c$ | $\nu_M$ |
|-------|---------|---------|---------|
| Coasting | -1 | -1 | 0 |
| Blandford-McKee | -3/2 | 1/2 | -1/4 |
| Sedov-Taylor | -1 | 1 | 0 |

### 3.3 Reverse Shock Dynamics (Thin Shell)

**ISM**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Crossing | 3/2 | 1 | 0 | 3/2 |
| Post-crossing | -1/4 | 1/4 | — | 0 |
| Sedov-Taylor | -2/5 | 2/5 | — | 0 |

**Wind**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Crossing | 1/2 | 1 | -1 | 1/2 |
| Post-crossing | — | 1/2 | — | 0 |
| Sedov-Taylor | — | 2/3 | — | 0 |

"—" indicates quantities not tested (insufficiently clean power-law for reliable fitting).

### 3.4 Reverse Shock Dynamics (Thick Shell)

**ISM**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Crossing | 1/4 | 1/2 | -1/4 | 1 |
| Post-crossing | -1/4 | 1/4 | — | 0 |
| Sedov-Taylor | -2/5 | 2/5 | — | 0 |

**Wind**

| Phase | $u$ | $r$ | $B$ | $N_p$ |
|-------|-----|-----|-----|-------|
| Crossing | 0 | 1 | -1 | 1 |
| Post-crossing | — | 1/2 | — | 0 |
| Sedov-Taylor | — | 2/3 | — | 0 |

### 3.5 Reverse Shock Characteristic Frequencies

Frequency scaling is validated during the crossing phase. Post-crossing and Sedov-Taylor frequency tests are not performed for reverse shocks since the reverse shock material is no longer freshly shocked.

**Thin Shell — Crossing Phase**

| Medium | $\nu_m$ | $\nu_c$ | $\nu_M$ |
|--------|---------|---------|---------|
| ISM | 0 | -2 | 0 |
| Wind | -1 | 1 | 0 |

**Thick Shell — Crossing Phase**

| Medium | $\nu_m$ | $\nu_c$ | $\nu_M$ |
|--------|---------|---------|---------|
| ISM | — | -1 | -1/4 |
| Wind | -1 | 1 | 0 |

---

## 4. Spectral Regimes

The synchrotron spectrum consists of power-law segments joined at break frequencies. The spectral index $\beta = d(\log F_\nu)/d(\log \nu)$ depends on the ordering of $\nu_a$, $\nu_m$, and $\nu_c$.

### Regime I: $\nu_a < \nu_m < \nu_c$ (Slow Cooling)

The standard slow-cooling spectrum where electrons cool on timescales longer than the dynamical time.

| Frequency Range | $\beta$ |
|-----------------|---------|
| $\nu < \nu_a$ | 2 |
| $\nu_a < \nu < \nu_m$ | 1/3 |
| $\nu_m < \nu < \nu_c$ | $-(p-1)/2$ |
| $\nu > \nu_c$ | $-p/2$ |

### Regime II: $\nu_m < \nu_a < \nu_c$

Self-absorption frequency lies between the injection and cooling frequencies.

| Frequency Range | $\beta$ |
|-----------------|---------|
| $\nu < \nu_m$ | 2 |
| $\nu_m < \nu < \nu_a$ | 5/2 |
| $\nu_a < \nu < \nu_c$ | $-(p-1)/2$ |
| $\nu > \nu_c$ | $-p/2$ |

### Regime III: $\nu_a < \nu_c < \nu_m$ (Fast Cooling)

Fast-cooling regime where electrons radiate most of their energy before the next dynamical time.

| Frequency Range | $\beta$ |
|-----------------|---------|
| $\nu < \nu_a$ | 2 |
| $\nu_a < \nu < \nu_c$ | 1/3 |
| $\nu_c < \nu < \nu_m$ | -1/2 |
| $\nu > \nu_m$ | $-p/2$ |

### Regime IV: $\nu_c < \nu_a < \nu_m$

Self-absorption occurs above the cooling frequency in fast-cooling conditions.

| Frequency Range | $\beta$ |
|-----------------|---------|
| $\nu < \nu_c$ | 2 |
| $\nu_c < \nu < \nu_a$ | 2 |
| $\nu_a < \nu < \nu_m$ | -1/2 |
| $\nu > \nu_m$ | $-p/2$ |

### Regime V: $\nu_c < \nu_m < \nu_a$

Heavily self-absorbed fast-cooling spectrum.

| Frequency Range | $\beta$ |
|-----------------|---------|
| $\nu < \nu_c$ | 2 |
| $\nu_c < \nu < \nu_m$ | 2 |
| $\nu_m < \nu < \nu_a$ | 5/2 |
| $\nu > \nu_a$ | $-p/2$ |

The electron index $p = 2.2$ is used, giving $-(p-1)/2 = -0.6$ and $-p/2 = -1.1$.

---

## 5. Check Summary Interpretation

The HTML report (`tests/report.py`) lists every regression check grouped by category (forward/reverse shock dynamics, frequencies, spectral shapes), allowing rapid identification of any deviations from expected behavior.

**Check Display**

Each check compares simulation and theory on an interval chart:

- Dot: Measured power-law exponent
- Tick: Expected theoretical value
- Band: Allowed tolerance around the expected value

A check passes when |measured - expected| < tolerance; checks with insufficient data are reported as not measurable.

**Tolerances**

- Shock dynamics: 0.1
- Characteristic frequencies: 0.1
- Spectral shapes: 0.15

---

## 6. Evolution and Spectrum Figures

The figures accompanying each category provide diagnostic information for understanding any discrepancies in the check summary.

### Shock Dynamics Figures

Log-log evolution panels of $u$, $r$, $B$, and $N_p$ vs observer time, for the forward shock and for the thin-shell and thick-shell reverse shock separately:

- Shaded bands mark the fitted phase time windows
- Dashed guides carry the expected power-law slope (hover for the fit numbers)

Agreement between the computed curve and the dashed guide inside each window confirms correct implementation.

### Characteristic Frequency Figures

Same layout as shock dynamics, tracking $\nu_m$, $\nu_c$, $\nu_a$, and $\nu_M$ evolution. The frequency ordering determines which spectral regime applies at each time.

### Spectrum Shape Figures (per regime)

Each regime has a dedicated panel showing $F_\nu$ vs $\nu$ with break frequencies ($\nu_a$, $\nu_m$, $\nu_c$) marked as dashed vertical lines. Dashed guides carry the expected spectral index for each frequency segment between breaks.
