# Dead code removed 2026-07-02

Functions and commented-out blocks with no production callers, removed from the
live tree during the consistency cleanup. `PowerLawSyn` and
`SynElectrons::compute_N_gamma` were deliberately **kept** in the tree:
`PowerLawSyn` is a selectable alternative photon model behind the `SynPhotons`
alias, and `compute_N_gamma`'s unit tests exercise the live `compute_spectrum`
internals.

Unit tests that existed solely to test the removed functions were removed with
them (from `tests/cpp/test_mesh.cpp`, `test_shock_physics.cpp`,
`test_synchrotron.cpp`).

---

## `get_post_cross_g` — src/dynamics/reverse-shock.tpp (no callers, no tests)

```cpp
/**
 * @internal
 * @brief Calculates the power-law index for post-crossing four-velocity evolution.
 * @details The index transitions from g_low=1.5 for low relative Lorentz factors to g_high=3.5
 *          for high relative Lorentz factors (Blandford-McKee limit).
 * @param gamma_rel Relative Lorentz factor
 * @param k Medium power law index (default: 0)
 * @return The power-law index for velocity evolution
 */
inline Real get_post_cross_g(Real gamma_rel, Real /*k*/ = 0) {
    constexpr Real g_low = 1.5;  // k is the medium power law index
    constexpr Real g_high = 3.5; // Blandford-McKee limit// TODO: need to be modified for non ISM medium
    const Real p = std::sqrt(std::sqrt(std::max(gamma_rel - 1, 0.0)));
    return g_low + (g_high - g_low) * p / (1 + p);
}
```

## `compute_adiabatic_cooling_rate` — src/dynamics/shock-physics.h (superseded by `compute_adiabatic_cooling_rate2`)

```cpp
/**
 * @brief Computes the adiabatic cooling rate.
 * @param ad_idx Adiabatic index
 * @param r Radius
 * @param Gamma Lorentz factor
 * @param u Internal energy density
 * @param drdt Rate of change of radius
 * @param dGammadt Rate of change of the Lorentz factor
 * @return The adiabatic cooling rate
 */
inline constexpr Real compute_adiabatic_cooling_rate(Real ad_idx, Real r, Real Gamma, Real u, Real drdt,
                                                     Real dGammadt) noexcept {
    return -(ad_idx - 1) * (3 * drdt / r - dGammadt / Gamma) * u;
}
```

## `logspace_center` — src/core/mesh.h (production uses `logspace_boundary_center`)

```cpp
template <typename Arr = Array>
void logspace_center(Real lg2_min, Real lg2_max, size_t size, Arr& center) {
    center = Arr::from_shape({size});
    if (size == 0) {
        return;
    }

    const Real dlg2 = (lg2_max - lg2_min) / static_cast<Real>(size);
    const Real r = std::exp2(dlg2);
    const Real s = std::sqrt(r);
    Real left = std::exp2(lg2_min);

    for (std::size_t i = 0; i < size; ++i) {
        center(i) = left * s;
        left *= r;
    }
}
```

## `compute_gamma_peak` (both overloads) — src/radiation/synchrotron.h / .cpp

```cpp
/**
 * @brief Determines the electron Lorentz factor at which the number density peaks.
 * @details Based on the relative ordering of absorption, minimum, and cooling Lorentz factors.
 * @param gamma_a Absorption Lorentz factor
 * @param gamma_m Minimum electron Lorentz factor
 * @param gamma_c Cooling electron Lorentz factor
 * @return Peak Lorentz factor
 */
Real compute_gamma_peak(Real gamma_a, Real gamma_m, Real gamma_c) {
    const Real gamma_peak = std::min(gamma_m, gamma_c);
    if (gamma_a > gamma_c) {
        return gamma_a;
    } else {
        return gamma_peak;
    }
}

/**
 * @brief Determines the peak Lorentz factor directly from a SynElectrons object.
 * @details Convenient wrapper around the three-parameter version.
 * @param e Synchrotron electron object
 * @return Peak Lorentz factor
 */
Real compute_gamma_peak(SynElectrons const& e) {
    return compute_gamma_peak(e.gamma_a, e.gamma_m, e.gamma_c);
}
```

---

## Commented-out blocks

### Alternate trapezoid-sum `estimate_t_dec` body — src/core/grid-refinement.h

```cpp
    /*const Real gamma = jet.Gamma0(phi, theta);
    if (gamma <= 1) {
        return con::inf;
    }
    const Real beta = physics::relativistic::gamma_to_beta(gamma);
    Real m_jet = jet.eps_k(phi, theta) / (gamma * con::c2);
    if constexpr (HasSigma<Ejecta>) {
        m_jet /= (1.0 + jet.sigma0(phi, theta));
    }
    const Real target = m_jet / gamma;

    constexpr size_t N = 256;
    const Real u_min = std::log(1e-3);
    const Real u_max = u_min + 40 * std::log(10.0);
    const Real du = (u_max - u_min) / N;

    Real mass = 0;
    Real r_prev = std::exp(u_min);
    Real f_prev = medium.rho(phi, theta, r_prev) * r_prev * r_prev;

    for (size_t i = 1; i <= N; ++i) {
        const Real r_i = std::exp(u_min + i * du);
        const Real f_i = medium.rho(phi, theta, r_i) * r_i * r_i;
        const Real dr = r_i - r_prev;
        mass += 0.5 * (f_prev + f_i) * dr;

        if (mass >= target) {
            const Real r_dec = r_prev + (target - (mass - 0.5 * (f_prev + f_i) * dr)) / f_i;
            return r_dec * (1 - beta) / (beta * con::c);
        }
        f_prev = f_i;
        r_prev = r_i;
    }
    return std::exp(u_max) * (1 - beta) / (beta * con::c);*/
```

### Disabled reverse-shock jet spreading — src/dynamics/reverse-shock.tpp (`generate_shock_pair`)

```cpp
        // Real theta_s =
        //     jet_spreading_edge(jet, medium, coord.phi(i), coord.theta.front(), coord.theta.back(), coord.t.front());
```

### Alternate bulirsch-stoer stepper — src/dynamics/reverse-shock.tpp (`grid_solve_shock_pair`)

```cpp
    // auto stepper = bulirsch_stoer_dense_out<typename Eqn::State>{rtol, rtol};
```

### Alternate SimpleShockEqn driver — src/dynamics/forward-shock.tpp (`generate_fwd_shock`)

```cpp
            //auto eqn = SimpleShockEqn(medium, jet, coord.phi(i), coord.theta(j), rad_params, theta_s);
```

### `MaskGrid required` remnants — src/dynamics/shock.h / shock.cpp

```cpp
    // shock.h (member)
    // MaskGrid required;        ///< Grid points actually required for final flux calculation

    // shock.cpp (constructor init + resize)
    // required({phi_size, theta_size, t_size}, 1),    // Initialize the required grid with all-true
    // required.resize({phi_size, theta_size, t_size});
    // required.fill(1);
```

## `MaskGrid` alias — src/core/mesh.h (zero uses anywhere; companion of the MaskGrid remnants above)

```cpp
using MaskGrid = xt::xtensor<int, 3>; ///< 3D grid alias for masks
```
