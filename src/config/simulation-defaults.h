//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cstddef>

/**
 * <!-- ************************************************************************************** -->
 * @file simulation-defaults.h
 * @brief Centralized simulation configuration defaults
 *
 * This header consolidates all tunable simulation parameters in one place.
 * These are numerical/algorithmic parameters, NOT physical constants.
 *
 * Benefits:
 * - Single source of truth for all defaults
 * - Easy to create different presets (fast vs. accurate)
 * - Clear documentation of what each parameter does
 * - No magic numbers scattered across the codebase
 * <!-- ************************************************************************************** -->
 */

namespace defaults {

    /**
     * <!-- ************************************************************************************** -->
     * @brief Physics cutoff thresholds
     *
     * These values determine when certain physical effects are negligible
     * and can be safely ignored for numerical efficiency.
     * <!-- ************************************************************************************** -->
     */
    namespace cutoffs {
        /// Lorentz factor below which the flow is considered non-relativistic
        /// Used to terminate shock evolution when Gamma approaches 1
        inline constexpr double gamma_cut = 1.0 + 1e-6;

        /// Magnetization parameter below which magnetic effects are ignored
        /// sigma = B^2 / (4*pi*rho*c^2) - ratio of magnetic to rest-mass energy
        inline constexpr double sigma_cut = 1e-6;

        inline constexpr double gamma_therm_cut = 1.0 + 1e-6;
    } // namespace cutoffs

    /**
     * <!-- ************************************************************************************** -->
     * @brief Grid generation parameters
     *
     * These control the adaptive mesh refinement and resolution
     * for the (phi, theta, t) computational grid.
     * <!-- ************************************************************************************** -->
     */
    namespace grid {
        // Default grid resolutions, calibrated 2026-07 against the validation
        // gates (per-component mean error < 5%, max < 15%, over the pinning
        // jet families across viewing angles). Forward-shock-only runs hold
        // the gates on a coarser grid than reverse-shock runs: the
        // deceleration-band time lattice and smoother emission allow
        // (0.06, 0.15, 6), while reverse-shock components need
        // (0.06, 0.2, 10) — theta for the wing rows driving both shocks,
        // time for the crossing-resolved light curves (t = 9 passes on-axis
        // but fails off-axis: 8%/34% at theta_obs = 2 theta_c). An explicit
        // `resolutions` argument overrides either set.

        /// Azimuthal (phi) angular resolution factor
        inline constexpr double phi_resolution = 0.06;

        /// Polar (theta) angular resolution factor (forward-shock only)
        inline constexpr double theta_resolution = 0.15;

        /// Temporal resolution factor (forward-shock only)
        inline constexpr double time_resolution = 6.0;

        /// Polar (theta) angular resolution for reverse-shock runs
        inline constexpr double rvs_theta_resolution = 0.2;

        /// Temporal resolution for reverse-shock runs
        inline constexpr double rvs_time_resolution = 10.0;

        /// Minimum number of theta grid points (forward shock)
        inline constexpr std::size_t min_theta_points = 36;

        /// Minimum polar angle (avoids singularity at theta=0)
        inline constexpr double theta_min = 1e-6;
    } // namespace grid

    /**
     * <!-- ************************************************************************************** -->
     * @brief Numerical solver tolerances
     *
     * These control the accuracy of ODE integration and root-finding.
     * <!-- ************************************************************************************** -->
     */
    namespace solver {
        /// Relative tolerance for ODE integration (Runge-Kutta)
        /// (grid-construction CDF sampling; the shock dynamics default is below)
        inline constexpr double ode_rtol = 1e-6;

        /// Default relative tolerance for the shock-dynamics ODE solves.
        /// Calibrated 2026-07: relaxing to 1e-5 looked safe on total-flux
        /// benchmarks but broke the validation suite -- the thin-shell
        /// reverse-shock crossing scalings drift (u slope 1.68 vs 1.5) and the
        /// wind reverse-shock convergence measures inherit a ~5-8% noise
        /// floor. 1e-6 is the loosest tolerance the reverse-shock physics
        /// validates at; the solve runs away entirely below ~1e-3.
        inline constexpr double dynamics_rtol = 1e-6;

        /// Tolerance tightening applied to magnetized (sigma > 0) shell solves,
        /// whose jump conditions are far more tolerance-sensitive than the
        /// sigma = 0 path (see grid_solve_shock_pair for the convergence data).
        inline constexpr double magnetized_rtol_factor = 0.1;

        /// Absolute tolerance for ODE integration
        inline constexpr double ode_atol = 1e-9;

        /// Tolerance for binary search (jet edge finding)
        inline constexpr double binary_search_eps = 1e-9;

        /// Maximum ODE integration steps before giving up
        /// Prevents MCMC from hanging on pathological parameter combinations
        inline constexpr std::size_t max_ode_steps = 100000;
    } // namespace solver

    /**
     * <!-- ************************************************************************************** -->
     * @brief Sampling parameters for integration
     * <!-- ************************************************************************************** -->
     */
    namespace sampling {
        /// Number of sample points for theta integration
        inline constexpr std::size_t theta_samples = 200;
    } // namespace sampling

} // namespace defaults
