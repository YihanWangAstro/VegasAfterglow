//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cstddef>

/**
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
 */

namespace defaults {

    /**
 * @brief Physics cutoff thresholds
 *
 * These values determine when certain physical effects are negligible
 * and can be safely ignored for numerical efficiency.
 */
    namespace cutoffs {
        /// Lorentz factor below which the flow is considered non-relativistic
        /// Used to terminate shock evolution when Gamma approaches 1
        inline constexpr double gamma_cut = 1.0 + 1e-6;

        /// Magnetization parameter below which magnetic effects are ignored
        /// sigma = B^2 / (4*pi*rho*c^2) - ratio of magnetic to rest-mass energy
        inline constexpr double sigma_cut = 1e-6;
    } // namespace cutoffs

    /**
 * @brief Grid generation parameters
 *
 * These control the adaptive mesh refinement and resolution
 * for the (phi, theta, t) computational grid.
 */
    namespace grid {
        /// Azimuthal (phi) angular resolution factor
        /// Higher values = more phi grid points near the viewing direction
        inline constexpr double phi_resolution = 0.3;

        /// Polar (theta) angular resolution factor
        /// Higher values = more theta grid points
        inline constexpr double theta_resolution = 1.0;

        /// Temporal resolution factor
        /// Higher values = more time grid points (logarithmically spaced)
        inline constexpr double time_resolution = 5.0;

        /// Minimum number of theta grid points
        /// Ensures adequate angular coverage even for narrow jets
        inline constexpr std::size_t min_theta_points = 56;

        /// Forward ratio for adaptive grid refinement
        /// Controls distribution of grid points relative to jet edge
        inline constexpr double forward_ratio = 0.3;

        /// Minimum polar angle (avoids singularity at theta=0)
        inline constexpr double theta_min = 1e-6;
    } // namespace grid

    /**
 * @brief Numerical solver tolerances
 *
 * These control the accuracy of ODE integration and root-finding.
 */
    namespace solver {
        /// Relative tolerance for ODE integration (Runge-Kutta)
        inline constexpr double ode_rtol = 1e-6;

        /// Absolute tolerance for ODE integration
        inline constexpr double ode_atol = 1e-9;

        /// Tolerance for binary search (jet edge finding)
        inline constexpr double binary_search_eps = 1e-9;

        /// Tolerance for checking linear/log scale arrays
        inline constexpr double scale_check_tol = 1e-6;
    } // namespace solver

    /**
 * @brief Observer frame parameters
 */
    namespace observer {
        /// Minimum observer time to avoid numerical issues at t=0
        /// Value in seconds (will be converted to code units)
        inline constexpr double min_obs_time_sec = 0.1;
    } // namespace observer

    /**
 * @brief Sampling parameters for integration
 */
    namespace sampling {
        /// Number of sample points for theta integration
        inline constexpr std::size_t theta_samples = 200;
    } // namespace sampling

} // namespace defaults
