//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <utility>

#include "../util/macros.h"
#include "../util/utilities.h"
/**
 * <!-- ************************************************************************************** -->
 * @class Medium
 * @brief Represents the generic medium or any user-defined surrounding medium that the GRB jet interacts with.
 * @details The class provides methods to compute the density (rho) as a function of position (phi, theta, r).
 * <!-- ************************************************************************************** -->
 */
class Medium {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with density function
     * @param rho Density function
     * <!-- ************************************************************************************** -->
     */
    explicit Medium(TernaryFunc rho) noexcept : rho(std::move(rho)) {}

    Medium() = default;

    /// Density function that returns the mass density at a given position (phi, theta, r)
    /// The function is initialized to zero by default
    TernaryFunc rho{func::zero_3d};
};

/**
 * <!-- ************************************************************************************** -->
 * @class ISM
 * @brief Implements a uniform interstellar medium (ISM) with constant density.
 * @details Provides methods to compute density at any position.
 *          The ISM is characterized by the particle number density n_ism.
 * <!-- ************************************************************************************** -->
 */
class ISM {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with particle number density in cm^-3
     * @param n_ism Particle number density in cm^-3
     * <!-- ************************************************************************************** -->
     */
    explicit ISM(Real n_ism) noexcept : rho_(n_ism * con::mp) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Return density at a given position (constant everywhere)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance (unused)
     * @return Constant density value
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] inline Real rho(Real phi, Real theta, Real r) const noexcept { return rho_; }

  private:
    Real rho_{0}; ///< Mass density (particle number density × proton mass)
};

/**
 * <!-- ************************************************************************************** -->
 * @class Wind
 * @brief Implements a stellar wind medium with density proportional to 1/r².
 * @details Provides methods to compute density at any position.
 *          The wind is characterized by the wind parameter A_star.
 * <!-- ************************************************************************************** -->
 */
class Wind {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor: Initialize with wind parameter A_star (in standard units)
     * @param A_star Wind density parameter in standard units
     * @param n_ism number density of ISM at large radii
     * @param n0 number density of wind at small radii
     * <!-- ************************************************************************************** -->
     */
    explicit Wind(Real A_star, Real n_ism = 0, Real n0 = con::inf) noexcept
        : A(A_star * 5e11 * unit::g / unit::cm), rho_ism(n_ism * con::mp), r02(A / (n0 * con::mp)) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Return density at given positions (proportional to 1/r²)
     * @param phi Azimuthal angle (unused)
     * @param theta Polar angle (unused)
     * @param r Radial distance
     * @return Density value at radius r (= A/r²)
     * <!-- ************************************************************************************** -->
     */

    [[nodiscard]] inline Real rho(Real phi, Real theta, Real r) const noexcept { return A / (r02 + r * r) + rho_ism; }

  private:
    Real A{0};       ///< Wind density parameter in physical units
    Real rho_ism{0}; ///< ISM density floor
    Real r02{0};     ///< Radius where ISM transitions to
};

/**
 * <!-- ************************************************************************************** -->
 * @namespace evn
 * @brief Provides functions to create different types of ambient medium profiles.
 * @details These functions return lambda functions that compute the density at any given position.
 * <!-- ************************************************************************************** -->
 */
namespace evn {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a uniform interstellar medium (ISM) profile
     * @param n_ism Number density of particles in cm^-3
     * @return function for density
     * <!-- ************************************************************************************** -->
     */
    inline auto ISM(Real n_ism) {
        const Real rho = n_ism * con::mp;

        return [=](Real phi, Real theta, Real r) noexcept { return rho; };
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Creates a stellar wind medium profile
     * @param A_star Wind parameter in standard units
     * @param n_ism Number density of the ISM [cm^-3]
     * @param n0 Number density of inner region [cm^-3]
     * @param k power-law index
     * @return functions for density calculation
     * @details Converts A_star to proper units (A_star * 5e11 g/cm) and returns functions that compute
     *          density = A/r² and mass = A*r, representing a steady-state stellar wind where density
     *          falls off as 1/r²
     * <!-- ************************************************************************************** -->
     */
    inline auto wind(Real A_star, Real n_ism = 0, Real n0 = con::inf, Real k = 2) {
        // Convert A_star to proper units: A_star * 5e11 g/cm
        constexpr Real r0 = 1e17 * unit::cm; // reference radius
        const Real A = A_star * 5e11 * unit::g / unit::cm * std::pow(r0, k - 2);
        const Real rho_ism = n_ism * con::mp;
        const Real r0k = A / (n0 * con::mp);

        // Return a function that computes density = A/r^k
        // This represents a steady-state stellar wind where density falls off as 1/r^k
        return [=](Real phi, Real theta, Real r) noexcept { return A / (r0k + std::pow(r, k)) + rho_ism; };
    }
} // namespace evn
