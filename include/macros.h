//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <limits>
#include <numeric>
/**
 * <!-- ************************************************************************************** -->
 * @namespace con
 * @brief Physical and astronomical constants used throughout the simulation.
 * @details Contains constant definitions used throughout the simulation. These constants define unit
 *          conversions (e.g., seconds, centimeters, grams), physical constants (e.g., speed of light, masses,
 *          Planck's constant), energy units, and cosmological parameters.
 * <!-- ************************************************************************************** -->
 */

/// Type alias for floating-point precision used throughout the code
using Real = double;

/// Literal operator for converting long double literals to Real
constexpr Real operator"" _r(long double x) {
    return static_cast<Real>(x);
}

/// Literal operator for converting integer literals to Real
constexpr Real operator"" _r(unsigned long long x) {
    return static_cast<Real>(x);
}

/**
 * <!-- ************************************************************************************** -->
 * @namespace unit
 * @brief Unit conversion constants
 * @details Contains physical unit conversion factors normalized to the internal code units
 * <!-- ************************************************************************************** -->
 */
namespace unit {
    constexpr Real len = 1.5e13;                        ///< Length unit in cm
    constexpr Real cm = 1 / len;                        ///< Centimeter in code units
    constexpr Real sec = 3e10 / len;                    ///< Second in code units
    constexpr Real cm2 = cm * cm;                       ///< Square centimeter in code units
    constexpr Real cm3 = cm * cm * cm;                  ///< Cubic centimeter in code units
    constexpr Real g = 1 / 2e33;                        ///< Gram in code units
    constexpr Real kg = 1000 * g;                       ///< Kilogram in code units
    constexpr Real Gauss = 8.66e-11 / sec;              ///< Gauss (magnetic field) in code units
    constexpr Real Hz = 1 / sec;                        ///< Hertz (frequency) in code units
    constexpr Real erg = g * cm * cm / sec / sec;       ///< Erg (energy) in code units
    constexpr Real M_sun = 2e33 * g;                    ///< Solar mass in code units
    constexpr Real eV = 1.60218e-12 * erg;              ///< Electron volt in code units
    constexpr Real keV = 1e3 * eV;                      ///< Kilo-electron volt in code units
    constexpr Real MeV = 1e6 * eV;                      ///< Mega-electron volt in code units
    constexpr Real GeV = 1e9 * eV;                      ///< Giga-electron volt in code units
    constexpr Real TeV = 1e12 * eV;                     ///< Tera-electron volt in code units
    constexpr Real deg = 3.14159265358979323846 / 180;  ///< Degree in radians
    constexpr Real flux_cgs = erg / cm2 / sec;          ///< Flux in CGS units
    constexpr Real flux_den_cgs = erg / cm2 / sec / Hz; ///< Flux density in CGS units
    constexpr Real Jy = 1e-23 * erg / cm2 / sec / Hz;   ///< Jansky (radio flux density) in code units
    constexpr Real mJy = 1e-3 * Jy;                     ///< Milli-Jansky in code units
    constexpr Real uJy = 1e-6 * Jy;                     ///< Micro-Jansky in code units
    constexpr Real m = 100 * cm;                        ///< Meter in code units
    constexpr Real km = 1000 * m;                       ///< Kilometer in code units
    constexpr Real au = 1.5e13 * cm;                    ///< Astronomical Unit in code units
    constexpr Real pc = 2.06265e5 * au;                 ///< Parsec in code units
    constexpr Real kpc = 1000 * pc;                     ///< Kiloparsec in code units
    constexpr Real Mpc = 1e6 * pc;                      ///< Megaparsec in code units
    constexpr Real hr = 3600 * sec;                     ///< Hour in code units
    constexpr Real day = 24 * hr;                       ///< Day in code units
    constexpr Real yr = 365.2425 * day;                 ///< Year in code units
    constexpr Real Myr = 1e6 * yr;                      ///< a Million years in code units
    constexpr Real Gyr = 1e9 * yr;                      ///< a Billion years in code units
} // namespace unit

/**
 * <!-- ************************************************************************************** -->
 * @namespace con
 * @brief Physical and astronomical constants
 * @details Contains fundamental physical constants and cosmological parameters
 * <!-- ************************************************************************************** -->
 */
namespace con {
    constexpr Real c = 1;                                ///< Speed of light (normalized to 1 in code units)
    constexpr Real c2 = c * c;                           ///< Speed of light squared
    constexpr Real mp = 1.67e-24 * unit::g;              ///< Proton mass in code units
    constexpr Real me = mp / 1836;                       ///< Electron mass in code units
    constexpr Real mec2 = me * c2;                       ///< Electron rest energy in code units
    constexpr Real mpc2 = mp * c2;                       ///< Proton rest energy in code units
    constexpr Real h = 6.63e-27 * unit::erg * unit::sec; ///< Planck constant in code units
    constexpr Real e = 4.8e-10 / 4.472136e16 / 5.809475e19 / unit::sec; ///< Elementary charge in code units
    constexpr Real e2 = e * e;                                          ///< Elementary charge squared
    constexpr Real e3 = e2 * e;                                         ///< Elementary charge cubed
    constexpr Real pi = 3.14159265358979323846;                         ///< Pi
    constexpr Real sigmaT = 6.65e-25 * unit::cm * unit::cm;             ///< Thomson cross-section in code units
    constexpr Real Omega_m = 0.27;                                      ///< Matter density parameter in ΛCDM cosmology
    constexpr Real Omega_L = 0.73;                                ///< Dark energy density parameter in ΛCDM cosmology
    constexpr Real H0 = 67.66 * unit::km / unit::sec / unit::Mpc; ///< Hubble constant in code units
    constexpr Real Gamma_cut = 1 + 1e-6;                          ///< Cutoff Lorentz factor value
    constexpr Real inf = std::numeric_limits<Real>::infinity();   ///< Infinity value
    constexpr Real sigma_cut = 1e-6;                              ///< Cutoff magnetization value
    constexpr Real min_obs_time = 0.1 * unit::sec;                ///< Minimum observer time
} // namespace con
