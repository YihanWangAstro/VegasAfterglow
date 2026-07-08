//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <span>
#include <vector>

#include "../util/fast-math.h"
#include "../util/macros.h"
#include "../util/traits.h"
#include "physics.h"
#include "xtensor/containers/xadapt.hpp"
#include "xtensor/generators/xbuilder.hpp"
#include "xtensor/views/xview.hpp"

constexpr Real log2_10 = std::numbers::ln10 / std::numbers::ln2;
/**
 * <!-- ************************************************************************************** -->
 * @defgroup Mesh_Utilities Array and Grid Utilities
 * @brief Functions for generating and processing arrays and grids.
 * @details Declares a set of functions for generating and processing arrays and grids. These functions include:
 *          - Converting boundary arrays to center arrays (linear and logarithmic),
 *          - Creating linearly and logarithmically spaced arrays,
 *          - Creating arrays with uniform spacing in cosine,
 *          - Generating arrays of zeros and ones,
 *          - Finding the minimum and maximum of grids,
 *          - Checking if an array is linearly or logarithmically scaled,
 *          - Creating 2D and 3D grids.
 * <!-- ************************************************************************************** -->
 */

/// Type alias for 1D arrays (e.g., time points)
using Array = xt::xtensor<Real, 1>;
/// Type alias for 2D grids (e.g., spatial coordinates)
using MeshGrid = xt::xtensor<Real, 2>;
/// Type alias for 3D grids (e.g., full spatial-temporal data)
using MeshGrid3d = xt::xtensor<Real, 3>;
/// Type alias for 3D boolean grids for masking operations
/// Type alias for 2D index grids
using IndexGrid = xt::xtensor<size_t, 2>;

/**
 * <!-- ************************************************************************************** -->
 * @brief Symmetry level of the computational grid, auto-detected from jet initial conditions.
 * <!-- ************************************************************************************** -->
 */
enum class Symmetry : uint8_t {
    structured,    ///< Full (phi, theta) computation needed
    phi_symmetric, ///< Uniform in phi; compute one phi slice, broadcast across phi
    piecewise,     ///< phi-symmetric + piecewise-constant in theta; compute K representative thetas
    isotropic      ///< Uniform in both phi and theta; compute one point, broadcast everywhere
};

/**
 * <!-- ************************************************************************************** -->
 * @class Coord
 * @brief Represents a coordinate system with arrays for phi, theta, and t.
 * @details This class is used to define the computational grid for GRB simulations.
 *          It stores the angular coordinates (phi, theta) and time (t) arrays,
 *          along with derived quantities needed for numerical calculations.
 *          Lifecycle: instances are built by auto_grid (grid-refinement.h), which fills
 *          theta/phi via the adaptive-grid algorithms, classifies symmetry via
 *          detect_symmetry, and populates t via build_time_grid.
 * <!-- ************************************************************************************** -->
 */
class Coord {
  public:
    /// Default constructor
    Coord() = default;

    Array phi;          ///< Array of azimuthal angles (phi) in radians
    Array theta;        ///< Array of polar angles (theta) in radians
    MeshGrid3d t;       ///< Array of engine time points
    Real theta_view{0}; ///< Viewing angle
    Real phi_view{0};   ///< Viewing angle

    /// True when phi covers [0, pi] only: axisymmetric jets viewed off-axis are
    /// mirror-symmetric about the jet-observer plane, so the phi integral over
    /// [0, 2pi] is exactly twice the integral over [0, pi]. The doubling lives
    /// in the observer's phi weights.
    bool phi_mirrored{false};

    Symmetry symmetry{Symmetry::structured}; ///< Auto-detected symmetry level
    bool spreading{false};                   ///< Jet lateral spreading (set by detect_symmetry)
    std::vector<size_t> theta_reps;          ///< Representative theta indices (contiguous groups)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the dimensions of the coordinate arrays
     * @return Tuple containing (n_phi, n_theta, n_t)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.shape()[2]); }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Detects the symmetry level of the grid based on jet and medium properties.
     * @details Compares adjacent theta cells on initial conditions (Gamma0, eps_k, sigma0,
     *          and time-dependent injection). Sets symmetry and theta_reps accordingly.
     * @param jet The jet/ejecta object
     * @param medium The medium object
     * @param t_min Minimum observer time
     * @param t_max Maximum observer time
     * @param z Redshift
     * @param t_resol Time resolution for sampling injection checks
     * <!-- ************************************************************************************** -->
     */
    template <typename Ejecta, typename Medium>
    void detect_symmetry(Ejecta const& jet, Medium const& medium, Real t_min, Real t_max, Real z, Real t_resol);
};

template <typename Ejecta, typename Medium>
void Coord::detect_symmetry(Ejecta const& jet, Medium const& medium, Real t_min, Real t_max, Real z, Real t_resol) {
    const size_t theta_size = theta.size();
    spreading = jet.spreading;

    if (jet.spreading || !medium.isotropic) {
        symmetry = Symmetry::structured;
        theta_reps.resize(theta_size);
        std::iota(theta_reps.begin(), theta_reps.end(), size_t(0));
        return;
    }

    const Real phi0 = phi(0);
    theta_reps.clear();
    theta_reps.reserve(theta_size);
    theta_reps.push_back(0);

    // Logspaced engine time for time-dependent injection checks
    const size_t t_check = std::max<size_t>(static_cast<size_t>(std::log10(t_max / t_min) * t_resol), 24);
    const Array temp_t = xt::logspace(std::log10(t_min / (1 + z)), std::log10(t_max / (1 + z)), t_check);

    auto jet_ic_differs = [&](size_t ja, size_t jb) {
        const Real theta_a = theta(ja);
        const Real theta_b = theta(jb);
        if (jet.eps_k(phi0, theta_a) != jet.eps_k(phi0, theta_b)) {
            return true;
        }
        if (jet.Gamma0(phi0, theta_a) != jet.Gamma0(phi0, theta_b)) {
            return true;
        }
        if constexpr (HasSigma<Ejecta>) {
            if (jet.sigma0(phi0, theta_a) != jet.sigma0(phi0, theta_b)) {
                return true;
            }
        }
        if constexpr (HasDedt<Ejecta>) {
            for (size_t k = 0; k < t_check; ++k) {
                if (jet.deps_dt(phi0, theta_a, temp_t(k)) != jet.deps_dt(phi0, theta_b, temp_t(k))) {
                    return true;
                }
            }
        }
        if constexpr (HasDmdt<Ejecta>) {
            for (size_t k = 0; k < t_check; ++k) {
                if (jet.dm_dt(phi0, theta_a, temp_t(k)) != jet.dm_dt(phi0, theta_b, temp_t(k))) {
                    return true;
                }
            }
        }
        return false;
    };

    for (size_t j = 1; j < theta_size; ++j) {
        if (jet_ic_differs(j - 1, j)) {
            theta_reps.push_back(j);
        }
    }

    if (theta_reps.size() == 1) {
        symmetry = Symmetry::isotropic;
    } else if (theta_reps.size() < theta_size) {
        symmetry = Symmetry::piecewise;
    } else {
        symmetry = Symmetry::phi_symmetric;
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using linear interpolation
 * @tparam Arr1 Type of the input array
 * @tparam Arr2 Type of the output array
 * @param boundary Input array of boundary values
 * @param center Output array of center values
 * @details This is used to compute cell-centered values from cell-boundary values.
 * <!-- ************************************************************************************** -->
 */
template <typename Arr1, typename Arr2>
void boundary_to_center(Arr1 const& boundary, Arr2& center);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using geometric mean (logarithmic interpolation)
 * @tparam Arr1 Type of the input array
 * @tparam Arr2 Type of the output array
 * @param boundary Input array of boundary values
 * @param center Output array of center values
 * @details This is used when working with logarithmically scaled quantities.
 * <!-- ************************************************************************************** -->
 */
template <typename Arr1, typename Arr2>
void boundary_to_center_log(Arr1 const& boundary, Arr2& center);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using linear interpolation
 * @param boundary Input array of boundary values
 * @return Array of center values
 * <!-- ************************************************************************************** -->
 */
Array boundary_to_center(Array const& boundary);

/**
 * <!-- ************************************************************************************** -->
 * @brief Converts boundary values to center values using geometric mean
 * @param boundary Input array of boundary values
 * @return Array of center values
 * <!-- ************************************************************************************** -->
 */
Array boundary_to_center_log(Array const& boundary);

template <typename Arr1, typename Arr2>
void boundary_to_center(Arr1 const& boundary, Arr2& center) {
    if (boundary.size() < 2 || center.size() + 1 != boundary.size()) {
        return;
    }
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = 0.5 * (boundary[i] + boundary[i + 1]);
    }
}

template <typename Arr1, typename Arr2>
void boundary_to_center_log(Arr1 const& boundary, Arr2& center) {
    if (boundary.size() < 2 || center.size() + 1 != boundary.size()) {
        return;
    }
    for (size_t i = 0; i < center.size(); ++i) {
        center[i] = std::sqrt(boundary[i] * boundary[i + 1]);
    }
}

template <typename Arr = Array>
void log2space(Real lg2_min, Real lg2_max, Real grid_per_order, Arr& array) {
    const Real decades = (lg2_max - lg2_min) / log2_10;
    const size_t n = std::max<size_t>(2, static_cast<size_t>(std::ceil(decades * grid_per_order)));

    // Direct exp2 lattice: xt::logspace(base 2) routes every node through libm
    // pow(2, x), which dominates the cost of building these small grids.
    array = Arr::from_shape({n + 1});
    const Real step = (lg2_max - lg2_min) / static_cast<Real>(n);
    for (size_t i = 0; i <= n; ++i) {
        array(i) = std::exp2(lg2_min + step * static_cast<Real>(i));
    }
}

template <typename Arr>
void logspace_boundary_center(Real lg2_min, Real lg2_max, size_t size, Arr& center, Array& bin_width) {
    center = Arr::from_shape({size});
    bin_width = Array::from_shape({size});
    if (size == 0) {
        return;
    }

    const Real dlg2 = (lg2_max - lg2_min) / static_cast<Real>(size);
    const Real r = fast_exp2(dlg2);
    const Real s = std::sqrt(r);
    const Real w = r - 1.;

    Real left = fast_exp2(lg2_min);

    for (std::size_t i = 0; i < size; ++i) {
        center(i) = left * s;
        bin_width(i) = left * w;
        left *= r;
    }
}
