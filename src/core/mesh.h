//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <numeric>
#include <span>
#include <vector>

#include "../util/macros.h"
#include "../util/traits.h"
#include "boost/numeric/odeint.hpp"
#include "physics.h"
#include "xtensor/containers/xadapt.hpp"
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
using MaskGrid = xt::xtensor<int, 3>;
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

/// Tunable coefficients for the adaptive theta grid PDF.
struct ThetaGridParams {
    Real core_beam_coeff = 55.0; ///< Extra beam points per log-decade for core beaming
    Real view_beam_coeff = 25.0; ///< Extra beam points per log-decade for view beaming
    Real doppler_alpha = 12.0;   ///< Doppler boost factor in structure term
    Real floor_fraction = 0.1;   ///< Uniform floor weight as fraction of peak_weight
};

/**
 * <!-- ************************************************************************************** -->
 * @class Coord
 * @brief Represents a coordinate system with arrays for phi, theta, and t.
 * @details This class is used to define the computational grid for GRB simulations.
 *          It stores the angular coordinates (phi, theta) and time (t) arrays,
 *          along with derived quantities needed for numerical calculations.
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

    Symmetry symmetry{Symmetry::structured}; ///< Auto-detected symmetry level
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

template <typename Arr = Array>
void log2space(Real lg2_min, Real lg2_max, Real grid_per_order, Arr& array) {
    const Real decades = (lg2_max - lg2_min) / log2_10;
    const size_t n = std::max<size_t>(2, static_cast<size_t>(std::ceil(decades * grid_per_order)));

    array = xt::logspace<Real>(lg2_min, lg2_max, n + 1, 2.0);
}

template <typename Arr>
void logspace_boundary_center(Real lg2_min, Real lg2_max, size_t size, Arr& center, Array& bin_width) {
    center = Arr::from_shape({size});
    bin_width = Array::from_shape({size});
    if (size == 0) {
        return;
    }

    const Real dlg2 = (lg2_max - lg2_min) / static_cast<Real>(size);
    const Real r = std::exp2(dlg2);
    const Real s = std::sqrt(r);
    const Real w = r - 1.;

    Real left = std::exp2(lg2_min);

    for (std::size_t i = 0; i < size; ++i) {
        center(i) = left * s;
        bin_width(i) = left * w;
        left *= r;
    }
}
