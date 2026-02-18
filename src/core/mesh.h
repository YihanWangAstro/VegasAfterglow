//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <algorithm>
#include <cassert>
#include <span>
#include <vector>

#include "../util/macros.h"
#include "boost/numeric/odeint.hpp"

constexpr Real log2_10 = std::numbers::ln10 / std::numbers::ln2;
#include "physics.h"
#include "xtensor/containers/xadapt.hpp"
#include "xtensor/views/xview.hpp"
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
// Type alias for 2D grids
using IndexGrid = xt::xtensor<size_t, 2>;

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

    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the dimensions of the coordinate arrays
     * @return Tuple containing (n_phi, n_theta, n_t)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] auto shape() const { return std::make_tuple(phi.size(), theta.size(), t.shape()[2]); }
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Checks if an array is linearly scaled within a given tolerance
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from linearity (default: 1e-6)
 * @return True if the array is linearly scaled, false otherwise
 * <!-- ************************************************************************************** -->
 */
bool is_linear_scale(Array const& arr, Real tolerance = defaults::solver::scale_check_tol);

/**
 * <!-- ************************************************************************************** -->
 * @brief Checks if an array is logarithmically scaled within a given tolerance
 * @param arr The array to check
 * @param tolerance Maximum allowed deviation from logarithmic scaling (default: 1e-6)
 * @return True if the array is logarithmically scaled, false otherwise
 * <!-- ************************************************************************************** -->
 */
bool is_log_scale(Array const& arr, Real tolerance = defaults::solver::scale_check_tol);

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

/**
 * <!-- ************************************************************************************** -->
 * @brief Constructs a coordinate grid (Coord) for shock evolution
 * @tparam Ejecta Type of the jet/ejecta class
 * @param jet The jet/ejecta object
 * @param t_obs Array of observation times
 * @param theta_cut Maximum theta value to include
 * @param theta_view Viewing angle
 * @param z Redshift
 * @param phi_resol
 * @param theta_resol
 * @param t_resol
 * @param is_axisymmetric Whether the jet is axisymmetric (default: true)
 * @param phi_view Viewing angle (default: 0)
 * @param min_theta_num Minimum number of theta grid points (default: 56)
 * @param fwd_ratio Forward ratio for grid resolution (default: 0.3)
 * @return A Coord object with the constructed grid
 * @details The grid is based on the observation times (t_obs), maximum theta value (theta_cut), and
 *          specified numbers of grid points in phi, theta, and t. The radial grid is logarithmically
 *          spaced between t_min and t_max, and the theta grid is generated linearly.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Coord auto_grid(Ejecta const& jet, Medium const& medium, Array const& t_obs, Real theta_cut, Real theta_view, Real z,
                Real phi_resol = defaults::grid::phi_resolution, Real theta_resol = defaults::grid::theta_resolution,
                Real t_resol = defaults::grid::time_resolution, bool is_axisymmetric = true, Real phi_view = 0,
                size_t min_theta_num = defaults::grid::min_theta_points,
                Real fwd_ratio = defaults::grid::forward_ratio);

/**
 * <!-- ************************************************************************************** -->
 * @brief Finds all sharp edges in the jet Gamma0 profile plus the outermost boundary
 * @tparam Ejecta Type of the jet/ejecta class
 * @param jet The jet/ejecta object
 * @param gamma_cut Lorentz factor cutoff value
 * @param is_axisymmetric Flag for axisymmetric jets
 * @return Sorted vector of edge theta positions; last element is the outermost jet edge
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
std::vector<Real> find_jet_jumps(Ejecta const& jet, Real gamma_cut, bool is_axisymmetric);

/**
 * <!-- ************************************************************************************** -->
 * @brief Determines the edge of the jet where the spreading is strongest
 * @tparam Ejecta Type of the jet/ejecta class
 * @tparam Medium Type of the ambient medium
 * @param jet The jet/ejecta object
 * @param medium The ambient medium object
 * @param phi Azimuthal angle
 * @param theta_min Minimum polar angle to consider
 * @param theta_max Maximum polar angle to consider
 * @param t0 Initial time
 * @return Angle (in radians) where the spreading is strongest
 * @details The spreading strength is measured by the derivative of the pressure with respect to theta,
 *          which is proportional to d((Gamma-1)Gamma rho)/dtheta
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& medium, Real phi, Real theta_min, Real theta_max, Real t0);
//========================================================================================================
//                                  template function implementation
//========================================================================================================
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
    if (size == 0)
        return;

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
    if (size == 0)
        return;

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

/**
 * <!-- ************************************************************************************** -->
 * @brief Builds an adaptive log-spaced grid with explicit break points.
 *
 * Places grid points at endpoints and break frequencies, then fills segments
 * between with approximately `pts_per_decade` log-spaced points. Nearby points
 * are merged automatically based on the local step size.
 *
 * @param lg2_min  Log2 of minimum grid value
 * @param lg2_max  Log2 of maximum grid value
 * @param breaks   Break frequencies (linear scale) to include
 * @param break_weights  Importance weights for each break (top values get refined regions)
 * @param pts_per_decade  Base density of fill points per decade
 * @param grid     Output array of grid point values (resized internally)
 * @param max_refined_breaks  Maximum number of breaks to refine (sorted by weight)
 * @param refine_radius_decades  Half-width of the refined region around each break (in decades)
 * @param refine_factor  Grid density multiplier inside refined regions
 * @return         Number of grid points
 * <!-- ************************************************************************************** -->
 */
inline size_t adaptive_grid_with_breaks(Real lg2_min, Real lg2_max, std::span<const Real> breaks,
                                        std::span<const Real> break_weights, Real pts_per_decade, Array& grid,
                                        size_t max_refined_breaks = 3, Real refine_radius_decades = 0.5,
                                        Real refine_factor = 3.0) {
    constexpr size_t MAX_PTS = 256;
    const Real lg2_per_decade = log2_10;
    const Real refine_radius = refine_radius_decades * lg2_per_decade;
    const Real coarse_step = lg2_per_decade / std::max(pts_per_decade, Real(1e-6));
    const Real fine_step = coarse_step / refine_factor;
    const Real merge_eps = 0.5 * fine_step;

    // Collect valid breaks within range with weights
    struct Break {
        Real lg2;
        Real weight;
    };
    std::vector<Break> valid_breaks;
    valid_breaks.reserve(breaks.size());
    for (size_t i = 0; i < breaks.size(); ++i) {
        const Real b = breaks[i];
        if (!(b > 0) || !std::isfinite(b)) {
            continue;
        }
        const Real lg2_b = std::log2(b);
        if (lg2_b > lg2_min && lg2_b < lg2_max) {
            Real w = 1;
            if (i < break_weights.size()) {
                const Real wi = break_weights[i];
                if (std::isfinite(wi) && wi >= 0) {
                    w = wi;
                }
            }
            valid_breaks.push_back({lg2_b, w});
        }
    }

    // Sort by weight descending; top N get refinement
    std::ranges::sort(valid_breaks, [](const auto& a, const auto& b) { return a.weight > b.weight; });
    const size_t n_refine = std::min(max_refined_breaks, valid_breaks.size());

    // Build grid points with break tracking for merge protection
    struct Pt {
        Real lg2;
        bool is_break;
    };
    std::vector<Pt> pts;
    pts.reserve(MAX_PTS * 2);

    // Coarse uniform grid
    const Real span = std::max(lg2_max - lg2_min, Real(0));
    const size_t n_coarse = (span > 0 && coarse_step > 0)
                                ? std::max<size_t>(1, static_cast<size_t>(std::ceil(span / coarse_step)))
                                : size_t(1);
    for (size_t i = 0; i <= n_coarse; ++i) {
        pts.push_back({lg2_min + span * static_cast<Real>(i) / static_cast<Real>(n_coarse), false});
    }

    // Break anchors
    for (const auto& br : valid_breaks) {
        pts.push_back({br.lg2, true});
    }

    // Refinement around top-weighted breaks
    for (size_t i = 0; i < n_refine; ++i) {
        const Real lo = std::max(lg2_min, valid_breaks[i].lg2 - refine_radius);
        const Real hi = std::min(lg2_max, valid_breaks[i].lg2 + refine_radius);
        const size_t n = std::max<size_t>(1, static_cast<size_t>((hi - lo) / fine_step));
        for (size_t k = 0; k <= n; ++k) {
            pts.push_back({lo + (hi - lo) * static_cast<Real>(k) / static_cast<Real>(n), false});
        }
    }

    // Sort and merge close points, protecting break positions
    std::ranges::sort(pts, [](const auto& a, const auto& b) { return a.lg2 < b.lg2; });

    std::vector<Pt> merged;
    merged.reserve(pts.size());
    for (const auto& p : pts) {
        if (!merged.empty() && p.lg2 - merged.back().lg2 < merge_eps) {
            if (p.is_break && !merged.back().is_break) {
                merged.back() = p;
            } else {
                merged.back().is_break = merged.back().is_break || p.is_break;
            }
        } else {
            merged.push_back(p);
        }
    }

    // Ensure endpoints
    if (merged.empty()) {
        merged.push_back({lg2_min, true});
        merged.push_back({lg2_max, true});
    } else {
        merged.front().lg2 = lg2_min;
        if (merged.size() == 1) {
            merged.push_back({lg2_max, true});
        } else {
            merged.back().lg2 = lg2_max;
        }
    }

    // Hard cap with uniform subsampling
    if (merged.size() > MAX_PTS) {
        std::vector<Pt> capped;
        capped.reserve(MAX_PTS);
        for (size_t i = 0; i < MAX_PTS; ++i) {
            capped.push_back(merged[i * (merged.size() - 1) / (MAX_PTS - 1)]);
        }
        merged.swap(capped);
    }

    grid = Array::from_shape({merged.size()});
    for (size_t i = 0; i < merged.size(); ++i) {
        grid(i) = std::exp2(merged[i].lg2);
    }
    return merged.size();
}

template <typename Ejecta>
std::vector<Real> find_jet_jumps(Ejecta const& jet, Real gamma_cut, [[maybe_unused]] bool is_axisymmetric) {
    // Find all edges in the active jet region (Gamma0 >= gamma_cut).
    // Sharp edges: large fractional jumps (> 50%) between adjacent scan points.
    // Smooth jets: the outermost theta where Gamma0 crosses gamma_cut.
    constexpr size_t n_scan = 512;
    constexpr Real eps = defaults::solver::binary_search_eps;
    const Real theta_lo = defaults::grid::theta_min;
    const Real theta_hi = con::pi / 2;

    if (jet.Gamma0(0, theta_hi) >= gamma_cut) {
        return {theta_hi};
    }

    std::vector<Real> edges;
    edges.reserve(100);
    Real prev_th = theta_lo;
    Real prev_G = jet.Gamma0(0, theta_lo);

    for (size_t j = 1; j < n_scan; ++j) {
        const Real cur_th = theta_lo + (theta_hi - theta_lo) * static_cast<Real>(j) / (n_scan - 1);
        const Real cur_G = jet.Gamma0(0, cur_th);

        if (prev_G >= gamma_cut || cur_G >= gamma_cut) {
            const Real dG = std::abs(cur_G - prev_G);
            const Real scale = std::max(prev_G - 1, cur_G - 1);
            if (scale > 0 && dG > 0.5 * scale) {
                Real lo = prev_th, hi = cur_th;
                while (hi - lo > eps) {
                    const Real mid = 0.5 * (lo + hi);
                    const Real G_mid = jet.Gamma0(0, mid);
                    if (std::abs(G_mid - prev_G) < std::abs(G_mid - cur_G)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                edges.push_back(lo);
            }
        }
        prev_th = cur_th;
        prev_G = cur_G;
    }
    return edges;
}

template <typename Ejecta>
Real find_theta_max(Ejecta const& jet, Real gamma_cut) {
    constexpr size_t n_scan = 512;
    const Real theta_lo = defaults::grid::theta_min;
    const Real theta_hi = con::pi / 2;

    const Real step = (theta_hi - theta_lo) / n_scan;
    for (Real th = theta_hi; th >= theta_lo; th -= step) {
        if (jet.Gamma0(0, th) >= gamma_cut) {
            return th;
        }
    }
    return theta_lo;
}

template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& /*medium*/, Real phi, Real theta_min, Real theta_max,
                        Real /*t0*/) {
    const Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (Real theta = theta_min; theta <= theta_max; theta += step) {
        // Real G = jet.Gamma0(phi, theta);
        // Real beta0 = physics::relativistic::gamma_to_beta(G);
        // Real r0 = beta0 * con::c * t0 / (1 - beta0);
        // Real rho = medium.rho(phi, theta, 0);
        Real th_lo = std::max(theta - step, theta_min);
        Real th_hi = std::min(theta + step, theta_max);
        const Real dG = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);
        // Real drho = (medium.rho(phi, th_hi, r0) - medium.rho(phi, th_lo, r0)) / (th_hi - th_lo);
        const Real dp = dG; //(2 * G - 1) * rho * dG + (G - 1) * G * drho;

        if (dp < dp_min) {
            dp_min = dp;
            theta_s = theta;
        }
    }
    if (dp_min == 0) {
        theta_s = theta_max;
    }

    return theta_s;
}

template <typename Func>
Array inverse_CFD_sampling(Func&& pdf, Real min, Real max, size_t num) {
    using namespace boost::numeric::odeint;
    constexpr Real rtol = defaults::solver::ode_rtol;
    constexpr size_t sample_num = defaults::sampling::theta_samples;
    Array x_i = xt::linspace(min, max, sample_num);
    Array CDF_i = xt::zeros<Real>({sample_num});

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<Real>());
    stepper.initialize(0, min, (max - min) / 1e3);

    for (size_t k = 1; stepper.current_time() <= max;) {
        stepper.do_step(pdf);
        while (k < x_i.size() && stepper.current_time() > x_i(k)) {
            stepper.calc_state(x_i(k), CDF_i(k));
            ++k;
        }
    }

    Array CDF_out = xt::linspace(CDF_i.front(), CDF_i.back(), num);
    Array x_out = xt::zeros<Real>({num});

    for (size_t k = 0; k < num; ++k) {
        for (size_t j = 0; j < sample_num; ++j) {
            if (CDF_out(k) <= CDF_i(j)) {
                if (j == 0) {
                    x_out(k) = x_i(j);
                } else {
                    Real slope = (x_i(j) - x_i(j - 1)) / (CDF_i(j) - CDF_i(j - 1));
                    x_out(k) = x_i(j - 1) + slope * (CDF_out(k) - CDF_i(j - 1));
                }
                break;
            }
        }
    }
    return x_out;
}

template <typename Ejecta>
Array adaptive_theta_grid(Ejecta const& jet, Real theta_min, Real theta_max, size_t theta_num, Real theta_v) {
    auto eqn = [=, &jet](Real const& /*cdf*/, Real& pdf, Real theta) {
        const Real Gamma = jet.Gamma0(0, theta);
        const Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
        // Real D = 1 / (Gamma * (1 - beta * std::cos(theta - theta_v)));
        const Real a = (1 - beta) / (1 - beta * std::cos(theta - theta_v));
        pdf = a * Gamma * std::sqrt((Gamma - 1) * Gamma) * std::sin(theta);
        // pdf = D * std::sin(theta);
    };

    return inverse_CFD_sampling(eqn, theta_min, theta_max, theta_num);
}

inline Array jump_refinement_grid(std::vector<Real> const& jet_jumps, Real theta_min, Real theta_max, Real avg_spacing,
                                  Real theta_resol) {
    std::vector<Real> points;
    points.reserve(150);
    for (auto& jump : jet_jumps) {
        if (jump >= con::pi / 2 - 0.01)
            continue;
        const Real half = 3 * avg_spacing;
        const size_t n = std::max<size_t>(static_cast<size_t>(10 * theta_resol), 2);
        for (size_t i = 1; i <= n; ++i) {
            // Log-spaced offsets: denser near the edge
            const Real frac = static_cast<Real>(i) / static_cast<Real>(n);
            const Real offset = half * frac * frac;
            const Real below = jump - offset;
            const Real above = jump + offset;
            if (below >= theta_min)
                points.push_back(below);
            if (above <= theta_max)
                points.push_back(above);
        }
        if (jump >= theta_min && jump <= theta_max)
            points.push_back(jump);
    }
    std::ranges::sort(points);
    points.erase(std::ranges::unique(points).begin(), points.end());
    return xt::adapt(points, {points.size()});
}

template <typename Ejecta>
Array adaptive_phi_grid(Ejecta const& jet, size_t phi_num, Real theta_v, Real theta_max, bool is_axisymmetric) {
    if ((theta_v == 0 && is_axisymmetric) || theta_v > theta_max) {
        return xt::linspace(0., 2 * con::pi, phi_num);
    } else {
        const Real cos_tv = std::cos(theta_v);
        const Real sin_tv = std::sin(theta_v);

        // Sample a few theta values to find the peak Doppler weight at each phi
        constexpr size_t n_theta_sample = 8;
        const Real theta_min = defaults::grid::theta_min;

        auto eqn = [=, &jet](Real const& /*cdf*/, Real& pdf, Real phi) {
            const Real cos_phi = std::cos(phi);
            Real max_weight = 0;
            for (size_t it = 0; it < n_theta_sample; ++it) {
                const Real theta = theta_min + (theta_max - theta_min) * static_cast<Real>(it) / (n_theta_sample - 1);
                const Real Gamma = jet.Gamma0(phi, theta);
                const Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
                // cos(angle to LOS) = cos(theta)cos(theta_v) + sin(theta)sin(theta_v)cos(phi)
                const Real cos_alpha = std::cos(theta) * cos_tv + std::sin(theta) * sin_tv * cos_phi;
                const Real a = (1 - beta) / (1 - beta * cos_alpha);
                const Real weight = a * Gamma * std::sqrt((Gamma - 1) * Gamma) * std::sin(theta);
                max_weight = std::max(max_weight, weight);
            }
            pdf = max_weight;
        };

        return inverse_CFD_sampling(eqn, 0, 2 * con::pi, phi_num);
    }
}

template <typename Arr>
Array merge_grids(Arr const& arr1, Arr const& arr2) {
    std::vector<Real> result;
    result.reserve(arr1.size() + arr2.size());

    size_t i = 0, j = 0;
    auto add_unique = [&](Real val) {
        if (result.empty() || result.back() != val)
            result.push_back(val);
    };

    while (i < arr1.size() && j < arr2.size()) {
        if (arr1[i] <= arr2[j]) {
            add_unique(arr1[i++]);
            if (arr1[i - 1] == arr2[j])
                j++;
        } else {
            add_unique(arr2[j++]);
        }
    }
    while (i < arr1.size())
        add_unique(arr1[i++]);
    while (j < arr2.size())
        add_unique(arr2[j++]);

    return xt::adapt(result);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Estimates deceleration time for a given angle from ejecta and medium properties.
 * @details Integrates swept mass outward until m_swept = m_jet/Gamma0(phi,theta). Returns lab-frame time.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Real estimate_t_dec(Ejecta const& jet, Medium const& medium, Real phi, Real theta) {
    const Real gamma = jet.Gamma0(phi, theta);
    if (gamma <= 1) {
        return con::inf;
    }
    const Real beta = physics::relativistic::gamma_to_beta(gamma);
    const Real m_jet = jet.eps_k(phi, theta) / (gamma * con::c2);
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
    return std::exp(u_max) * (1 - beta) / (beta * con::c);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Create a refined logarithmic time grid with higher density around a specified time.
 *
 * Constructs a three-segment logspace grid:
 * - Segment 1: [t_start, t_lo] with base resolution
 * - Segment 2: [t_lo, t_hi] with higher resolution (around t_refine)
 * - Segment 3: [t_hi, t_end] with base resolution
 *
 * Falls back to simple logspace if refinement is not applicable.
 *
 * @param t_start Start time of the grid
 * @param t_end End time of the grid
 * @param t_refine Time around which to refine (e.g., deceleration time)
 * @param t_num Total number of grid points
 * @param t_resol Resolution parameter (points per decade)
 * @param refine_lo_factor Lower bound factor for refinement region (t_lo = refine_lo_factor * t_refine)
 * @param refine_hi_factor Upper bound factor for refinement region (t_hi = refine_hi_factor * t_refine)
 * @return Array of time grid points
 * <!-- ************************************************************************************** -->
 */
inline Array refined_time_grid(Real t_start, Real t_end, Real t_refine, size_t t_num, Real t_resol,
                               Real refine_lo_factor = 0.3, Real refine_hi_factor = 10.0) {
    Array result = xt::zeros<Real>({t_num});

    const Real log_s = std::log10(t_start);
    const Real log_e = std::log10(t_end);

    if (t_refine <= t_start || t_refine >= t_end) {
        return xt::logspace(log_s, log_e, t_num);
    }

    const Real t_lo = std::max(refine_lo_factor * t_refine, t_start);
    const Real t_hi = std::min(refine_hi_factor * t_refine, t_end);
    const Real log_lo = std::log10(t_lo);
    const Real log_hi = std::log10(t_hi);
    const Real log_before = log_lo - log_s;
    const Real log_after = log_e - log_hi;

    size_t n1 = (log_before > 0) ? std::max<size_t>(static_cast<size_t>(t_resol * log_before), 2) : 0;
    size_t n3 = (log_after > 0) ? std::max<size_t>(static_cast<size_t>(t_resol * log_after), 2) : 0;
    size_t shared = (n1 > 0 ? 1 : 0) + (n3 > 0 ? 1 : 0);
    size_t n2_raw = (n1 + n3 + 2 <= t_num + shared) ? (t_num + shared - n1 - n3) : 0;

    if (n2_raw < 2) {
        return xt::logspace(log_s, log_e, t_num);
    }

    size_t n2 = n2_raw;
    size_t idx = 0;

    // Segment 1: [t_start, t_lo]
    if (n1 > 0) {
        auto seg1 = xt::logspace(log_s, log_lo, n1);
        for (size_t k = 0; k < n1; ++k)
            result(idx++) = seg1(k);
    }

    // Segment 2: [t_lo, t_hi] - skip first point if it overlaps with seg1
    {
        auto seg2 = xt::logspace(log_lo, log_hi, n2);
        size_t k0 = (n1 > 0) ? 1 : 0;
        for (size_t k = k0; k < n2; ++k)
            result(idx++) = seg2(k);
    }

    // Segment 3: [t_hi, t_end] - skip first point (overlaps with seg2)
    if (n3 > 0) {
        auto seg3 = xt::logspace(log_hi, log_e, n3);
        for (size_t k = 1; k < n3; ++k)
            result(idx++) = seg3(k);
    }

    return result;
}

template <typename Ejecta, typename Medium>
Coord auto_grid(Ejecta const& jet, Medium const& medium, Array const& t_obs, Real theta_cut, Real theta_view, Real z,
                Real phi_resol, Real theta_resol, Real t_resol, bool is_axisymmetric, Real /*phi_view*/,
                size_t min_theta_num, Real fwd_ratio) {
    Coord coord;
    coord.theta_view = theta_view;

    const auto jet_jumps = find_jet_jumps(jet, con::Gamma_cut, is_axisymmetric);
    const Real jet_edge = find_theta_max(jet, con::Gamma_cut);
    if (!jet_jumps.empty()) {
        jet_edge = std::max(jet_edge, jet_jumps.back());
    }
    Real theta_min = defaults::grid::theta_min;
    Real theta_max = std::min(jet_edge, theta_cut);

    size_t theta_num = min_theta_num + static_cast<size_t>((theta_max - theta_min) * 180 / con::pi * theta_resol);
    const size_t uniform_theta_num = static_cast<size_t>(static_cast<Real>(theta_num) * fwd_ratio);
    size_t adaptive_theta_num = theta_num - uniform_theta_num;

    const Array uniform_theta = xt::linspace(theta_min, theta_max, uniform_theta_num);
    const Array adaptive_theta = adaptive_theta_grid(jet, theta_min, theta_max, adaptive_theta_num, theta_view);

    const Real avg_spacing = (theta_max - theta_min) / theta_num;
    const Array jump_theta = jump_refinement_grid(jet_jumps, theta_min, theta_max, avg_spacing, theta_resol);
    coord.theta = merge_grids(merge_grids(uniform_theta, adaptive_theta), jump_theta);

    // coord.theta = uniform_theta;
    const size_t phi_num = std::max<size_t>(static_cast<size_t>(360 * phi_resol), 1);

    coord.phi = adaptive_phi_grid(jet, phi_num, theta_view, theta_max, is_axisymmetric);

    if (t_obs.size() == 0) {
        assert(false && "auto_grid: t_obs is empty");
        return coord;
    }

    const Real t_max = *std::ranges::max_element(t_obs);
    const Real t_min = *std::ranges::min_element(t_obs);

    // Compute base time grid size + extra points for refinement around t_dec
    const size_t base_t_num = std::max<size_t>(static_cast<size_t>(std::log10(t_max / t_min) * t_resol), 24);
    constexpr Real refine_ratio = 2.0;
    constexpr Real refine_lo_factor = 0.3;
    constexpr Real refine_hi_factor = 10.0;
    const size_t refine_extra =
        static_cast<size_t>((refine_ratio - 1.0) * t_resol * std::log10(refine_hi_factor / refine_lo_factor));
    const size_t t_num = base_t_num + refine_extra;

    const size_t theta_size = coord.theta.size();
    const size_t phi_size_needed = is_axisymmetric ? 1 : coord.phi.size();
    const Real theta_v_max = coord.theta.back() + theta_view;

    coord.t = xt::zeros<Real>({phi_size_needed, theta_size, t_num});
    for (size_t i = 0; i < phi_size_needed; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            const Real b = physics::relativistic::gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
            const Real t_start =
                std::max<Real>(0.99 * t_min * (1 - b) / (1 - std::cos(theta_v_max) * b) / (1 + z), 1e-2 * unit::sec);
            const Real t_end = 1.01 * t_max / (1 + z);
            const Real t_dec = estimate_t_dec(jet, medium, coord.phi(i), coord.theta(j));

            xt::view(coord.t, i, j, xt::all()) =
                refined_time_grid(t_start, t_end, t_dec, t_num, t_resol, refine_lo_factor, refine_hi_factor);
        }
    }

    return coord;
}
