//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/
#pragma once

#include <algorithm>
#include <array>
#include <span>
#include <stdexcept>
#include <vector>

#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "mesh.h"
#include "quadrature.h"
#include "xtensor/core/xnoalias.hpp"

/**
 * @struct SkyImageResult
 * @brief Result of a batched sky image computation: 3D surface brightness map with shared extent and pixel solid angle.
 */
struct SkyImageResult {
    MeshGrid3d image;           ///< [n_frames, npixel, npixel] surface brightness (code units)
    std::array<Real, 4> extent; ///< {x_min, x_max, y_min, y_max} angular extent (rad)
    Real pixel_solid_angle{0};  ///< Pixel solid angle (sr)
};

/// Internal cell grid for Gaussian splatting (sky position, luminosity, kernel width per cell).
struct SplatGrid {
    MeshGrid x, y, L, sigma;
    size_t n_phi{0}, n_theta{0};
    Real xmin{}, xmax{}, ymin{}, ymax{};

    SplatGrid() = default;
    SplatGrid(size_t n_phi, size_t n_theta);

    /// Reset all arrays and bounds for reuse across frames.
    void reset();

    /// Compute per-cell Gaussian σ from neighbor distances.
    void compute_sigma();

    /// Render the splat grid into a frame of the 3D result image.
    void render(MeshGrid3d& image, size_t frame, Real pix_size, bool azimuthal_avg, bool mirror_y) const;
};

/// One cell's per-node rows along the source-time lattice (all log2). Raw row pointers:
/// interleaved output stores keep the compiler from hoisting member-tensor strides.
/// Every flux path evaluates node fluxes through this struct (single-sourced math).
template <typename Photon>
struct CellRow {
    Real const* t;    ///< log2 arrival time per node
    Real const* dop;  ///< log2 Doppler factor per node
    Real const* geom; ///< log2 geometry factor per node
    Photon* ph;       ///< comoving spectrum per node

    /// log2 observed luminosity of node `node` at source-frame log2 frequency `lg2_nu_src`.
    [[nodiscard]] Real node_flux(size_t node, Real lg2_nu_src) const noexcept {
        return ph[node].compute_log2_I_nu(lg2_nu_src - dop[node]) + geom[node];
    }

    /// node_flux across a whole frequency row into dst; dst.size() == lg2_nu_src.size().
    /// Per-node terms are hoisted to stay in registers across the opaque spectrum calls.
    void node_flux_row(size_t node, std::span<const Real> lg2_nu_src, std::span<Real> dst) const noexcept {
        const Real lg2_dop = dop[node];
        const Real lg2_geom = geom[node];
        auto& spectrum = ph[node];
        for (size_t l = 0; l < dst.size(); l++) {
            dst[l] = spectrum.compute_log2_I_nu(lg2_nu_src[l] - lg2_dop) + lg2_geom;
        }
    }
};

/**
 * <!-- ************************************************************************************** -->
 * @class Observer
 * @brief Represents an observer in the GRB afterglow simulation.
 * @details This class handles the calculation of observed quantities such as specific flux, integrated flux,
 *          and spectra. It accounts for relativistic effects (Doppler boosting), cosmological effects (redshift),
 *          and geometric effects (solid angle). The observer can be placed at any viewing angle relative to the
 *          jet axis.
 * <!-- ************************************************************************************** -->
 */
class Observer {
  public:
    /// Default constructor
    Observer() = default;

    /// Grid of observation times
    MeshGrid3d time;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Sets up the Observer for flux calculation.
     * @details Initializes the observation time and Doppler factor grids, as well as the emission surface.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * @param luminosity_dist Luminosity distance to the source
     * @param redshift Redshift of the source
     * <!-- ************************************************************************************** -->
     */
    void observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at a single observed frequency
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Observed frequency
     * @param photons Parameter pack of photon grid objects
     * @details Returns the specific flux (as an Array) for a single observed frequency (nu_obs) by computing the
     *          specific flux over the observation times.
     * @return Array of specific flux values at each observation time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at multiple observed frequencies
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Array of observed frequencies
     * @param photons Parameter pack of photon grid objects
     * @return 2D grid of specific flux values (frequency × time)
     * @details Returns the specific flux (as a MeshGrid) for multiple observed frequencies (nu_obs) by computing
     *          the specific flux for each frequency and assembling the results into a grid. This method accounts for
     *          relativistic beaming and cosmological effects.
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    MeshGrid specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the specific flux at a single observed frequency for multiple observation times
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param nu_obs Observed frequency
     * @param photons Parameter pack of photon grid objects
     * @return Array of specific flux values at each observation time for a single observed frequency
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array specific_flux_series(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the integrated flux over a frequency band
     * @tparam PhotonGrid Types of photon grid objects
     * @param t_obs Array of observation times
     * @param band_freq Array of frequency band boundaries
     * @param photons Parameter pack of photon grid objects
     * @details Computes the integrated flux over a frequency band specified by band_freq.
     *          It converts band boundaries to center frequencies, computes the specific flux at each frequency,
     *          and integrates (sums) the flux contributions weighted by the frequency bin widths.
     * @return Array of integrated flux values at each observation time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array flux(Array const& t_obs, Array const& band_freq, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the spectrum at multiple observation times
     * @tparam PhotonGrid Types of photon grid objects
     * @param freqs Array of frequencies
     * @param t_obs Array of observation times
     * @param photons Parameter pack of photon grid objects
     * @details Returns the spectra (as a MeshGrid) for multiple observation times by computing the specific flux
     *          for each frequency and transposing the result to get freq × time format.
     * @return 2D grid of spectra (frequency × time)
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    MeshGrid spectra(Array const& freqs, Array const& t_obs, PhotonGrid& photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the spectrum at a single observation time
     * @tparam PhotonGrid Types of photon grid objects
     * @param freqs Array of frequencies
     * @param t_obs Single observation time
     * @param photons Parameter pack of photon grid objects
     * @details Returns the spectrum (as an Array) at a single observation time by computing the specific flux
     *          for each frequency in the given array.
     * @return Array containing the spectrum at the given time
     * <!-- ************************************************************************************** -->
     */
    template <typename PhotonGrid>
    Array spectrum(Array const& freqs, Real t_obs, PhotonGrid& photons);

    /**
     * @brief Computes resolved sky images at multiple observer times and a single frequency.
     * @tparam PhotonGrid Type of photon grid object
     * @param coord Coordinate grid (needed for phi, theta_view projection)
     * @param shock Shock object (needed for r, theta at interpolated time)
     * @param t_obs Array of observer times (internal units)
     * @param nu_obs Observer frequency (internal units)
     * @param photons Photon grid object
     * @param npixel Number of pixels per side
     * @param pix_size Pixel size in radians (FOV = pix_size × npixel)
     * @return SkyImageResult with [n_frames, npixel, npixel] surface brightness, shared extent and pixel_solid_angle
     */
    template <typename PhotonGrid>
    SkyImageResult sky_image(Coord const& coord, Shock const& shock, Array const& t_obs, Real nu_obs,
                             PhotonGrid& photons, size_t npixel, Real pix_size);

    MeshGrid3d lg2_t;           ///< Log2 of observation time grid
    MeshGrid3d lg2_doppler;     ///< Log2 of Doppler factor grid
    MeshGrid3d lg2_geom_factor; ///< Log2 of observe frame geometric factor (solid angle * r^2 * D^3)
    Real one_plus_z{1};         ///< 1 + redshift
  private:
    Real lumi_dist{1}; ///< Luminosity distance
    // Grid dimensions
    size_t jet_3d{0};             ///< Flag indicating if the jet is non-axis-symmetric (non-zero if true)
    size_t eff_phi_grid{1};       ///< Effective number of phi grid points
    size_t theta_grid{0};         ///< Number of theta grid points
    size_t t_grid{0};             ///< Number of time grid points
    bool jet_spreading_{false};   ///< Jet lateral spreading (from Coord, set by detect_symmetry)
    bool geom_pre_logged_{false}; ///< lg2_geom_factor already holds log2(dOmega r^2) (factorized axisymmetric path)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Builds the time grid and related structures for observation.
     * @details Initializes grid dimensions based on the shock evolution data and sets up arrays for
     *          time, Doppler factor, and emission surface.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * @param luminosity_dist Luminosity distance to the source
     * @param redshift Redshift of the source
     * <!-- ************************************************************************************** -->
     */
    void build_time_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the observation time grid and Doppler factor grid.
     * @details For each grid point, computes the Doppler factor based on the Lorentz factor and
     *          calculates the observed time taking redshift into account.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_t_obs(Coord const& coord, Shock const& shock);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the observe frame solid angle for each grid point.
     * @details Computes the observe frame solid angle as the product of the differential
     *          cosine of theta and either 2π (if the effective phi size is 1) or the differential phi value.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_solid_angle(Coord const& coord, Shock const& shock);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Fused EAT grid + geometric factor computation for non-spreading jets.
     * @details Computes time, lg2_t, lg2_doppler, and lg2_geom_factor in a single pass,
     *          avoiding redundant memory reads of shock.r and lg2_doppler.
     * @param coord Coordinate grid containing angular information
     * @param shock Shock object containing the evolution data
     * <!-- ************************************************************************************** -->
     */
    void calc_eat_non_spreading(Coord const& coord, Shock const& shock);

    /// Converts the linear-valued grids filled by the calc_* methods to log2 space
    /// in vectorized whole-tensor passes (one libm-free pass per tensor).
    void finalize_log_grids();

    /// Fills an existing splat grid (sky positions, luminosities, Gaussian widths) for a single frame.
    template <typename PhotonGrid>
    void build_splat_grid(SplatGrid& grid, Coord const& coord, Shock const& shock, Real t_obs, Real nu_obs,
                          PhotonGrid& photons, std::vector<Real> const& cos_phi_tab,
                          std::vector<Real> const& sin_phi_tab);

    /**
     * <!-- ************************************************************************************** -->
     * @struct InterpState
     * @brief Log-log chord interpolation state for single-time sampling (sky-image splatting only).
     * @details The flux paths interpolate with HermiteSegment; the chord is sufficient here because
     *          each frame samples one target time per cell and the splat kernel's own smoothing
     *          dominates.
     * <!-- ************************************************************************************** -->
     */
    struct InterpState {
        Real slope{0};       ///< Slope for logarithmic interpolation
        Real lg2_L_nu_lo{0}; ///< Lower boundary of specific luminosity (log2 scale)
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Interpolates the specific luminosity in log-log space, returning the result in linear scale.
     * @param state The interpolation state
     * @param lg2_t_obs Observation time (in log2 scale)
     * @param lg2_t_lo Lower boundary of observation time (in log2 scale)
     * @return The interpolated specific luminosity (linear scale)
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] static Real loglog_interpolate(InterpState const& state, Real lg2_t_obs, Real lg2_t_lo) noexcept;

    /**
     * @tparam Photon Photon type of the cell row
     * @param state The interpolation state to update
     * @param row The cell's per-node rows
     * @param k Time grid index
     * @param lg2_nu_src Log2 of the observed frequency
     * @details Computes the log2 luminosity at the two time-grid endpoints of interval k and the
     *          interpolation slope between them.
     * @return True if the slope is finite (both boundaries valid for interpolation), false otherwise
     * <!-- ************************************************************************************** -->
     */
    template <typename Photon>
    bool set_boundaries(InterpState& state, CellRow<Photon> const& row, size_t k, Real lg2_nu_src) noexcept;
};

//========================================================================================================
//                                  template function implementation
//========================================================================================================

/**
 * <!-- ************************************************************************************** -->
 * @brief Helper function that advances an iterator through an array until the value exceeds the target
 * @param value Target value to exceed
 * @param arr Array to iterate through
 * @param it Iterator position (updated by this function)
 * @details Used for efficiently finding the appropriate position in a sorted array without binary search.
 * <!-- ************************************************************************************** -->
 */
inline void iterate_to(Real value, Array const& arr, size_t& it) noexcept {
    while (it < arr.size() && arr(it) < value) {
        it++;
    }
}

/// Same as iterate_to, but also consumes entries equal to the target value.
inline void iterate_through(Real value, Array const& arr, size_t& it) noexcept {
    while (it < arr.size() && arr(it) <= value) {
        it++;
    }
}

template <typename Photon>
bool Observer::set_boundaries(InterpState& state, CellRow<Photon> const& row, size_t k, Real lg2_nu_src) noexcept {
    state.lg2_L_nu_lo = row.node_flux(k, lg2_nu_src);
    const Real lg2_L_nu_hi = row.node_flux(k + 1, lg2_nu_src);
    state.slope = (lg2_L_nu_hi - state.lg2_L_nu_lo) / (row.t[k + 1] - row.t[k]);
    return std::isfinite(state.slope);
}

/**
 * <!-- ************************************************************************************** -->
 * @struct HermiteSegment
 * @brief Local cubic Hermite (log-log, weighted-secant slopes) over one source-time interval.
 * @details Tracks the light-curve bend across spectral-break crossings that a straight chord
 *          cuts; C1 across nodes, exact chord on pure power-law segments and when a neighbor
 *          node is missing (lattice ends, non-finite values).
 * <!-- ************************************************************************************** -->
 */
struct HermiteSegment {
    Real t_lo;  ///< log2 arrival time of the lower node
    Real lo;    ///< log2 flux at the lower node
    Real slope; ///< interpolant slope at the lower node (the chord slope when degenerate)
    Real c2;    ///< quadratic coefficient
    Real c3;    ///< cubic coefficient

    [[nodiscard]] Real operator()(Real t) const noexcept {
        const Real x = t - t_lo;
        return lo + x * (slope + x * (c2 + x * c3));
    }
};

/// Unified "no usable value" marker: absent stencil neighbor, unevaluated window slot,
/// or no previous band. Known-dead nodes are -inf, not NaN, so they are never re-evaluated.
inline constexpr Real missing_node = std::numeric_limits<Real>::quiet_NaN();

/// Node times and reciprocals of one source-time interval, shared by every frequency
/// evaluated on it (the divisions run once per interval instead of once per frequency).
struct SegmentGeometry {
    Real t_lo{0};                                   ///< log2 time of the lower node
    Real dt{0}, dt_prev{0}, dt_next{0};             ///< interval widths: [k, k+1], [k-1, k], [k+1, k+2]
    Real inv_dt{0}, inv_dt_prev{0}, inv_dt_next{0}; ///< their reciprocals
    Real inv_span_prev{0}, inv_span_next{0};        ///< 1 / (dt_prev + dt), 1 / (dt + dt_next)
};

inline SegmentGeometry make_segment_geometry(Real const* t_row, size_t k, size_t t_grid) noexcept {
    SegmentGeometry g;
    g.t_lo = t_row[k];
    g.dt = t_row[k + 1] - g.t_lo;
    g.inv_dt = 1 / g.dt;
    if (k >= 1) {
        g.dt_prev = g.t_lo - t_row[k - 1];
        g.inv_dt_prev = 1 / g.dt_prev;
        g.inv_span_prev = 1 / (g.dt_prev + g.dt);
    }
    if (k + 2 < t_grid) {
        g.dt_next = t_row[k + 2] - t_row[k + 1];
        g.inv_dt_next = 1 / g.dt_next;
        g.inv_span_next = 1 / (g.dt + g.dt_next);
    }
    return g;
}

/// Node slopes are the interval-weighted mean of adjacent secants, then passed through a
/// Fritsch--Carlson monotonicity limiter (sign-matched to the interval chord, magnitude
/// capped at 3|chord|), so a segment can never overshoot its nodes. Unlimited slopes ring
/// around sharp per-cell features (e.g. the Doppler edge as 1/Gamma sweeps a theta ring)
/// and imprint ring-count oscillations on summed light curves. Missing neighbors (NaN)
/// degrade to the chord; non-finite endpoints yield the -inf sentinel, so callers evaluate
/// unconditionally.
inline HermiteSegment make_hermite_segment(SegmentGeometry const& g, Real lo_prev, Real lo, Real hi,
                                           Real hi_next) noexcept {
    const Real chord = (hi - lo) * g.inv_dt;
    if (!std::isfinite(chord)) {
        return {g.t_lo, -con::inf, 0.0, 0.0, 0.0};
    }
    if (!std::isfinite(lo_prev) || !std::isfinite(hi_next)) {
        return {g.t_lo, lo, chord, 0.0, 0.0};
    }
    const Real chord_prev = (lo - lo_prev) * g.inv_dt_prev;
    const Real chord_next = (hi_next - hi) * g.inv_dt_next;

    const auto limit = [chord](Real s) noexcept -> Real {
        if (s * chord <= 0.0) {
            return 0.0; // node is a local extremum of the data: flatten
        }
        const Real cap = 3.0 * chord;
        return std::abs(s) <= std::abs(cap) ? s : cap;
    };
    const Real slope_lo = limit((g.dt_prev * chord + g.dt * chord_prev) * g.inv_span_prev);
    const Real slope_hi = limit((g.dt * chord_next + g.dt_next * chord) * g.inv_span_next);
    const Real c2 = (3 * chord - 2 * slope_lo - slope_hi) * g.inv_dt;
    const Real c3 = (slope_lo + slope_hi - 2 * chord) * g.inv_dt * g.inv_dt;
    return {g.t_lo, lo, slope_lo, c2, c3};
}

/// Rolling three-node flux window keyed by band: after building interval k it holds nodes
/// (k, k+1, k+2), so the next interval costs one fresh evaluation.
struct StencilWindow {
    size_t first_node{0};       ///< node index of vals[0]
    Real nu{missing_node};      ///< band the held fluxes belong to (NaN: nothing held yet)
    std::array<Real, 3> vals{}; ///< fluxes at nodes first_node + {0, 1, 2}; NaN = not evaluated

    /// Window flux at `node` for `band`, or the marker when not held.
    [[nodiscard]] Real held(size_t node, Real band) const noexcept {
        return (nu == band && node - first_node < 3) ? vals[node - first_node] : missing_node;
    }

    /// Window flux at `node` for `band`, evaluated fresh when not held.
    template <typename Photon>
    [[nodiscard]] Real fetch(CellRow<Photon> const& row, size_t node, Real band) const noexcept {
        const Real f = held(node, band);
        return std::isnan(f) ? row.node_flux(node, band) : f;
    }

    /// Builds interval k's segment at band lg2_nu and re-anchors the window to (k, k+1, k+2).
    /// Nodes are fetched only when usable: endpoints always, lo_prev only for a live
    /// interval, hi_next only when the chord upgrades to the cubic.
    template <typename Photon>
    HermiteSegment build(CellRow<Photon> const& row, SegmentGeometry const& seg_geom, size_t k, size_t t_grid,
                         Real lg2_nu) noexcept {
        const Real f_lo = fetch(row, k, lg2_nu);
        const Real f_hi = fetch(row, k + 1, lg2_nu);

        Real f_lo_prev = missing_node;
        Real f_hi_next = missing_node;
        if (std::isfinite(f_lo) && std::isfinite(f_hi)) { // dead intervals stay a 2-node sentinel
            if (k >= 1) {
                f_lo_prev = fetch(row, k - 1, lg2_nu);
            }
            if (std::isfinite(f_lo_prev) && k + 2 < t_grid) { // chord-degenerate: hi_next unused
                f_hi_next = fetch(row, k + 2, lg2_nu);
            }
        }

        first_node = k;
        nu = lg2_nu;
        vals = {f_lo, f_hi, f_hi_next};
        return make_hermite_segment(seg_geom, f_lo_prev, f_lo, f_hi, f_hi_next);
    }
};

/// Lazily-filled node-flux table for the grid path (whole frequency rows per node) --
/// the grid analogue of StencilWindow. Returned rows are views valid until the next
/// begin_cell(); absent or gated-off stencil rows read missing().
template <typename Photon>
struct NodeTable {
    std::span<const Real> lg2_nu_src; ///< requested source-frame frequencies (borrowed)
    MeshGrid vals;                    ///< [t_grid, nu_len] log2 node fluxes
    Array missing_row;                ///< missing_node markers served for absent rows
    std::vector<char> filled;         ///< per-node fill flag for the current cell
    CellRow<Photon> row{};            ///< current cell

    NodeTable(size_t t_grid, Array const& lg2_nu_src)
        : lg2_nu_src(lg2_nu_src.data(), lg2_nu_src.size()),
          vals(MeshGrid::from_shape({t_grid, lg2_nu_src.size()})),
          missing_row({lg2_nu_src.size()}, missing_node),
          filled(t_grid, 0) {}
    NodeTable(size_t, Array&&) = delete; // the borrowed frequencies must outlive the table

    void begin_cell(CellRow<Photon> const& cell) noexcept {
        row = cell;
        std::fill(filled.begin(), filled.end(), char{0});
    }

    /// The node's flux row across all frequencies, evaluated on first use.
    [[nodiscard]] std::span<const Real> node(size_t k) {
        const std::span<Real> dst{&vals(k, 0), lg2_nu_src.size()};
        if (!filled[k]) {
            filled[k] = 1;
            row.node_flux_row(k, lg2_nu_src, dst);
        }
        return dst;
    }

    /// The row served for stencil rows outside the lattice or gated off.
    [[nodiscard]] std::span<const Real> missing() const noexcept { return {&missing_row(0), missing_row.size()}; }

    /// Interval k's stencil rows plus liveness, fetched only when usable: endpoints always;
    /// lo_prev only if some band is live; hi_next only if some live band has finite lo_prev.
    struct Stencil {
        std::span<const Real> lo_prev, lo, hi, hi_next;
        bool any_live; ///< some band has both endpoints finite
    };

    [[nodiscard]] Stencil stencil(size_t k) {
        Stencil st{missing(), node(k), node(k + 1), missing(), false};
        const size_t nu_len = lg2_nu_src.size();
        for (size_t l = 0; l < nu_len && !st.any_live; l++) {
            st.any_live = std::isfinite(st.lo[l]) && std::isfinite(st.hi[l]);
        }
        if (st.any_live && k >= 1) {
            st.lo_prev = node(k - 1);
            bool any_cubic = false; // some band can use the full four-node stencil
            for (size_t l = 0; l < nu_len && !any_cubic; l++) {
                any_cubic = std::isfinite(st.lo_prev[l]) && std::isfinite(st.lo[l]) && std::isfinite(st.hi[l]);
            }
            if (any_cubic && k + 2 < filled.size()) {
                st.hi_next = node(k + 2);
            }
        }
        return st;
    }
};

template <typename PhotonGrid>
MeshGrid Observer::specific_flux(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons) {
    const size_t t_obs_len = t_obs.size();
    const size_t nu_len = nu_obs.size();

    const Array lg2_t_obs = xt::log2(t_obs);
    const Array lg2_nu_src = xt::log2(nu_obs) + fast_log2(one_plus_z);

    MeshGrid F_nu({nu_len, t_obs_len}, 0);

    using Photon = std::remove_reference_t<decltype(photons(0, 0, 0))>;
    NodeTable<Photon> nodes(t_grid, lg2_nu_src);
    MeshGrid col_lg2 = MeshGrid::from_shape({nu_len, t_obs_len});

    for (size_t i = 0; i < eff_phi_grid; i++) {
        const size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            const CellRow row{&lg2_t(i, j, 0), &lg2_doppler(i, j, 0), &lg2_geom_factor(i, j, 0), &photons(eff_i, j, 0)};

            // Skip rows whose arrival times never intersect the requested window.
            if (row.t[t_grid - 1] < lg2_t_obs(0) || row.t[0] > lg2_t_obs(t_obs_len - 1)) {
                continue;
            }

            nodes.begin_cell(row);

            size_t t_idx = 0;
            iterate_to(row.t[0], lg2_t_obs, t_idx);
            const size_t col_first = t_idx;

            // Walk the intervals front to back; each one processes the block of
            // observation times it brackets.
            for (size_t k = 0; t_idx < t_obs_len && k < t_grid - 1; k++) {
                if (row.t[k + 1] < lg2_t_obs(t_idx)) {
                    continue; // no observation time falls inside this interval
                }
                const size_t idx_start = t_idx;
                iterate_to(row.t[k + 1], lg2_t_obs, t_idx);

                const auto stencil = nodes.stencil(k);
                if (!stencil.any_live) { // every band dead: sentinel segments emit -inf
                    xt::view(col_lg2, xt::all(), xt::range(idx_start, t_idx)) = -con::inf;
                    continue;
                }

                // One cubic segment per frequency, evaluated over the bracketed times.
                const SegmentGeometry seg_geom = make_segment_geometry(row.t, k, t_grid);
                for (size_t l = 0; l < nu_len; l++) {
                    const HermiteSegment seg = make_hermite_segment(seg_geom, stencil.lo_prev[l], stencil.lo[l],
                                                                    stencil.hi[l], stencil.hi_next[l]);
                    for (size_t idx = idx_start; idx < t_idx; idx++) {
                        col_lg2(l, idx) = seg(lg2_t_obs(idx));
                    }
                }
            }

            if (col_first < t_idx) {
                xt::noalias(xt::view(F_nu, xt::all(), xt::range(col_first, t_idx))) +=
                    xt::exp2(xt::view(col_lg2, xt::all(), xt::range(col_first, t_idx)));
            }
        }
    }

    F_nu *= one_plus_z / (lumi_dist * lumi_dist);

    return F_nu;
}

template <typename PhotonGrid>
Array Observer::specific_flux_series(Array const& t_obs, Array const& nu_obs, PhotonGrid& photons) {
    const size_t t_obs_len = t_obs.size();

    if (nu_obs.size() != t_obs_len) {
        throw std::invalid_argument("specific_flux_series requires nu_obs and t_obs to have the same length");
    }

    const Array lg2_t_obs = xt::log2(t_obs);
    const Array lg2_nu_src = xt::log2(nu_obs) + fast_log2(one_plus_z);

    Array F_nu = xt::zeros<Real>({t_obs_len});
    // Per-column log2 contributions; exp2 runs once per column. Never reset between cells:
    // every flushed column is written by the current cell (dead bands write -inf).
    Array col_lg2 = Array::from_shape({t_obs_len});

    for (size_t i = 0; i < eff_phi_grid; i++) {
        const size_t eff_i = i * jet_3d;
        for (size_t j = 0; j < theta_grid; j++) {
            const CellRow row{&lg2_t(i, j, 0), &lg2_doppler(i, j, 0), &lg2_geom_factor(i, j, 0), &photons(eff_i, j, 0)};

            size_t idx = 0;
            iterate_to(row.t[0], lg2_t_obs, idx);
            const size_t col_first = idx;

            StencilWindow window;

            // Walk the intervals front to back; each one processes the block of
            // observation points it brackets.
            for (size_t k = 0; idx < t_obs_len && k < t_grid - 1; k++) {
                if (row.t[k + 1] < lg2_t_obs(idx)) {
                    continue; // no observation point falls inside this interval
                }
                const size_t block_first = idx;
                iterate_through(row.t[k + 1], lg2_t_obs, idx);

                const SegmentGeometry seg_geom = make_segment_geometry(row.t, k, t_grid);
                Real prev_nu = missing_node;
                HermiteSegment seg{};
                for (size_t s = block_first; s < idx; s++) {
                    const Real lg2_nu = lg2_nu_src(s);
                    if (lg2_nu != prev_nu) { // same-band runs share the whole segment
                        seg = window.build(row, seg_geom, k, t_grid, lg2_nu);
                        prev_nu = lg2_nu;
                    }
                    col_lg2(s) = seg(lg2_t_obs(s));
                }
            }

            if (col_first < idx) {
                xt::noalias(xt::view(F_nu, xt::range(col_first, idx))) +=
                    xt::exp2(xt::view(col_lg2, xt::range(col_first, idx)));
            }
        }
    }

    F_nu *= one_plus_z / (lumi_dist * lumi_dist);

    return F_nu;
}

template <typename PhotonGrid>
Array Observer::specific_flux(Array const& t_obs, Real nu_obs, PhotonGrid& photons) {
    return xt::view(specific_flux(t_obs, Array({nu_obs}), photons), 0);
}

template <typename PhotonGrid>
Array Observer::spectrum(Array const& freqs, Real t_obs, PhotonGrid& photons) {
    return xt::view(spectra(freqs, Array({t_obs}), photons), 0);
}

template <typename PhotonGrid>
MeshGrid Observer::spectra(Array const& freqs, Array const& t_obs, PhotonGrid& photons) {
    return xt::transpose(specific_flux(t_obs, freqs, photons));
}

template <typename PhotonGrid>
Array Observer::flux(Array const& t_obs, Array const& band_freq, PhotonGrid& photons) {
    MeshGrid F_nu = specific_flux(t_obs, band_freq, photons);
    Array weights;
    compute_boole_weights(band_freq, weights);
    Array flux({t_obs.size()}, 0);
    for (size_t i = 0; i < band_freq.size(); ++i) {
        for (size_t j = 0; j < flux.size(); ++j) {
            flux(j) += F_nu(i, j) * weights(i);
        }
    }
    return flux;
}

//========================================================================================================
//                                  sky_image template implementation
//========================================================================================================

template <typename PhotonGrid>
void Observer::build_splat_grid(SplatGrid& grid, Coord const& coord, Shock const& shock, Real t_obs, Real nu_obs,
                                PhotonGrid& photons, std::vector<Real> const& cos_phi_tab,
                                std::vector<Real> const& sin_phi_tab) {
    const Real lg2_t_target = fast_log2(t_obs);
    const Real lg2_nu_src = fast_log2(nu_obs * one_plus_z);
    const Real inv_d_A = one_plus_z * one_plus_z / lumi_dist;
    const Real cos_obs = std::cos(coord.theta_view);
    const Real sin_obs = std::sin(coord.theta_view);
    const bool is_axisym = (eff_phi_grid == 1);
    const size_t copies_per_phys = grid.n_phi / eff_phi_grid;

    for (size_t i_phys = 0; i_phys < eff_phi_grid; ++i_phys) {
        const size_t eff_i = i_phys * jet_3d;
        const Real cp_fixed = is_axisym ? 0 : std::cos(coord.phi[i_phys] - coord.phi_view);
        const Real sp_fixed = is_axisym ? 0 : std::sin(coord.phi[i_phys] - coord.phi_view);

        for (size_t j = 0; j < theta_grid; ++j) {
            const CellRow row{&lg2_t(i_phys, j, 0), &lg2_doppler(i_phys, j, 0), &lg2_geom_factor(i_phys, j, 0),
                              &photons(eff_i, j, 0)};
            if (lg2_t_target < row.t[0] || lg2_t_target >= row.t[t_grid - 1]) {
                continue;
            }

            // Find time bracket
            size_t k = 0;
            while (k + 2 < t_grid && row.t[k + 1] < lg2_t_target) {
                ++k;
            }

            // Interpolate luminosity
            InterpState state_L;
            if (!set_boundaries(state_L, row, k, lg2_nu_src)) {
                continue;
            }
            const Real L_val = loglog_interpolate(state_L, lg2_t_target, row.t[k]);
            if (!std::isfinite(L_val) || L_val <= 0) {
                continue;
            }

            // Interpolate position
            const Real w = (lg2_t_target - lg2_t(i_phys, j, k)) / (lg2_t(i_phys, j, k + 1) - lg2_t(i_phys, j, k));
            const Real r_ang = (shock.r(eff_i, j, k) + w * (shock.r(eff_i, j, k + 1) - shock.r(eff_i, j, k))) * inv_d_A;
            const Real theta_val = jet_spreading_ ? shock.theta(eff_i, j, k) +
                                                        w * (shock.theta(eff_i, j, k + 1) - shock.theta(eff_i, j, k))
                                                  : shock.theta(eff_i, j, 0);
            const Real ct = std::cos(theta_val);
            const Real st = std::sin(theta_val);
            const Real L_per_copy = L_val / copies_per_phys;

            // Project onto sky plane for each azimuthal copy
            for (size_t c = 0; c < copies_per_phys; ++c) {
                const Real cp = is_axisym ? cos_phi_tab[c] : cp_fixed;
                const Real sp = is_axisym ? sin_phi_tab[c] : sp_fixed;
                const Real x_ang = r_ang * (st * cp * cos_obs - ct * sin_obs);
                const Real y_ang = r_ang * st * sp;

                const size_t i_splat = i_phys * copies_per_phys + c;
                grid.x(i_splat, j) = x_ang;
                grid.y(i_splat, j) = y_ang;
                grid.L(i_splat, j) = L_per_copy;
                grid.xmin = std::min(grid.xmin, x_ang);
                grid.xmax = std::max(grid.xmax, x_ang);
                grid.ymin = std::min(grid.ymin, y_ang);
                grid.ymax = std::max(grid.ymax, y_ang);
            }
        }
    }

    // Convert luminosity to surface brightness
    grid.L *= one_plus_z / (lumi_dist * lumi_dist);
    grid.compute_sigma();
}

template <typename PhotonGrid>
SkyImageResult Observer::sky_image(Coord const& coord, Shock const& shock, Array const& t_obs, Real nu_obs,
                                   PhotonGrid& photons, size_t npixel, Real pix_size) {
    const size_t n_frames = t_obs.size();
    SkyImageResult result;
    result.image = MeshGrid3d({n_frames, npixel, npixel}, 0);

    if (npixel == 0 || t_grid < 2) {
        result.extent = {0, 0, 0, 0};
        return result;
    }

    // Compute extent and pixel_solid_angle once (depends only on pix_size and npixel)
    const Real half_fov = 0.5 * pix_size * static_cast<Real>(npixel);
    result.extent = {-half_fov, half_fov, -half_fov, half_fov};
    result.pixel_solid_angle = pix_size * pix_size;

    // 2D Gaussian splatting
    const size_t n_splat_phi = (eff_phi_grid == 1) ? 4 * theta_grid : eff_phi_grid;
    SplatGrid grid(n_splat_phi, theta_grid);

    const size_t copies_per_phys = n_splat_phi / eff_phi_grid;
    std::vector<Real> cos_phi_tab, sin_phi_tab;
    if (eff_phi_grid == 1) {
        cos_phi_tab.resize(copies_per_phys);
        sin_phi_tab.resize(copies_per_phys);
        for (size_t c = 0; c < copies_per_phys; ++c) {
            const Real phi = 2.0 * con::pi * c / n_splat_phi;
            cos_phi_tab[c] = std::cos(phi);
            sin_phi_tab[c] = std::sin(phi);
        }
    }

    const bool is_axisym = (jet_3d == 0);
    const bool azimuthal_avg = is_axisym && (coord.theta_view == 0);
    const bool mirror_y = is_axisym && !azimuthal_avg;

    for (size_t f = 0; f < n_frames; ++f) {
        grid.reset();
        build_splat_grid(grid, coord, shock, t_obs(f), nu_obs, photons, cos_phi_tab, sin_phi_tab);
        grid.render(result.image, f, pix_size, azimuthal_avg, mirror_y);
    }
    return result;
}
