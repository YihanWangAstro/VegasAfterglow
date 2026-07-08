//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <limits>

#include "../core/grid-refinement.h"
#include "../core/physics.h"
#include "../core/quadrature.h"
#include "../dynamics/shock.h"
#include "../util/macros.h"
#include "../util/utilities.h"
#include "syn-concepts.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct InverseComptonY
 * @brief Handles Inverse Compton Y parameter calculations and related threshold values.
 * <!-- ************************************************************************************** -->
 */
struct InverseComptonY {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes an InverseComptonY object with frequency thresholds, magnetic field and Y parameter.
     * @details Computes characteristic gamma values and corresponding frequencies, then determines cooling regime.
     * @param gamma_m Characteristic Lorentz factor for the minimum Lorentz factor
     * @param gamma_c Characteristic Lorentz factor for the cooling Lorentz factor
     * @param p Spectral index of electron distribution
     * @param B Magnetic field strength
     * @param Y_T Thomson scattering Y parameter
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY(Real gamma_m, Real gamma_c, Real p, Real B, Real Y_T, bool is_KN) noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Default constructor that initializes all member variables to zero.
     * <!-- ************************************************************************************** -->
     */
    InverseComptonY() noexcept;

    // Member variables
    Real gamma_m_hat{1}; ///< Lorentz factor threshold for minimum energy electrons
    Real gamma_c_hat{1}; ///< Lorentz factor threshold for cooling electrons
    Real gamma_self{1};  ///< Lorentz factor threshold for gamma_self = gamma_self_hat
    Real gamma0{1};      ///< Lorentz factor threshold for Y(gamma0) = 1x
    Real Y_T{0};         ///< Thomson scattering Y parameter
    size_t regime{0};    ///< Indicator for the operating regime (1=fast IC cooling, 2=slow IC cooling, 3=special case)

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter for a given electron Lorentz factor.
     * @details Different scaling relations apply depending on the cooling regime and gamma value.
     * @param gamma Electron Lorentz factor
     * @return The effective Y parameter at the given gamma
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real gamma_spectrum(Real gamma) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the effective Y parameter at the electron Lorentz factor corresponding to a given synchrotron frequency.
     * @details Converts frequency to electron Lorentz factor via the synchrotron relation, then evaluates Y(gamma).
     * @param nu Synchrotron frequency
     * @return The effective Y parameter at the corresponding electron Lorentz factor
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Real nu_spectrum(Real nu) const noexcept;

    void update_cooling_breaks(Real gamma_c, Real Y_T) noexcept;

  private:
    Real gamma_m_{1};    ///< Characteristic Lorentz factor for the minimum Lorentz factor
    Real B_{0};          ///< Magnetic field strength
    Real p_{2.3};        ///< Spectral index of electron distribution
    Real gamma_self3{1}; ///< Precomputed gamma_self^3 for efficiency
    void update_gamma0(Real gamma_c) noexcept;

    Real compute_gamma_hat(Real gamma) const noexcept;

    BrokenPowerLaw<5> segments_;

    void build_segments() noexcept;
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast KN correction lookup: returns sigma_KN/sigmaT at frequency nu.
 * @param nu Photon frequency.
 * @return Dimensionless correction factor to Thomson cross-section.
 * <!-- ************************************************************************************** -->
 */
Real compton_correction(Real nu);

/**
 * <!-- ************************************************************************************** -->
 * @brief Fast KN correction lookup returning both sigma_KN/sigmaT and its base-2 log.
 * @details The log2 value lets callers form local power-law slopes of KN-corrected
 *          spectra without a log call per evaluation.
 * <!-- ************************************************************************************** -->
 */
void compton_correction_pair(Real nu, Real& corr, Real& lg2_corr);

/**
 * <!-- ************************************************************************************** -->
 * @struct ICPhoton
 * @tparam Electrons Type of the electron distribution
 * @tparam Photons Type of the photon distribution
 * @brief Represents a single inverse Compton (IC) photon.
 * @details Contains methods to compute the photon intensity I_nu and to generate an IC photon spectrum based
 *          on electron and seed photon properties.
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
struct ICPhoton {
  public:
    static_assert(SynElectronModel<Electrons>,
                  "ICPhoton: Electrons must satisfy the SynElectronModel concept (see syn-concepts.h)");
    static_assert(SynPhotonModel<Photons>,
                  "ICPhoton: Photons must satisfy the SynPhotonModel concept (see syn-concepts.h)");

    /// Default constructor
    ICPhoton() = default;

    ICPhoton(Electrons const& electrons, Photons const& photons, bool KN, Real nu_eval_min = 0,
             Real nu_eval_max = std::numeric_limits<Real>::infinity()) noexcept;
    /**
     * <!-- ************************************************************************************** -->
     * @brief Returns the photon-specific intensity.
     * @param nu The frequency at which to compute the specific intensity
     * @return The specific intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_I_nu(Real nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the base-2 logarithm of the photon-specific intensity at a given frequency.
     * @param log2_nu The base-2 logarithm of the frequency
     * @return The base-2 logarithm of the photon specific intensity at the given frequency
     * <!-- ************************************************************************************** -->
     */
    Real compute_log2_I_nu(Real log2_nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Generates the IC spectrum grid now instead of lazily on the first query.
     * <!-- ************************************************************************************** -->
     */
    void build() {
        if (!generated) {
            generate_spectrum();
        }
    }

    Photons photons;

    Electrons electrons;

  private:
    void generate_spectrum();

    /**
     * <!-- ************************************************************************************** -->
     * @brief Grid parameter structure holding ranges and sizes for gamma and nu
     * <!-- ************************************************************************************** -->
     */
    struct GridParams {
        Real gamma_min;
        Real gamma_max;

        Real nu_min;
        Real nu_max;
        Real nu_IC_min;
        Real nu_IC_max;
    };

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes grid parameters based on electron and photon distributions
     * @return GridParams structure containing all computed parameters
     * <!-- ************************************************************************************** -->
     */
    GridParams compute_grid_params() const;

    void initialize_grids(GridParams const& params, Array& gamma, Array& nu_seed, Array& lg2_nu_seed, Array& dnu_seed);

    void sample_distributions(Array const& gamma, Array const& lg2_nu_seed, Array& dN_e_boost, Array& I_nu_seed);

    void compute_IC_spectrum(Array const& gamma, Array const& dN_e_boost, Array const& nu_seed,
                             Array const& lg2_nu_seed, Array const& dnu_seed, Array const& I_nu_seed);

    static void build_cdf_thomson(Array const& fv_th, Array const& lg2fv_th, Array const& nu_seed,
                                  Array const& dnu_seed, Array const& lg2r, Array const& inv_lg2r, Array& cdf_buf,
                                  Array& ratio_buf);

    static void build_cdf_KN(size_t i_gamma, Real gamma_i, Array const& corr_lat, Array const& lg2corr_lat,
                             Array const& fv_th, Array const& lg2fv_th, Array const& nu_seed, Array const& dnu_seed,
                             Array const& lg2r, Array const& inv_lg2r, Array const& cdf_th, Array const& ratio_th,
                             Array& fv_buf, Array& cdf_buf, Array& ratio_buf);

    static void accumulate_IC(Real dN_e_boost, long n_off, Array const& nu_seed, Array const& dnu_seed,
                              Array const& fv_buf, Array const& cdf_buf, Array const& ratio_buf, Array& I_buf,
                              size_t spec_size);

    /// The gamma and seed grids are commensurate lattices over a fundamental
    /// log2 quantum (1/20 decade): gamma advances gamma_mult quanta per node
    /// (5 points per decade), the seed nu_mult quanta (4 per decade). Every
    /// product gamma_i nu_j then sits on the quantum lattice at index
    /// gamma_mult*i + nu_mult*j, so the KN correction is evaluated once per
    /// lattice node instead of once per (i, j) pair.
    static constexpr Real lattice_quantum = 3.321928094887362 / 8; // log2(10)/8
    static constexpr size_t gamma_mult = 2;                        // 4 points per decade
    static constexpr size_t nu_mult = 2;                           // 4 points per decade
    static constexpr size_t ic_mult = 2;                           // 4 points per decade (output grid)

    /// Quantum index of the first nu_IC node in x-space (phase-locked to the
    /// seed lattice: for gamma node i and output node k the seed-lattice
    /// position of x = nu_IC(k)/(4 x0 gamma_i^2) is ic_idx0_ + ic_mult k - 2 gamma_mult i).
    long ic_idx0_{0};

    Array log2_nu_IC;

    Array log2_I_nu_IC;

    Array interp_slope;

    bool KN{false}; // Klein-Nishina flag

    Real nu_eval_min_{0}; ///< Lowest frequency the observer will sample (comoving frame)

    Real nu_eval_max_{std::numeric_limits<Real>::infinity()}; ///< Highest frequency the observer will sample

    mutable Real log2_nu_theory_max_{con::inf}; ///< log2 of the un-clamped upper bound (band-contract tripwire)

    mutable Real log2_nu_theory_min_{-con::inf}; ///< log2 of the un-clamped lower bound (band-contract tripwire)

    bool generated{false};
};

/// Defines a 3D grid (using xt::xtensor) for storing ICPhoton objects.
template <typename Electrons, typename Photons>
using ICPhotonGrid = xt::xtensor<ICPhoton<Electrons, Photons>, 3>;

template <typename Electrons>
using ElectronGrid = xt::xtensor<Electrons, 3>;

template <typename Photons>
using PhotonGrid = xt::xtensor<Photons, 3>;

inline constexpr Real IC_x0 = 0.47140452079103166;
inline constexpr Real inv_4x0 = 1 / (4 * IC_x0);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

template <typename Electrons, typename Photons>
ICPhoton<Electrons, Photons>::ICPhoton(Electrons const& electrons, Photons const& photons, bool KN, Real nu_eval_min,
                                       Real nu_eval_max) noexcept
    : photons(photons), electrons(electrons), KN(KN), nu_eval_min_(nu_eval_min), nu_eval_max_(nu_eval_max) {}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::generate_spectrum() {
    GridParams params = compute_grid_params();

    auto positive_finite = [](Real x) { return std::isfinite(x) && x > 0; };
    if (!(positive_finite(params.gamma_min) && positive_finite(params.gamma_max) && positive_finite(params.nu_min) &&
          positive_finite(params.nu_max) && positive_finite(params.nu_IC_min) && positive_finite(params.nu_IC_max))) {
        generated = true;
        return;
    }

    // Reused across cells: only log2_nu_IC/log2_I_nu_IC/interp_slope persist per object.
    static thread_local Array gamma, nu_seed, lg2_nu_seed, dnu_seed, dN_e_boost, I_nu_seed;
    initialize_grids(params, gamma, nu_seed, lg2_nu_seed, dnu_seed);

    sample_distributions(gamma, lg2_nu_seed, dN_e_boost, I_nu_seed);
    compute_IC_spectrum(gamma, dN_e_boost, nu_seed, lg2_nu_seed, dnu_seed, I_nu_seed);

    generated = true;
}

template <typename Electrons, typename Photons>
typename ICPhoton<Electrons, Photons>::GridParams ICPhoton<Electrons, Photons>::compute_grid_params() const {
    GridParams params;

    constexpr Real eps_tail = 1e-2;     // target suppression level for single-cutoff tail
    constexpr Real safety_margin = 2.0; // retain headroom beyond estimated tail
    constexpr Real min_tail_factor = 5.0;

    const Real tail_factor = std::max(-std::log(eps_tail), min_tail_factor);

    // A soft tail below the lowest electron break; not node-aligned.
    params.gamma_min = std::min(electrons.gamma_m, electrons.gamma_c) / 30;
    params.gamma_max = std::max(electrons.gamma_M * tail_factor, params.gamma_min);

    // A whole-decade margin lands the seed peak min(nu_a, nu_m) exactly on a grid
    // node in every cell (the seed lattice steps a quarter-decade), so its in-bin
    // position cannot drift cell-to-cell — that drift is the dominant source of
    // quadrature noise at the low-frequency end.
    params.nu_min = std::min(photons.nu_a, photons.nu_m) / 10;
    params.nu_max = std::max(photons.nu_M * tail_factor, params.nu_min);

    params.nu_IC_min = 4 * IC_x0 * params.nu_min * params.gamma_min * params.gamma_min;

    const Real nu_ic_base = 4 * IC_x0 * photons.nu_M * electrons.gamma_M * electrons.gamma_M;
    const Real nu_ic_single_cut = std::max(nu_ic_base * tail_factor * tail_factor, nu_ic_base * tail_factor);
    params.nu_IC_max = nu_ic_single_cut * safety_margin;

    // Clamp the output grid to the band the observer will actually sample
    // (with margin). The theoretical range spans tens of decades while the
    // observer only ever evaluates a few; points outside the band are never
    // read. Keep a minimal slice so the grid stays well-formed when the band
    // falls outside the theoretical range.
    constexpr Real band_margin = 4.0;
    constexpr Real min_span = 16.0;
    log2_nu_theory_max_ = fast_log2(params.nu_IC_max);
    log2_nu_theory_min_ = fast_log2(params.nu_IC_min);
    params.nu_IC_min = std::max(params.nu_IC_min, std::min(nu_eval_min_ / band_margin, params.nu_IC_max / min_span));
    params.nu_IC_max = std::min(params.nu_IC_max, std::max(nu_eval_max_ * band_margin, params.nu_IC_min * min_span));

    return params;
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::initialize_grids(GridParams const& params, Array& gamma, Array& nu_seed,
                                                    Array& lg2_nu_seed, Array& dnu_seed) {

    // Seed and gamma grids share the exact lattice_step so gamma-nu products of
    // the two lattices land on one 1D log2 lattice (see build_cdf_KN).
    auto build_lattice = [](Real lg2_lo, Real lg2_hi, Real step, Array& grid, Array& lg2_grid) {
        const size_t n = std::max<size_t>(2, static_cast<size_t>(std::ceil((lg2_hi - lg2_lo) / step)) + 1);
        grid.resize({n});
        lg2_grid.resize({n});
        for (size_t i = 0; i < n; ++i) {
            lg2_grid(i) = lg2_lo + step * static_cast<Real>(i);
            grid(i) = fast_exp2(lg2_grid(i));
        }
    };
    static thread_local Array lg2_gamma_buf;
    build_lattice(fast_log2(params.nu_min), fast_log2(params.nu_max), nu_mult * lattice_quantum, nu_seed, lg2_nu_seed);
    compute_bin_widths(nu_seed, dnu_seed);
    build_lattice(fast_log2(params.gamma_min), fast_log2(params.gamma_max), gamma_mult * lattice_quantum, gamma,
                  lg2_gamma_buf);

    // Output grid on the same quantum lattice, phase-locked per cell so every
    // x = nu_IC(k)/(4 x0 gamma_i^2) lands exactly on a seed-lattice quantum:
    // accumulate_IC then uses pure integer indexing (no scans, no float
    // compares). The phase drifts continuously with the cell's break
    // frequencies, so interior nodes move smoothly (no discrete jumps).
    const Real q = lattice_quantum;
    const Real phase = lg2_nu_seed(0) + 2 * lg2_gamma_buf(0) + fast_log2(4 * IC_x0);
    const long n_lo = static_cast<long>(std::floor((fast_log2(params.nu_IC_min) - phase) / (q * ic_mult)));
    const long n_hi = static_cast<long>(std::ceil((fast_log2(params.nu_IC_max) - phase) / (q * ic_mult)));
    const size_t n_ic = static_cast<size_t>(std::max<long>(n_hi - n_lo, 1)) + 1;
    ic_idx0_ = n_lo * static_cast<long>(ic_mult);
    log2_nu_IC.resize({n_ic});
    for (size_t k = 0; k < n_ic; ++k) {
        log2_nu_IC(k) = phase + q * static_cast<Real>(ic_idx0_ + static_cast<long>(ic_mult * k));
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::sample_distributions(Array const& gamma, Array const& lg2_nu_seed, Array& dN_e_boost,
                                                        Array& I_nu_seed) {
    const size_t g_size = gamma.size();
    dN_e_boost.resize({g_size});

    // Trapezoidal weights on the log-spaced gamma grid — deliberately NOT a
    // higher-order or power-law-exact scheme. The trapezoid's mass excess on
    // steep segments cancels against the piecewise-linear kernel deficit at
    // the emission edge; exact bin masses break that cancellation and degrade
    // the highest frequencies (where SSC matters most) far more than they
    // help the low-frequency end. At this grid density no nodal scheme
    // resolves the edge ridge, so the cancellation is the better trade.
    for (size_t i = 0; i < g_size; ++i) {
        const Real gamma_i = gamma(i);
        const Real dgamma_i = 0.5 * ((i + 1 < g_size ? gamma(i + 1) : gamma(i)) - (i > 0 ? gamma(i - 1) : gamma(i)));
        dN_e_boost(i) = electrons.compute_column_den(gamma_i) / (gamma_i * gamma_i) * dgamma_i;
    }

    const size_t nu_size = lg2_nu_seed.size();
    I_nu_seed.resize({nu_size});
    for (size_t j = 0; j < nu_size; ++j) {
        I_nu_seed(j) = fast_exp2(photons.compute_log2_I_nu(lg2_nu_seed(j)));
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Integrates one seed-spectrum bin as the local power law through its endpoints.
 * @details The seed spectrum is piecewise power-law, so the closed form
 *          (f_hi nu_hi - f_lo nu_lo)/(s+1) is exact per segment. A linear trapezoid on the
 *          coarse log-spaced seed grid overestimates every bin of a steep segment and the
 *          bias accumulates through the CDF (tens of percent at the Compton peak and above).
 *          Falls back to the trapezoid when an endpoint is zero.
 * <!-- ************************************************************************************** -->
 */
inline Real power_law_bin_integral(Real f_lo, Real f_hi, Real nu_lo, Real nu_hi, Real lg2f_lo, Real lg2f_hi, Real lg2r,
                                   Real inv_lg2r, Real trap) {
    if (!(f_lo > 0) || !(f_hi > 0)) {
        return trap;
    }
    const Real s1 = 1 + (lg2f_hi - lg2f_lo) * inv_lg2r; // local slope + 1
    if (std::abs(s1) > 1e-3) {
        return (f_hi * nu_hi - f_lo * nu_lo) / s1;
    }
    constexpr Real ln2 = 0.6931471805599453;
    return f_lo * nu_lo * lg2r * ln2; // s -> -1: integral of f nu^-1 = f_lo nu_lo ln(nu_hi/nu_lo)
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::build_cdf_thomson(Array const& fv_th, Array const& lg2fv_th, Array const& nu_seed,
                                                     Array const& dnu_seed, Array const& lg2r, Array const& inv_lg2r,
                                                     Array& cdf_buf, Array& ratio_buf) {
    const int nu_last = static_cast<int>(nu_seed.size()) - 1;

    cdf_buf(nu_last) = 0;
    for (int j = nu_last - 1; j >= 0; --j) {
        const Real trap = 0.5 * (fv_th(j) + fv_th(j + 1)) * dnu_seed(j);
        const Real exact = power_law_bin_integral(fv_th(j), fv_th(j + 1), nu_seed(j), nu_seed(j + 1), lg2fv_th(j),
                                                  lg2fv_th(j + 1), lg2r(j), inv_lg2r(j), trap);
        cdf_buf(j) = cdf_buf(j + 1) + exact;
        ratio_buf(j) = (trap > 0) ? exact / trap : 1;
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::build_cdf_KN(size_t i_gamma, Real gamma_i, Array const& corr_lat,
                                                Array const& lg2corr_lat, Array const& fv_th, Array const& lg2fv_th,
                                                Array const& nu_seed, Array const& dnu_seed, Array const& lg2r,
                                                Array const& inv_lg2r, Array const& cdf_th, Array const& ratio_th,
                                                Array& fv_buf, Array& cdf_buf, Array& ratio_buf) {
    const int nu_last = static_cast<int>(nu_seed.size()) - 1;

    // Below x = h gamma nu / (me c^2) = 1e-4 the KN correction is within 2e-5
    // of unity, so those bins are the precomputed Thomson ones shifted by a
    // constant: only the tail above nu_split needs per-gamma corrections.
    constexpr Real x_thomson = 1e-4;
    const Real nu_split = x_thomson * (con::me * con::c2 / con::h) / gamma_i;

    int j_split = 0;
    while (j_split < nu_last && nu_seed(j_split) < nu_split) {
        ++j_split;
    }

    // gamma_i nu_j sits at lattice index i_gamma + j (shared log2 grid step).
    Real corr = corr_lat(i_gamma + nu_mult * static_cast<size_t>(nu_last));
    Real lg2_corr = lg2corr_lat(i_gamma + nu_mult * static_cast<size_t>(nu_last));
    fv_buf(nu_last) = fv_th(nu_last) * corr;
    Real lg2f_hi = lg2fv_th(nu_last) + lg2_corr;
    cdf_buf(nu_last) = 0;

    for (int j = nu_last - 1; j >= j_split; --j) {
        corr = corr_lat(i_gamma + nu_mult * static_cast<size_t>(j));
        lg2_corr = lg2corr_lat(i_gamma + nu_mult * static_cast<size_t>(j));
        fv_buf(j) = fv_th(j) * corr;
        const Real lg2f_lo = lg2fv_th(j) + lg2_corr;

        const Real trap = 0.5 * (fv_buf(j) + fv_buf(j + 1)) * dnu_seed(j);
        const Real exact = power_law_bin_integral(fv_buf(j), fv_buf(j + 1), nu_seed(j), nu_seed(j + 1), lg2f_lo,
                                                  lg2f_hi, lg2r(j), inv_lg2r(j), trap);
        cdf_buf(j) = cdf_buf(j + 1) + exact;
        ratio_buf(j) = (trap > 0) ? exact / trap : 1;
        lg2f_hi = lg2f_lo;
    }

    if (j_split > 0) {
        const Real delta = cdf_buf(j_split) - cdf_th(j_split);
        for (int j = j_split - 1; j >= 0; --j) {
            fv_buf(j) = fv_th(j);
            ratio_buf(j) = ratio_th(j);
            cdf_buf(j) = cdf_th(j) + delta;
        }
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::accumulate_IC(Real dN_e_boost, long n_off, Array const& nu_seed,
                                                 Array const& dnu_seed, Array const& fv_buf, Array const& cdf_buf,
                                                 Array const& ratio_buf, Array& I_buf, size_t spec_size) {
    const long ns_top = (static_cast<long>(nu_seed.size()) - 1) * static_cast<long>(nu_mult);

    if (cdf_buf(0) <= 0)
        return;

    // x(k) = nu_IC(k)/(4 x0 gamma_i^2) sits at seed-lattice quantum index
    // n = n_off + ic_mult k: the containing bin is n / nu_mult and the in-bin
    // offset n % nu_mult takes one of nu_mult discrete values, so the partial
    // bin uses the same linear formula as before with x = nu_lo * 2^{fi q}.
    static const std::array<Real, nu_mult> expq = [] {
        std::array<Real, nu_mult> e{};
        for (size_t f = 0; f < nu_mult; ++f) {
            e[f] = fast_exp2(lattice_quantum * static_cast<Real>(f));
        }
        return e;
    }();

    // Plateau: x below the seed range, all seed photons scatter (CDF = cdf_buf[0]).
    const Real plateau = dN_e_boost * cdf_buf(0);
    size_t k = 0;
    long n = n_off;
    for (; k < spec_size && n < 0; ++k, n += static_cast<long>(ic_mult)) {
        I_buf(k) += plateau;
    }

    for (; k < spec_size && n < ns_top; ++k, n += static_cast<long>(ic_mult)) {
        const size_t j = static_cast<size_t>(n) / nu_mult;
        const size_t fi = static_cast<size_t>(n) % nu_mult;
        const Real nu_lo = nu_seed(j);
        const Real dnu = dnu_seed(j);
        const Real f_lo = fv_buf(j);
        const Real f_hi = fv_buf(j + 1);
        const Real frac = (nu_lo * expq[fi] - nu_lo) / dnu;
        const Real rem = 1.0 - frac;
        const Real f_seed = f_lo * rem + f_hi * frac;
        I_buf(k) += dN_e_boost * (cdf_buf(j + 1) + 0.5 * (f_seed + f_hi) * rem * dnu * ratio_buf(j));
    }
}

template <typename Electrons, typename Photons>
void ICPhoton<Electrons, Photons>::compute_IC_spectrum(Array const& gamma, Array const& dN_e_boost,
                                                       Array const& nu_seed, Array const& lg2_nu_seed,
                                                       Array const& dnu_seed, Array const& I_nu_seed) {
    const size_t spec_size = log2_nu_IC.size();
    const size_t nu_size = nu_seed.size();
    const int gamma_size = static_cast<int>(std::min(gamma.size(), dN_e_boost.size()));

    log2_I_nu_IC = Array::from_shape({spec_size});

    static thread_local Array I_buf, cdf_buf, fv_buf, ratio_buf, fv_th, lg2fv_th, lg2r, inv_lg2r;
    I_buf.resize({spec_size});
    cdf_buf.resize({nu_size});
    fv_buf.resize({nu_size});
    ratio_buf.resize({nu_size});
    fv_th.resize({nu_size});
    lg2fv_th.resize({nu_size});
    lg2r.resize({nu_size});
    inv_lg2r.resize({nu_size});
    I_buf.fill(0);

    // Seed tables shared by every gamma: f = I_nu/nu^2 with its log2, and the
    // log2 width of each bin (for local power-law slopes).
    for (size_t j = 0; j < nu_size; ++j) {
        const Real f = I_nu_seed(j) / (nu_seed(j) * nu_seed(j));
        fv_th(j) = f;
        lg2fv_th(j) = (f > 0) ? fast_log2(f) : -con::inf;
        if (j + 1 < nu_size) {
            lg2r(j) = lg2_nu_seed(j + 1) - lg2_nu_seed(j);
            inv_lg2r(j) = (lg2r(j) != 0) ? 1 / lg2r(j) : 0;
        }
    }

    if (KN) {
        static thread_local Array cdf_th, ratio_th, corr_lat, lg2corr_lat;
        cdf_th.resize({nu_size});
        ratio_th.resize({nu_size});
        build_cdf_thomson(fv_th, lg2fv_th, nu_seed, dnu_seed, lg2r, inv_lg2r, cdf_th, ratio_th);

        // All gamma-nu products live on one log2 lattice (shared grid step):
        // one KN-correction evaluation per lattice node serves every (i, j).
        const size_t n_lat = gamma_mult * (static_cast<size_t>(gamma_size) - 1) + nu_mult * (nu_size - 1) + 1;
        corr_lat.resize({n_lat});
        lg2corr_lat.resize({n_lat});
        const Real lg2_base = fast_log2(gamma(0)) + lg2_nu_seed(0);
        for (size_t k = 0; k < n_lat; ++k) {
            compton_correction_pair(fast_exp2(lg2_base + lattice_quantum * static_cast<Real>(k)), corr_lat(k),
                                    lg2corr_lat(k));
        }

        for (int i = 0; i < gamma_size; ++i) {
            if (dN_e_boost(i) <= 0)
                continue;
            const Real gamma_i = gamma(i);
            build_cdf_KN(gamma_mult * static_cast<size_t>(i), gamma_i, corr_lat, lg2corr_lat, fv_th, lg2fv_th, nu_seed,
                         dnu_seed, lg2r, inv_lg2r, cdf_th, ratio_th, fv_buf, cdf_buf, ratio_buf);
            accumulate_IC(dN_e_boost(i), ic_idx0_ - 2 * static_cast<long>(gamma_mult) * i, nu_seed, dnu_seed, fv_buf,
                          cdf_buf, ratio_buf, I_buf, spec_size);
        }
    } else {
        build_cdf_thomson(fv_th, lg2fv_th, nu_seed, dnu_seed, lg2r, inv_lg2r, cdf_buf, ratio_buf);
        for (int i = 0; i < gamma_size; ++i) {
            if (dN_e_boost(i) <= 0)
                continue;
            accumulate_IC(dN_e_boost(i), ic_idx0_ - 2 * static_cast<long>(gamma_mult) * i, nu_seed, dnu_seed, fv_th,
                          cdf_buf, ratio_buf, I_buf, spec_size);
        }
    }

    const Real log2_scale = fast_log2(0.25 * con::sigmaT);
    interp_slope = Array::from_shape({spec_size - 1});
    for (size_t i = 0; i < spec_size; ++i) {
        log2_I_nu_IC(i) = fast_log2(I_buf[i]) + log2_nu_IC(i) + log2_scale;
        if (i > 0) {
            const Real dl = log2_nu_IC(i) - log2_nu_IC(i - 1);
            interp_slope(i - 1) = (dl != 0) ? (log2_I_nu_IC(i) - log2_I_nu_IC(i - 1)) / dl : 0;
        }
    }
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_I_nu(Real nu) {
    return fast_exp2(compute_log2_I_nu(fast_log2(nu)));
}

template <typename Electrons, typename Photons>
Real ICPhoton<Electrons, Photons>::compute_log2_I_nu(Real log2_nu) {
    if (!generated) {
        generate_spectrum();
    }

    size_t n = log2_nu_IC.size();

    // Band contract: evaluations must come from the same nu_obs/Doppler set the
    // clamp was derived from. A query outside the clamped grid but inside the
    // physical spectrum range would otherwise read a silent zero (high side) or
    // extrapolate a mid-spectrum slope in the wrong regime (low side). The clamp
    // is only an optimization, so on breach drop it and rebuild the full-range
    // spectrum: correctness self-heals and the cost degrades to pre-clamp
    // behavior (at most one rebuild per cell).
    if (n >= 2 && ((log2_nu > log2_nu_IC(n - 1) && log2_nu < log2_nu_theory_max_) ||
                   (log2_nu < log2_nu_IC(0) && log2_nu > log2_nu_theory_min_))) {
        nu_eval_min_ = 0;
        nu_eval_max_ = con::inf; // compute_grid_params then reproduces the full pre-clamp physical range
        generate_spectrum();
        n = log2_nu_IC.size();
    }

    if (n < 2 || log2_nu > log2_nu_IC(n - 1)) {
        return -con::inf;
    }

    size_t idx = 0;
    while (idx + 2 < n && log2_nu_IC(idx + 1) <= log2_nu) {
        ++idx;
    }

    return log2_I_nu_IC(idx) + (log2_nu - log2_nu_IC(idx)) * interp_slope(idx);
}

/**
 * @brief Build the IC photon grid, restricting each cell's output band to what the observer samples.
 * @details ``nu_eval_min_k`` / ``nu_eval_max_k`` give the comoving evaluation band per time index
 *          (the observer band shifted by that k-slice's Doppler extrema). Broadcast groups share a
 *          cell across angles but always preserve k, so per-k bands cover every query the copies
 *          receive. Empty arrays disable the clamp (full theoretical range, e.g. for details()).
 */
template <SynElectronModel Electrons, SynPhotonModel Photons>
ICPhotonGrid<Electrons, Photons>
generate_IC_photons(ElectronGrid<Electrons> const& electrons, PhotonGrid<Photons> const& photons, bool KN,
                    Coord const& coord, Array const& nu_eval_min_k = {}, Array const& nu_eval_max_k = {}) noexcept {
    size_t phi_size = electrons.shape()[0];
    size_t theta_size = electrons.shape()[1];
    size_t t_size = electrons.shape()[2];

    ICPhotonGrid<Electrons, Photons> IC_ph({phi_size, theta_size, t_size});

    const size_t phi_compute = (coord.symmetry != Symmetry::structured) ? 1 : phi_size;

    const bool precompute = (coord.symmetry != Symmetry::structured);
    assert((nu_eval_min_k.size() == 0 || nu_eval_min_k.size() == t_size) &&
           (nu_eval_max_k.size() == 0 || nu_eval_max_k.size() == t_size) &&
           "per-k band arrays must match the grid time dimension");
    const bool clamped = (nu_eval_min_k.size() == t_size) && (nu_eval_max_k.size() == t_size);

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : coord.theta_reps) {
            for (size_t k = 0; k < t_size; ++k) {
                const Real nu_lo = clamped ? nu_eval_min_k(k) : 0;
                const Real nu_hi = clamped ? nu_eval_max_k(k) : std::numeric_limits<Real>::infinity();
                IC_ph(i, j, k) = ICPhoton(electrons(i, j, k), photons(i, j, k), KN, nu_lo, nu_hi);
                if (precompute) {
                    IC_ph(i, j, k).build();
                }
            }
        }
    }

    broadcast_symmetry(IC_ph, coord);

    return IC_ph;
}

inline Real eta_rad_Thomson(Real gamma_m, Real gamma_c, Real p) {
    if (gamma_c < gamma_m) {
        return 1;
    } else {
        return fast_pow(gamma_c / gamma_m, 2 - p);
    }
}

inline Real compute_Thomson_Y(RadParams const& rad, Real gamma_m, Real gamma_c) {
    Real eta_e = eta_rad_Thomson(gamma_m, gamma_c, rad.p);
    Real b = eta_e * rad.eps_e / rad.eps_B;
    return 0.5 * (std::sqrt(1. + 4. * b) - 1.);
}

Real compute_gamma_c(Real t_comv, Real B, Real Y);

Real compute_syn_gamma_M(Real B, Real Y, Real p);

Real compute_syn_I_peak(Real B, Real p, Real column_den);

Real compute_syn_gamma_a(Real B, Real I_nu_peak, Real gamma_m, Real gamma_c, Real gamma_M, Real p,
                         InverseComptonY const& Ys, Real Y_c);

size_t determine_regime(Real gamma_a, Real gamma_c, Real gamma_m);

Real compute_gamma_0(Real Y0, Real gamma_m, Real gamma_m_hat);

void update_gamma_c_Thomson(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                            Real gamma_c_last);

void update_gamma_c_KN(Real& gamma_c, InverseComptonY& Ys, RadParams const& rad, Real B, Real t_com, Real gamma_m,
                       Real gamma_c_last);

void update_gamma_M(Real& gamma_M, InverseComptonY const& Ys, Real p, Real B);

Real compute_syn_gamma(Real nu, Real B);

template <SynElectronModel Electrons, SynPhotonModel Photons, typename Updater>
void IC_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock,
                Coord const& coord, Updater&& update_gamma_c) {
    const size_t phi_size = electrons.shape()[0];
    const size_t t_size = electrons.shape()[2];

    const size_t phi_compute = (coord.symmetry != Symmetry::structured) ? 1 : phi_size;

    for (size_t i = 0; i < phi_compute; ++i) {
        for (size_t j : coord.theta_reps) {
            for (size_t k = 0; k < t_size; ++k) {
                const Real t_com = shock.t_comv(i, j, k);
                const Real B = shock.B(i, j, k);

                auto& elec = electrons(i, j, k);
                auto& Ys = elec.Ys;
                const Real p = elec.p;
                const Real gamma_c_last = electrons(i, j, k > 0 ? k - 1 : 0).gamma_c;

                update_gamma_c(elec.gamma_c, Ys, shock.rad, B, t_com, elec.gamma_m, gamma_c_last);

                update_gamma_M(elec.gamma_M, Ys, p, B);

                cool_relic_electrons(electrons, shock, i, j, k);

                const Real I_nu_peak = compute_syn_I_peak(B, p, elec.column_den);
                elec.Y_c = Ys.gamma_spectrum(elec.gamma_c);
                elec.gamma_a =
                    compute_syn_gamma_a(B, I_nu_peak, elec.gamma_m, elec.gamma_c, elec.gamma_M, p, Ys, elec.Y_c);
                elec.regime = determine_regime(elec.gamma_a, elec.gamma_c, elec.gamma_m);
            }
        }
    }
    broadcast_symmetry(electrons, coord);
    generate_syn_photons(photons, shock, electrons, coord);
}

template <SynElectronModel Electrons, SynPhotonModel Photons>
void Thomson_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock,
                     Coord const& coord) {
    IC_cooling(electrons, photons, shock, coord, update_gamma_c_Thomson);
}

template <SynElectronModel Electrons, SynPhotonModel Photons>
void KN_cooling(ElectronGrid<Electrons>& electrons, PhotonGrid<Photons>& photons, Shock const& shock,
                Coord const& coord) {
    IC_cooling(electrons, photons, shock, coord, update_gamma_c_KN);
}

/// Whether IC cooling shaped this photon's spectrum. When false, the spectrum
/// segments are empty and Y_c = 0, so inverse_compton_correction is exactly 1 —
/// callers can skip it (and any log/exp wrapping around it).
template <typename Photon>
bool has_IC_correction(Photon const& ph) {
    return ph.Y_c > 0 || ph.Ys.Y_T > 0;
}

template <typename Photon>
Real inverse_compton_correction(Photon const& ph, Real nu) {
    if (!has_IC_correction(ph)) {
        return 1.0;
    }
    return (1. + ph.Y_c) / (1 + ph.Ys.nu_spectrum(nu));
}
