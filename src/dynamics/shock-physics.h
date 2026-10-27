//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "lanes.h"
#include "shock.h"
#include "util/macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @defgroup ShockUtilities Shock Utilities
 * @brief Inline functions used in shock calculations.
 * @details This section defines a set of inline functions used in shock calculations. These functions compute
 *          various physical quantities such as the comoving magnetic field (via the Weibel instability),
 *          thermal energy density, time derivatives, jet width derivative, downstream number density, fluid
 *          velocities, and update the shock state.
 * <!-- ************************************************************************************** -->
 */

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the downstream four-velocity for a given relative Lorentz factor and magnetization parameter.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The downstream four-velocity in the shock frame
 * @details The calculation handles both magnetized (sigma > 0) and non-magnetized (sigma = 0) cases using
 *          different equations based on jump conditions across the shock front.
 * <!-- ************************************************************************************** -->
 */
Real compute_downstr_4vel(Real gamma_rel, Real sigma) noexcept;

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the upstream four-velocity from downstream four-velocity and relative Lorentz factor.
 * @param u_down Downstream four-velocity
 * @param gamma_rel Relative Lorentz factor
 * @return The upstream four-velocity in the shock frame
 * <!-- ************************************************************************************** -->
 */
inline Real compute_upstr_4vel(Real u_down, Real gamma_rel) noexcept {
    return std::sqrt((1 + u_down * u_down) * std::max((gamma_rel - 1) * (gamma_rel + 1), 0.0)) + u_down * gamma_rel;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the ratio of upstream to downstream four-velocity across the shock front.
 * @param gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @param sigma_upstr Magnetization parameter (ratio of magnetic to rest-mass energy density)
 * @return The ratio of upstream to downstream four-velocity
 * @details This ratio is a key parameter in determining various shock properties, such as compression ratio
 *          and jump conditions for density, pressure, and magnetic field.
 * <!-- ************************************************************************************** -->
 */
inline Real compute_4vel_jump(Real gamma_rel, Real sigma_upstr) noexcept {
    const Real u_down_s_ = compute_downstr_4vel(gamma_rel, sigma_upstr);
    const Real u_up_s_ = compute_upstr_4vel(u_down_s_, gamma_rel);
    Real ratio_u = u_up_s_ / u_down_s_;
    if (u_down_s_ == 0.) {
        ratio_u = 4 * gamma_rel; // (g_hat*gamma_rel+1)/(g_hat-1)
    }
    return ratio_u;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the sound speed in the shocked medium based on the relative Lorentz factor.
 * @param Gamma_rel Relative Lorentz factor between upstream and downstream regions
 * @return The sound speed as a fraction of light speed
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline Value compute_sound_speed(Value const& Gamma_rel) noexcept {
    using std::max, std::sqrt;
    const Value ad_idx = physics::thermo::adiabatic_idx(Gamma_rel);
    return sqrt(max(ad_idx * (ad_idx - 1) * (Gamma_rel - 1) / (1 + (Gamma_rel - 1) * ad_idx), Value(0))) * con::c;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the effective Lorentz factor accounting for the adiabatic index.
 * @param adx Adiabatic index of the medium
 * @param Gamma Bulk Lorentz factor
 * @return The effective Lorentz factor
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline constexpr Value compute_effective_Gamma(Value const& adx, Value const& Gamma) noexcept {
    return (adx * Gamma * Gamma - adx + 1) / Gamma;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the derivative of the effective Lorentz factor with respect to the bulk Lorentz factor.
 * @param adx Adiabatic index of the medium
 * @param Gamma Bulk Lorentz factor
 * @return The derivative of the effective Lorentz factor
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline constexpr Value compute_effective_Gamma_dGamma(Value const& adx, Value const& Gamma) noexcept {
    const Value Gamma2 = Gamma * Gamma;
    return (adx * Gamma2 + adx - 1) / Gamma2;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the comoving magnetic field using the Weibel instability mechanism.
 * @param eps_B Fraction of thermal energy in magnetic fields
 * @param e_thermal Thermal energy density
 * @return The comoving magnetic field strength
 * <!-- ************************************************************************************** -->
 */
inline Real compute_comv_weibel_B(Real eps_B, Real e_thermal) noexcept {
    return std::sqrt(8 * con::pi * eps_B * e_thermal);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of radius (dr/dt) based on the shock velocity.
 * @param beta Shock velocity as a fraction of light speed
 * @return The rate of change of radius with respect to observer time
 * <!-- ************************************************************************************** -->
 */
inline constexpr Real compute_dr_dt(Real beta) noexcept {
    return (beta * con::c) / (1 - beta);
}

/// Overload using four-velocity u and Lorentz factor Gamma directly.
/// dr/dt = βc/(1-β) = uc(Γ+u), exploiting Γ²-u²=1 to avoid divisions and 1-β cancellation.
template <typename Value = Real>
inline constexpr Value compute_dr_dt(Value const& Gamma, Value const& u) noexcept {
    return u * (Gamma + u) * con::c;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of theta (dθ/dt) for jet spreading.
 * @param theta_s Typical spreading angle of the jet
 * @param theta Theta of the current grid point
 * @param drdt Time derivative of radius
 * @param r Current radius
 * @param Gamma Current bulk Lorentz factor
 * @return The rate of change of the half-opening angle
 * <!-- ************************************************************************************** -->
 */
inline Real compute_dtheta_dt(Real theta_s, Real /*theta*/, Real drdt, Real r, Real Gamma) noexcept {
    constexpr Real Q = 7;
    const Real u2 = (Gamma - 1) * (Gamma + 1);
    const Real u = std::sqrt(u2);
    const Real f = 1 / (1 + u * theta_s * Q);
    return drdt / (2 * Gamma * r) * std::sqrt((2 * u2 + 3) / (4 * u2 + 3)) * f;
}

/// Overload using precomputed four-velocity terms to avoid redundant sqrt in hot loops.
template <typename Value = Real>
inline Value compute_dtheta_dt(Real theta_s, Value const& /*theta*/, Value const& drdt, Value const& r,
                               Value const& Gamma, Value const& u, Value const& u2) noexcept {
    using std::sqrt;
    constexpr Real Q = 7;
    const Value f = 1 / (1 + u * theta_s * Q);
    return drdt / (2 * Gamma * r) * sqrt((2 * u2 + 3) / (4 * u2 + 3)) * f;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the time derivative of comoving time (dt_comv/dt) based on the Lorentz factor.
 * @param Gamma Bulk Lorentz factor
 * @param beta Bulk velocity as a fraction of light speed
 * @return The rate of change of comoving time with respect to observer time
 * <!-- ************************************************************************************** -->
 */
inline constexpr Real compute_dt_dt_comv(Real Gamma, Real beta) noexcept {
    return 1 / (Gamma * (1 - beta));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the upstream magnetic field.
 * @param rho_up Upstream density
 * @param sigma Magnetization parameter
 * @return The upstream magnetic field
 * <!-- ************************************************************************************** -->
 */
inline Real compute_upstr_B(Real rho_up, Real sigma) noexcept {
    return std::sqrt((4 * con::pi * con::c2) * sigma * rho_up);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the relative Lorentz factor between two frames with given Lorentz factors.
 * @param gamma1 First frame's Lorentz factor
 * @param gamma2 Second frame's Lorentz factor
 * @return The relative Lorentz factor between the two frames
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline Value compute_rel_Gamma(Value const& gamma1, Value const& gamma2) noexcept {
    using std::max, std::sqrt;
    // Cancellation-free form. The naive gamma1*gamma2 - u1*u2 subtracts two
    // O(gamma^2) numbers to produce 1 + O(1e-6) near equal velocities; the
    // exact identity Gamma_rel - 1 = (gamma1 - gamma2)^2 / (gamma1*gamma2 - 1
    // + u1*u2) keeps full precision there (both denominator terms positive).
    const Value u1u2 = sqrt(max((gamma1 - 1) * (gamma1 + 1) * (gamma2 - 1) * (gamma2 + 1), Value(0)));
    const Value d = gamma1 - gamma2;
    const Value denom = gamma1 * gamma2 - 1 + u1u2;
    return where(denom <= 0, 1., 1 + d * d / denom);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the shock heating rate.
 * @param Gamma_rel Relative Lorentz factor
 * @param mdot Mass accretion rate
 * @return The shock heating rate
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline constexpr Value compute_shock_heating_rate(Value const& Gamma_rel, Value const& mdot) noexcept {
    return mdot * (Gamma_rel - 1) * con::c2;
}

/// Medium density with compile-time lane dispatch: the native lane overload when the
/// medium provides one, else per-lane scalar evaluation (std::function / python media).
template <typename Value, typename Medium>
inline Value lane_rho(Medium const& medium, Real phi, Value const& theta, Value const& r) noexcept {
    if constexpr (requires { medium.rho(phi, theta, r); }) {
        return medium.rho(phi, theta, r);
    } else {
        return lane::from_lanes<Value>([&](size_t w) { return medium.rho(phi, lane::get(theta, w), lane::get(r, w)); });
    }
}

/// Clamp-form smoothstep (identical values to the branch form for every input).
template <typename Value>
inline Value smoothstep(Real edge0, Real edge1, Value const& x) noexcept {
    using std::clamp;
    const Value t = clamp((x - edge0) / (edge1 - edge0), Value(0), Value(1));
    return t * t * (3.0 - 2.0 * t);
}

template <typename Value = Real>
inline Value compute_adiabatic_cooling_rate2(Value const& ad_idx, Value const& r, Value const& x, Value const& u,
                                             Value const& drdt, Value const& dxdt) noexcept {
    const Value dlnvdt = 2 * drdt / r;
    const Value with_shell = where(x > 0, dlnvdt + dxdt / x, dlnvdt);
    return -(ad_idx - 1) * with_shell * u;
}
/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the sound-speed expansion rate of the comoving shell width.
 * @details Radial widening of the shell at the local sound speed — distinct from
 *          lateral jet spreading (compute_dtheta_dt / Ejecta::spreading).
 * @param Gamma_rel Relative Lorentz factor
 * @param dtdt_comv Rate of change of comoving time with respect to burst time
 * @return The sound-speed expansion rate of the comoving shell width
 * <!-- ************************************************************************************** -->
 */
template <typename Value = Real>
inline Value compute_shell_sound_expansion_rate(Value const& Gamma_rel, Value const& dtdt_comv) noexcept {
    return compute_sound_speed(Gamma_rel) * dtdt_comv;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Radiative efficiency of a shock, precomputed from the radiation parameters.
 * @details Single owner of the gamma_m/gamma_c efficiency formula shared by all shock equations
 *          (forward, reverse pair, and simple). Construct once per equation, evaluate per ODE step.
 * <!-- ************************************************************************************** -->
 */
class RadiativeEfficiency {
  public:
    RadiativeEfficiency() = default;

    explicit RadiativeEfficiency(RadParams const& rad) noexcept
        : gamma_m_coeff_((rad.p - 2) / (rad.p - 1) * rad.eps_e * con::mp / con::me / rad.xi_e),
          // gamma_c = 6π me c / (σ_T B'² t') with B'² = 8π eps_B e_th (the synchrotron-side
          // convention, compute_gamma_c). The SSC (1+Y) term is deliberately omitted here to
          // keep the dynamics decoupled from the radiation-side Compton solve.
          gamma_c_coeff_(6 * con::pi * con::me * con::c / con::sigmaT / (8 * con::pi * rad.eps_B)),
          eps_e_(rad.radiative ? rad.eps_e : 0),
          p_(rad.p) {}

    /**
     * @brief Evaluates the radiative efficiency at the current shock state.
     * @param t_comv Comoving time
     * @param Gamma_th Thermal Lorentz factor
     * @param e_th Internal energy density
     * @return The radiative efficiency
     */
    template <typename Value>
    Value operator()(Value const& t_comv, Value const& Gamma_th, Value const& e_th) const noexcept {
        using std::sqrt;
        if (eps_e_ == 0) { // adiabatic fireball (radiative = false) or no accelerated electrons
            return 0.;
        }
        const Value gamma_m = gamma_m_coeff_ * (Gamma_th - 1) + 1;
        const Value gamma_bar = gamma_c_coeff_ / (e_th * t_comv);
        // Trans-relativistic correction (same as compute_gamma_c): cooling rate ∝ γ²β², so
        // gamma_c -> 1 smoothly as gamma_bar -> 0 instead of a hard clamp.
        const Value gamma_c = 0.5 * (gamma_bar + sqrt(gamma_bar * gamma_bar + 4));
        const Value ratio = gamma_m / gamma_c;
        if (p_ > 2) { // slow cooling where gamma_m < gamma_c, per lane
            return lane::from_lanes<Value>([&](size_t w) {
                const Real r = lane::get(ratio, w);
                return r < 1 ? eps_e_ * fast_pow(r, p_ - 2) : eps_e_;
            });
        }
        return eps_e_; // fast cooling or p<=2
    }

  private:
    Real gamma_m_coeff_{0}; ///< Precomputed: (p-2)/(p-1) * eps_e * mp/me / xi_e
    Real gamma_c_coeff_{0}; ///< Precomputed: 6π me c / (σ_T * 8π * eps_B)
    Real eps_e_{0};         ///< Fraction of shock energy in electrons
    Real p_{0};             ///< Electron spectral index
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the thermal Lorentz factor from the thermal energy and mass.
 * @param U_th Thermal energy
 * @param mass Mass
 * @param limiter
 * @return The thermal Lorentz factor
 * <!-- ************************************************************************************** -->
 */
inline Real compute_Gamma_therm(Real U_th, Real mass, bool limiter = false) noexcept {
    if (mass == 0) [[unlikely]] {
        return 1;
    } else [[likely]] {
        const Real Gamma_th = U_th / (mass * con::c2) + 1;
        if (limiter && Gamma_th < con::gamma_therm_cut) {
            return 1;
        } else {
            return Gamma_th;
        }
    }
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Saves the current state of the shock.
 * @param shock Reference to the Shock object to be updated
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param t_comv Comoving time
 * @param r Radius
 * @param theta Angle
 * @param Gamma Lorentz factor
 * @param Gamma_th Thermal Lorentz factor
 * @param B Magnetic field strength
 * @param mass Mass
 * <!-- ************************************************************************************** -->
 */

inline void write_shock_state(Shock& shock, size_t i, size_t j, size_t k, Real t_comv, Real r, Real theta, Real Gamma,
                              Real Gamma_th, Real B, Real mass) {
    shock.t_comv(i, j, k) = t_comv;
    shock.r(i, j, k) = r;
    shock.theta(i, j, k) = theta;
    shock.Gamma(i, j, k) = Gamma;
    shock.Gamma_th(i, j, k) = Gamma_th;
    shock.B(i, j, k) = B;
    shock.N_p(i, j, k) = mass / con::mp;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the compression ratio across the shock.
 * @param Gamma_upstr Lorentz factor upstream
 * @param Gamma_downstr Lorentz factor downstream
 * @param sigma_upstr Upstream magnetization
 * @return The compression ratio
 * <!-- ************************************************************************************** -->
 */
inline Real compute_compression(Real Gamma_upstr, Real Gamma_downstr, Real sigma_upstr) noexcept {
    const Real Gamma_rel = compute_rel_Gamma(Gamma_upstr, Gamma_downstr);
    return compute_4vel_jump(Gamma_rel, sigma_upstr);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Computes the downstream magnetic field strength.
 * @param eps_B Magnetic energy fraction
 * @param rho_upstr Upstream rest mass density
 * @param B_upstr Upstream magnetic field strength
 * @param Gamma_th Thermal Lorentz factor
 * @param comp_ratio Compression ratio
 * @return The downstream magnetic field strength
 * <!-- ************************************************************************************** -->
 */
inline Real compute_downstr_B(Real eps_B, Real rho_upstr, Real B_upstr, Real Gamma_th, Real comp_ratio) noexcept {
    const Real rho_downstr = rho_upstr * comp_ratio;

    const Real e_th = (Gamma_th - 1) * rho_downstr * con::c2;

    return compute_comv_weibel_B(eps_B, e_th) + B_upstr * comp_ratio;
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Sets a stopping shock state when the Lorentz factor drops below a threshold.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param shock Reference to the Shock object to be updated
 * @param state0 Initial state to be used for some parameters
 * <!-- ************************************************************************************** -->
 */
template <typename State>
inline void set_stopping_shock(size_t i, size_t j, Shock& shock, State const& state0);

//========================================================================================================
//                                  template function implementation
//========================================================================================================
template <typename State>
inline void set_stopping_shock(size_t i, size_t j, Shock& shock, State const& state0) {
    xt::view(shock.t_comv, i, j, xt::all()) = state0.t_comv;
    xt::view(shock.r, i, j, xt::all()) = state0.r;
    xt::view(shock.theta, i, j, xt::all()) = state0.theta;
    xt::view(shock.Gamma, i, j, xt::all()) = 1;
    xt::view(shock.Gamma_th, i, j, xt::all()) = 1;
    xt::view(shock.B, i, j, xt::all()) = 0;
    xt::view(shock.N_p, i, j, xt::all()) = 0;
}

/// @brief Runs `body(js, thetas, n_live)` once per lockstep batch of theta rows at phi index i.
/// @details Owns the batch-formation policy: lone rows get width-1 arrays (so the scalar
///          Value = Real instantiation skips the padding overhead), and trailing lanes of a
///          partial batch pad with the last live row.
template <typename F>
void for_each_theta_batch(Coord const& coord, F&& body) {
    constexpr size_t lanes = lane::width;
    auto const& reps = coord.theta_reps;
    for (size_t g = 0; g < reps.size(); g += lanes) {
        const size_t n_live = std::min(lanes, reps.size() - g);
        if (n_live == 1) {
            body(std::array<size_t, 1>{reps[g]}, std::array<Real, 1>{coord.theta(reps[g])}, size_t{1});
        } else {
            std::array<size_t, lanes> js;
            std::array<Real, lanes> thetas;
            for (size_t w = 0; w < lanes; ++w) {
                js[w] = reps[g + std::min(w, n_live - 1)];
                thetas[w] = coord.theta(js[w]);
            }
            body(js, thetas, n_live);
        }
    }
}

/// @brief Integration start shared by every theta row at phi index i: below each row's first
///        output time and a fraction of its deceleration time. One value per phi keeps a row's
///        solution independent of which rows share its lockstep batch.
template <typename Ejecta, typename Medium>
Real shock_start_time(Coord const& coord, size_t i, Ejecta const& jet, Medium const& medium, Real t_floor) {
    Real t0 = t_floor;
    for (size_t j : coord.theta_reps) {
        const Real t_dec = estimate_t_dec(jet, medium, coord.phi(i), coord.theta(j));
        t0 = std::min({t0, coord.t(i, j, 0), 0.1 * t_dec});
    }
    return t0;
}

/// @brief Output-time rows of one lockstep batch and the last output time of its live rows.
template <size_t K>
auto batch_time_span(Coord const& coord, size_t i, std::array<size_t, K> const& js, size_t n_live) {
    struct Span {
        std::array<Real const*, K> t_row; ///< each lane's output-time row
        Real t_end;                       ///< integration end
    } span;

    const size_t t_size = coord.t.shape()[2];
    for (size_t w = 0; w < K; ++w) {
        span.t_row[w] = &coord.t(i, js[w], 0);
    }
    span.t_end = 0;
    for (size_t w = 0; w < n_live; ++w) {
        span.t_end = std::max(span.t_end, span.t_row[w][t_size - 1]);
    }
    return span;
}

/// @brief Marks lanes whose initial Lorentz factor is at/below `Gamma_stop` as dead, invoking
///        `on_stopped(w)` for each, and parks dead/padding lanes on the first live lane.
/// @return Per-lane liveness; all-false means no lane is live and the batch can be skipped.
template <size_t K, typename State, typename F>
std::array<bool, K> mark_live_lanes(State& state, size_t n_live, Real Gamma_stop, F&& on_stopped) {
    std::array<bool, K> live{};
    size_t first_live = K;
    for (size_t w = 0; w < n_live; ++w) {
        if (lane::get(state.Gamma, w) <= Gamma_stop) {
            on_stopped(w);
        } else {
            live[w] = true;
            if (first_live == K) {
                first_live = w;
            }
        }
    }
    if (first_live < K) {
        park_dead_lanes(state, live, first_live);
    }
    return live;
}

/// @brief Simpson's rule integration in log-space: ∫ f(r') dr' via u = ln(r').
/// Integrates from r * e^{-18} to r with N=32 panels.
template <typename Func>
Real simpson_logspace(Func const& f, Real r) {
    constexpr size_t N = 32; // must be even
    const Real u_max = std::log(r);
    const Real u_min = u_max - 18; // start from r * e^{-18} ≈ r * 1.5e-8
    const Real h = (u_max - u_min) / N;

    Real sum = f(u_min) + f(u_max);
    for (size_t i = 1; i < N; i += 2)
        sum += 4 * f(u_min + i * h);
    for (size_t i = 2; i < N; i += 2)
        sum += 2 * f(u_min + i * h);

    return sum * h / 3;
}

/// @brief Computes enclosed mass per solid angle.
template <typename Func>
Real enclosed_mass(Func const& rho, Real r) {
    return simpson_logspace(
        [&](Real u) {
            const Real ri = std::exp(u);
            return rho(ri) * ri * ri * ri; // ρ(r') * r'³ (Jacobian from u = ln r')
        },
        r);
}

/// @brief Swept mass weighted by adiabatic cooling: ∫ ρ(r') r'² (r'/r)^{3(ad-1)} dr'.
template <typename Func>
Real enclosed_thermal_integral(Func const& rho, Real r, Real ad_idx) {
    const Real cooling_exp = 3 * (ad_idx - 1);
    return simpson_logspace(
        [&](Real u) {
            const Real ri = std::exp(u);
            return rho(ri) * ri * ri * ri * std::pow(ri / r, cooling_exp);
        },
        r);
}

/// @brief Compute enclosed thermal energy per solid angle,
/// accounting for adiabatic cooling and radiative losses of material swept at earlier radii.
template <typename Func>
Real enclosed_thermal_energy(Func const& rho, Real r, Real Gamma, Real ad_idx, Real eps_e) {
    return (1 - eps_e) * (Gamma - 1) * con::c2 * enclosed_thermal_integral(rho, r, ad_idx);
}

/// @brief Fast path for enclosed mass with compile-time medium detection.
template <typename Medium>
Real enclosed_mass_medium(Medium const& medium, Real phi, Real theta, Real r) {
    if constexpr (HasMass<Medium>) {
        return medium.mass(r);
    } else {
        return enclosed_mass([&](Real r_) { return medium.rho(phi, theta, r_); }, r);
    }
}

/// @brief Fast path for the cooling-weighted swept mass with compile-time medium detection.
template <typename Medium>
Real enclosed_thermal_integral_medium(Medium const& medium, Real phi, Real theta, Real r, Real ad_idx) {
    using MediumT = std::remove_cvref_t<Medium>;

    if constexpr (std::is_same_v<MediumT, ISM>) {
        const Real rho = medium.rho(phi, theta, r);
        const Real cooling_exp = 3 * (ad_idx - 1);
        const Real pow_exp = 3 + cooling_exp;
        const Real x0 = std::exp(-18.0);
        const Real attenuation = 1 - std::pow(x0, pow_exp);
        return rho * r * r * r * attenuation / pow_exp;
    } else {
        return enclosed_thermal_integral([&](Real r_) { return medium.rho(phi, theta, r_); }, r, ad_idx);
    }
}

/// @brief Fast path for enclosed thermal energy with compile-time medium detection.
template <typename Medium>
Real enclosed_thermal_energy_medium(Medium const& medium, Real phi, Real theta, Real r, Real Gamma, Real ad_idx,
                                    Real eps_e) {
    return (1 - eps_e) * (Gamma - 1) * con::c2 * enclosed_thermal_integral_medium(medium, phi, theta, r, ad_idx);
}

/// Start-radius state that conserves the blast-wave energy with the mass already swept there.
struct InitEnergyState {
    Real Gamma; ///< Lorentz factor (Gamma0 would overcount the swept energy by ~Gamma0 m2/m_jet: 10% in winds)
    Real U_th;  ///< thermal energy per solid angle at the start radius
};

/// @brief Solves Gamma (m_jet + m2) + Gamma_eff U/c² + E_rad/c² = Gamma0 m_jet + m2 for the start state.
/// @param thermal_integral Cooling-weighted swept mass at the adiabatic index of Gamma0 (its Gamma dependence
///        is second order, so one evaluation serves every Newton step: 33 density calls per row, not per step).
/// @details The pre-start history (thermal and radiated shares) is evaluated at the mean of Gamma0 and Gamma.
inline InitEnergyState initial_state_energy_conserving(Real Gamma0, Real m_jet, Real m2, Real eps_e,
                                                       Real thermal_integral) noexcept {
    Real Gamma = Gamma0;
    Real U = 0; // per c²
    const Real integral = thermal_integral;
    for (int it = 0; it < 16; ++it) {
        const Real Gamma_h = 0.5 * (Gamma0 + Gamma);
        const Real Gamma_eff_h = compute_effective_Gamma(physics::thermo::adiabatic_idx(Gamma_h), Gamma_h);
        U = (1 - eps_e) * (Gamma_h - 1) * integral;
        const Real E_rad = eps_e * Gamma_eff_h * (Gamma_h - 1) * m2; // per c²

        const Real ad_idx = physics::thermo::adiabatic_idx(Gamma);
        const Real Gamma_eff = compute_effective_Gamma(ad_idx, Gamma);
        const Real g = Gamma * (m_jet + m2) + Gamma_eff * U + E_rad - Gamma0 * m_jet - m2;
        const Real dg = (m_jet + m2) + compute_effective_Gamma_dGamma(ad_idx, Gamma) * U +
                        0.5 * (Gamma_eff * (1 - eps_e) * integral + eps_e * Gamma_eff_h * m2);
        const Real step = g / dg;
        Gamma = std::max(Gamma - step, 1.0);
        if (std::abs(step) <= 1e-12 * Gamma) {
            break;
        }
    }
    return {Gamma, U * con::c2};
}
