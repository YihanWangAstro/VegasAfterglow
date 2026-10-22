//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "../util/macros.h"
#include "lanes.h"
#include "shock-physics.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct ReverseState
 * @brief State vector for the coupled forward/reverse shock system, one shock pair per lane.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Value = Real>
struct ReverseState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;   ///< Whether ejecta has mass injection
    static constexpr bool energy_inject = HasDedt<Ejecta>; ///< Whether ejecta has energy injection
    static constexpr size_t array_size = 12;

    MAKE_THIS_ODEINT_STATE_T(ReverseState, data, array_size, Value)

    union {
        struct {
            Value Gamma;  ///< Lorentz factor of the shocked region
            Value x4;     ///< Comoving frame width of region 4
            Value x3;     ///< Comoving frame width of region 3
            Value m2;     ///< Shocked medium mass per solid angle
            Value m3;     ///< Shocked ejecta mass per solid angle
            Value U2_th;  ///< internal energy per solid angle in region 2
            Value U3_th;  ///< internal energy per solid angle in region 3
            Value r;      ///< Radius
            Value t_comv; ///< Comoving time
            Value theta;  ///< Angular coordinate theta
            Value eps4;   ///< energy  per solid angle in region 4
            Value m4;     ///< mass per solid angle in region 4
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class FRShockEqn
 * @brief The forward/reverse shock pair equation for a given jet and medium; integrates
 *        `lanes` theta rows in lockstep when Value is a lane batch.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium, typename Value = Real>
class FRShockEqn {
  public:
    using State = ReverseState<Ejecta, Value>;
    using ScalarState = ReverseState<Ejecta, Real>;
    static constexpr size_t lanes = lane::width_of<Value>; ///< theta rows integrated in lockstep

    FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, std::array<Real, lanes> const& thetas,
               RadParams const& rad_fwd, RadParams const& rad_rvs); // rad_rvs reserved (reverse Shock carries its own)

    /// Scalar convenience constructor (single theta row).
    FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, RadParams const& rad_fwd,
               RadParams const& rad_rvs)
        requires(lanes == 1)
        : FRShockEqn(medium, ejecta, phi, std::array<Real, 1>{theta}, rad_fwd, rad_rvs) {}

    Medium const& medium;                       ///< Reference to the medium properties
    Ejecta const& ejecta;                       ///< Reference to the jet properties
    RadParams const rad_fwd;                    ///< Radiation parameters for forward shock
    Real const phi{0};                          ///< Angular coordinate phi
    std::array<Real, lanes> const theta0_lanes; ///< Initial theta of each lane's row
    Value Gamma4{1};                            ///< Initial Lorentz factor of the jet (per lane)

    // Crossing snapshots, written per lane at that lane's crossing time.
    std::array<Real, lanes> B3_ordered_x{}; ///< Ordered magnetic field in region 3 at crossing
    std::array<Real, lanes> V3_comv_x{};    ///< Comoving volume in region 3 at crossing
    std::array<Real, lanes> rho3_x{};       ///< Density in region 3 at crossing

    /// The reverse shock ODE system (all lanes at once).
    void operator()(State const& state, State& diff, Real t) const noexcept;

    /// Initial conditions for all lanes at engine time t0.
    void set_init_state(State& state, Real t0) const noexcept;

    /// Shell magnetization sigma = eps / (Gamma m c^2) - 1, per lane.
    Value compute_shell_sigma(State const& state) const noexcept;

    /// Shell magnetization of one lane's scalar state.
    Real compute_shell_sigma(ScalarState const& state, size_t w) const noexcept;

    /// Saves lane w's crossing snapshot.
    void save_cross_state(ScalarState const& state, size_t w);

    /// Fraction of injection still active per lane: 0 after injection, 1 during full injection.
    inline Value injection_efficiency(State const& diff) const noexcept;

    /// Smooth injection shutdown window over 50% of T0 (t is shared: lane-uniform).
    Real injection_window(Real t) const noexcept { return smoothstep(ejecta.T0 * 1.5, ejecta.T0 * 0.5, t); }

    /// Whether lane w's reverse shock has finished crossing the ejecta.
    inline bool crossing_complete(Real m3, Real m4, Real t, size_t w) const noexcept;

  private:
    /// Initial conditions of a single lane (the scalar body assembled by set_init_state).
    void set_init_state_lane(ScalarState& state, Real t0, size_t w) const noexcept;

    /// Evaluates a scalar profile at each row's initial theta: lane w of the result is f(theta0_lanes[w]).
    /// (theta is frozen in the pair solve — diff.theta = 0 — so this is also the current angle.)
    template <typename F>
    Value at_theta0(F&& f) const noexcept {
        return lane::from_lanes<Value>([&](size_t w) { return f(theta0_lanes[w]); });
    }

    inline Value compute_dGamma_dt(State const& state, State const& diff, Real t, Value const& Gamma34) const noexcept;

    inline Value compute_dU2_dt(State const& state, State const& diff, Value const& rho) const noexcept;

    inline Value compute_dU3_dt(State const& state, State const& diff, Value const& Gamma34) const noexcept;

    inline Value compute_dx3_dt(State const& state, State const& diff, Value const& f, Value const& Gamma34,
                                Value const& sigma, Value const& comp_ratio) const noexcept;

    inline Value compute_dx4_dt(State const& state, State const& diff, Value const& f) const noexcept;

    inline Value compute_dm3_dt(State const& state, State const& diff, Value const& f,
                                Value const& comp_ratio) const noexcept;

    inline Value compute_deps4_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Value compute_dm4_dt(State const& state, State const& diff, Real t) const noexcept;

    Value deps0_dt_{0};               ///< Ejecta energy injection rate (per lane)
    Value dm0_dt_{0};                 ///< Ejecta mass injection rate (per lane)
    Value u4_{0};                     ///< Four-velocity of the unshocked ejecta (per lane)
    Value beta4_{0};                  ///< Velocity fraction of the unshocked ejecta (per lane)
    Value cs4_{0};                    ///< Sound speed of the unshocked ejecta (per lane)
    RadiativeEfficiency eps_rad_fwd_; ///< Radiative efficiency evaluator for the forward component
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Function templates for shock generation
 * @details Declares interfaces to generate forward shocks (2D and 3D) and forward/reverse shock pairs.
 * <!-- ************************************************************************************** -->
 */
using ShockPair = std::pair<Shock, Shock>;

/**
 * <!-- ************************************************************************************** -->
 * @brief Generates a pair of forward and reverse shocks.
 * @details Creates two Shock objects and solves the shock evolution for each grid point.
 * @param coord Coordinate system definition
 * @param medium The medium through which the shock propagates
 * @param jet The jet (ejecta) driving the shock
 * @param rad_fwd Radiation params for forward shock
 * @param rad_rvs Radiation params for reverse shock
 * @param rtol Relative tolerance for ODE solver
 * @return A pair of Shock objects {forward_shock, reverse_shock}
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
ShockPair generate_shock_pair(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_fwd,
                              RadParams const& rad_rvs, Real rtol = 1e-5);

//========================================================================================================
//                                  template function implementation
//========================================================================================================

#include "../src/dynamics/reverse-shock.tpp"
