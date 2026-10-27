//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "lanes.h"
#include "shock-physics.h"
#include "shock.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct ForwardState
 * @brief State vector structure for forward shock calculations.
 * @details Template structure that defines the state vector for forward shock calculations, adapting its size
 *          based on if the Ejecta class supports mass and energy injection methods.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Value = Real>
struct ForwardState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;   ///< whether Ejecta class has dmdt method
    static constexpr bool energy_inject = HasDedt<Ejecta>; ///< whether Ejecta class has dedt method
    /// use the least fixed array size for integrator efficiency
    static constexpr size_t array_size = 6 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0);

    MAKE_THIS_ODEINT_STATE_T(ForwardState, data, array_size, Value)

    union {
        struct {
            Value Gamma;  ///< Lorentz factor
            Value m2;     ///< swept mass
            Value U2_th;  ///< internal energy per solid angle
            Value r;      ///< radius
            Value t_comv; ///< comoving time
            Value theta;  ///< angle

            /// shell energy density per solid angle
            [[no_unique_address]] std::conditional_t<energy_inject, Value, class Empty> eps_jet;

            /// shell mass per solid angle
            [[no_unique_address]] std::conditional_t<mass_inject, Value, class Empty> m_jet;
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class ForwardShockEqn
 * @brief Represents the forward shock equation for a given jet and medium.
 * @details It defines a state vector (with variable size based on template parameters) and overloads operator()
 *          to compute the derivatives of the state with respect to source time. It also declares helper functions for
 *          the derivatives. This class implements the physical equations governing the forward shock evolution.
 *
 * @tparam Ejecta The ejecta class template parameter
 * @tparam Medium The medium class template parameter
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium, typename Value = Real>
class ForwardShockEqn {
  public:
    using State = ForwardState<Ejecta, Value>;
    using ScalarState = ForwardState<Ejecta, Real>;
    static constexpr size_t lanes = lane::width_of<Value>; ///< theta rows integrated in lockstep

    /**
     * <!-- ************************************************************************************** -->
     * @brief ForwardShockEqn constructor
     * @details Initializes the forward shock equation with the given medium, ejecta, and parameters.
     * @param medium The medium through which the shock propagates
     * @param ejecta The ejecta driving the shock
     * @param phi Azimuthal angle
     * @param thetas Initial theta of each lane's row
     * @param rad_params Radiation parameters
     * @param theta_s Critical angle for jet spreading
     * <!-- ************************************************************************************** -->
     */
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, std::array<Real, lanes> const& thetas,
                    RadParams const& rad_params, Real theta_s);

    /// Scalar convenience constructor (single theta row).
    ForwardShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, RadParams const& rad_params,
                    Real theta_s)
        requires(lanes == 1)
        : ForwardShockEqn(medium, ejecta, phi, std::array<Real, 1>{theta}, rad_params, theta_s) {}

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivatives of the state variables with respect to engine time t.
     * @details Implements the system of ODEs that describe the evolution of the forward shock.
     * @param state Current state of the system
     * @param diff Output derivatives to be populated
     * @param t Current time
     * <!-- ************************************************************************************** -->
     */
    void operator()(State const& state, State& diff, Real t) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Set the initial conditions for the forward shock ODE solver.
     * @details Computes the initial state of the shock based on the ejecta properties and ambient medium.
     * @param state State vector to initialize
     * @param t0 Initial time
     * <!-- ************************************************************************************** -->
     */
    void set_init_state(State& state, Real t0) const noexcept;

    Medium const& medium;                       ///< Reference to the ambient medium properties
    Ejecta const& ejecta;                       ///< Reference to the ejecta properties
    Real const phi{0};                          ///< Angular coordinate phi in the jet frame
    std::array<Real, lanes> const theta0_lanes; ///< Initial theta of each lane's row
    RadParams const rad;                        ///< Radiation parameters

  private:
    /// Initial conditions of a single lane (the scalar body assembled by set_init_state).
    void set_init_state_lane(ScalarState& state, Real t0, size_t w) const noexcept;

    /// Evaluates a scalar profile at each row's initial theta: lane w of the result is f(theta0_lanes[w]).
    template <typename F>
    Value at_theta0(F&& f) const noexcept {
        return lane::from_lanes<Value>([&](size_t w) { return f(theta0_lanes[w]); });
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivative of Gamma with respect to engine time t.
     * @details Calculates the rate of change of the Lorentz factor based on various physical factors.
     * @param state Current state of the system
     * @param diff Current derivatives
     * @param ad_idx Adiabatic index
     * @return The time derivative of Gamma
     * <!-- ************************************************************************************** -->
     */
    inline Value compute_dGamma_dt(State const& state, State const& diff, Value const& ad_idx, Value const& sin_theta,
                                   Value const& cos_theta) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivative of internal energy with respect to time t.
     * @details Calculates the rate of change of the internal energy density considering adiabatic expansion
     *          and energy injection from newly swept-up material.
     * @param eps_rad radiative efficiency
     * @param state Current state of the system
     * @param diff Current derivatives
     * @param ad_idx Adiabatic index
     * @return The time derivative of internal energy
     * <!-- ************************************************************************************** -->
     */
    inline Value compute_dU_dt(Value const& eps_rad, State const& state, State const& diff, Value const& ad_idx,
                               Value const& sin_theta, Value const& cos_theta) const noexcept;

    Value const dOmega0_;         ///< Initial solid angle element (per lane)
    Real const theta_s_{0};       ///< Critical angle for jet spreading
    Value m_jet0_;                ///< Ejecta mass per solid angle (per lane)
    RadiativeEfficiency eps_rad_; ///< Radiative efficiency evaluator (shared formula in shock-physics.h)
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Updates the forward shock state at a grid point using the current ODE solution.
 * @details Computes physical quantities for the shock state and updates the shock object.
 * @param i Grid index for phi
 * @param j Grid index for theta
 * @param k Grid index for time
 * @param eqn The equation system containing physical parameters
 * @param state Current state of the system
 * @param shock Shock object to update
 * <!-- ************************************************************************************** -->
 */
template <typename Eqn, typename State>
void save_fwd_shock_state(size_t i, size_t j, size_t k, Eqn const& eqn, State const& state, Shock& shock);

/**
 * <!-- ************************************************************************************** -->
 * @brief Integrates a batch of theta rows in lockstep from the phi-wide start `t0` and stores
 *        states on each row's t grid.
 * @details Lanes [0, n_live) map to rows js[0..n_live); trailing lanes are padding (duplicates
 *          of the last live row) and never save output. A lane whose initial Lorentz factor is
 *          below the cutoff is recorded via set_stopping_shock and parked on a live lane.
 * <!-- ************************************************************************************** -->
 */
template <typename FwdEqn>
void grid_solve_fwd_shock(size_t i, std::array<size_t, FwdEqn::lanes> const& js, size_t n_live, Coord const& coord,
                          Shock& shock, FwdEqn const& eqn, Real rtol, Real t0);

/**
 * <!-- ************************************************************************************** -->
 * @brief Generates a forward shock model for given coordinates, medium, and jet parameters.
 * @details Creates a Shock object and solves the shock evolution for each grid point.
 * @param coord Coordinate system definition
 * @param medium The medium through which the shock propagates
 * @param jet The jet (ejecta) driving the shock
 * @param rad_params Radiation parameters
 * @param rtol Relative tolerance for the ODE solver (default: 1e-5)
 * @return A Shock object containing the evolution data
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Shock generate_fwd_shock(Coord const& coord, Medium const& medium, Ejecta const& jet, RadParams const& rad_params,
                         Real rtol = 1e-5);

//========================================================================================================
//                                  template function implementation
//========================================================================================================
#include "../src/dynamics/forward-shock.tpp"
