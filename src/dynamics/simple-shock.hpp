//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "shock.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct SimpleState
 * @brief Represents the state vector for the simple shock equation.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
struct SimpleState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;   ///< whether Ejecta class has dmdt method
    static constexpr bool energy_inject = HasDedt<Ejecta>; ///< whether Ejecta class has dedt method
    /// use least fixed array size for integrator efficiency
    static constexpr size_t array_size = 5 + (mass_inject ? 1 : 0) + (energy_inject ? 1 : 0);

    MAKE_THIS_ODEINT_STATE(SimpleState, data, array_size)

    union {
        struct {
            Real Gamma;  ///< Lorentz factor
            Real m2;     ///< swept mass
            Real r;      ///< radius
            Real t_comv; ///< comoving time
            Real theta;  ///< angle

            // shell energy density per solid angle
            [[no_unique_address]] std::conditional_t<energy_inject, Real, class Empty> eps_jet;

            // shell mass per solid angle
            [[no_unique_address]] std::conditional_t<mass_inject, Real, class Empty> m_jet;
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class SimpleShockEqn
 * @brief Represents the forward shock equation for a given Jet.
 * @details It defines a state vector and overloads operator() to compute the derivatives of the state with
 *          respect to time t. It also declares helper functions for the derivatives. Simple version from
 *          Huang et al. 2000
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
class SimpleShockEqn {
  public:
    using State = SimpleState<Ejecta>;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes a SimpleShockEqn object with medium, ejecta, and other parameters.
     * @details Creates a new shock equation object with references to the medium and ejecta
     *          along with the angular coordinates and energy fraction.
     * @param medium The medium through which the shock propagates
     * @param ejecta The ejecta driving the shock
     * @param phi Azimuthal angle
     * @param theta Polar angle
     * @param rad_params Radiation parameters
     * @param theta_s Critical angle for jet spreading
     * <!-- ************************************************************************************** -->
     */
    SimpleShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, RadParams const& rad_params,
                   Real theta_s);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivatives of the state variables with respect to engine time t.
     * @details Implements the system of ODEs for the simple shock model.
     * @param state Current state of the system
     * @param diff Output derivatives to be populated
     * @param t Current time
     * <!-- ************************************************************************************** -->
     */
    void operator()(State const& state, State& diff, Real t) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Initializes the state vector at time t0.
     * @details Sets appropriate initial values based on ejecta and medium properties.
     * @param state State vector to initialize
     * @param t0 Initial time
     * <!-- ************************************************************************************** -->
     */
    void set_init_state(State& state, Real t0) const noexcept;

    Medium const& medium; ///< Reference to the medium properties
    Ejecta const& ejecta; ///< Reference to the ejecta properties
    Real const phi{0};    ///< Angular coordinate phi
    Real const theta0{0}; ///< Angular coordinate theta
    RadParams const rad;  ///< Radiation parameters

  private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Computes the derivative of Gamma with respect to engine time t.
     * @details Calculates the rate of change of the Lorentz factor based on swept-up mass and energy injection.
     * @param eps_rad radiative efficiency
     * @param state Current state of the system
     * @param diff Current derivatives
     * @return The time derivative of Gamma
     * <!-- ************************************************************************************** -->
     */
    Real dGamma_dt(Real eps_rad, State const& state, State const& diff) const noexcept;

    Real const dOmega0{0}; ///< Initial solid angle
    Real const theta_s{0}; ///< Critical angle for jet spreading
    Real m_jet0{0};        ///< Ejecta mass per solid angle
};

#include "../src/dynamics/simple-shock.tpp"
