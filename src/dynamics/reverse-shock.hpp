//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once
#include <array>

#include "../util/macros.h"
/**
 * <!-- ************************************************************************************** -->
 * @struct ReverseState
 * @brief Represents the state variables for the reverse shock simulation.
 * @details It defines a state vector containing properties like shell width, mass, radius, time, and energy.
 *          The struct dynamically adapts its size based on template parameters to include mass/energy
 *          injection capabilities.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta>
struct ReverseState {
    static constexpr bool mass_inject = HasDmdt<Ejecta>;   ///< Whether ejecta has mass injection
    static constexpr bool energy_inject = HasDedt<Ejecta>; ///< Whether ejecta has energy injection
    static constexpr size_t array_size = 12;

    MAKE_THIS_ODEINT_STATE(ReverseState, data, array_size)

    union {
        struct {
            Real Gamma;  ///< Lorentz factor of the shocked region
            Real x4;     ///< Comoving frame width of region 4
            Real x3;     ///< Comoving frame width of region 3
            Real m2;     ///< Shocked medium mass per solid angle
            Real m3;     ///< Shocked ejecta mass per solid angle
            Real U2_th;  ///< internal energy per solid angle in region 3
            Real U3_th;  ///< internal energy per solid angle in region 2
            Real r;      ///< Radius
            Real t_comv; ///< Comoving time
            Real theta;  ///< Angular coordinate theta
            Real eps4;   ///< energy  per solid angle in region 4
            Real m4;     ///< mass per solid angle in region 4
        };
        array_type data;
    };
};

/**
 * <!-- ************************************************************************************** -->
 * @class FRShockEqn
 * @brief Represents the reverse shock (or forward-reverse shock) equation for a given Jet and medium.
 * @details It defines a state vector (an array of 8 Reals) and overloads operator() to compute the
 *          derivatives of the state with respect to radius r.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
class FRShockEqn {
  public:
    using State = ReverseState<Ejecta>;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Constructor for the FRShockEqn class.
     * @details Initializes the forward-reverse shock equation with the given medium, ejecta, and parameters.
     * @param medium The medium through which the shock propagates
     * @param ejecta The ejecta driving the shock
     * @param phi Azimuthal angle
     * @param theta Polar angle
     * @param rad_fwd Radiation params for forward shock
     * @param rad_rvs Radiation params for reverse shock
     * <!-- ************************************************************************************** -->
     */
    FRShockEqn(Medium const& medium, Ejecta const& ejecta, Real phi, Real theta, RadParams const& rad_fwd,
               RadParams const& rad_rvs);

    Medium const& medium;    ///< Reference to the medium properties
    Ejecta const& ejecta;    ///< Reference to the jet properties
    RadParams const rad_fwd; ///< Radiation parameters for forward shock
    RadParams const rad_rvs; ///< Radiation parameters for reverse shock
    Real const phi{0};       ///< Angular coordinate phi
    Real const theta0{0};    ///< Angular coordinate theta
    Real Gamma4{1};          ///< Initial Lorentz factor of the jet
    Real u_x{0};             ///< Reverse shock crossed four-velocity
    Real r_x{0};             ///< Reverse shock crossed radius
    Real B3_ordered_x{0};    ///< Ordered magnetic field in region 3 at crossing
    Real V3_comv_x{0};       ///< Comoving Volume in region 3 at crossing
    Real rho3_x{0};          ///< Density in region 3 at crossing

    /**
     * <!-- ************************************************************************************** -->
     * @brief Implements the reverse shock ODE system.
     * @details Computes the derivatives of state variables with respect to time.
     * @param state Current state of the system
     * @param diff Output derivatives to be populated
     * @param t Current time
     * <!-- ************************************************************************************** -->
     */
    void operator()(State const& state, State& diff, Real t);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Set the initial conditions for the reverse shock ODE.
     * @details Sets up initial state values and determines if the shock has already crossed.
     * @param state State vector to initialize
     * @param t0 Initial time
     * <!-- ************************************************************************************** -->
     */
    void set_init_state(State& state, Real t0) const noexcept;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculates the magnetization parameter of the shell.
     * @details Sigma is defined as (ε/Γmc²) - 1, where ε is the energy per solid angle.
     * @param state Current state of the system
     * @return The magnetization parameter of the shell
     * <!-- ************************************************************************************** -->
     */
    Real compute_shell_sigma(State const& state) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Saves the state of the system at the crossing point.
     * @param state Current state of the system
     * <!-- ************************************************************************************** -->
     */
    void save_cross_state(State const& state);

  private:
    inline Real compute_dGamma_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dU2_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dU3_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dx3_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dx4_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dm2_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dm3_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_deps4_dt(State const& state, State const& diff, Real t) const noexcept;

    inline Real compute_dm4_dt(State const& state, State const& diff, Real t) const noexcept;

    Real deps0_dt{0}; ///< Ejecta energy injection rate
    Real dm0_dt{0};   ///< Ejecta mass injection rate
    Real u4{0};       ///< Four-velocity of the unshocked ejecta
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
