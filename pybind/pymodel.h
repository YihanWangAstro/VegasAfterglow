//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <optional>
#include <utility>
#include <vector>

#include "../include/afterglow.h"
#include "pybind.h"
#include "util/macros.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct PyMagnetar
 * @brief Magnetar model for gamma-ray burst central engine energy injection.
 * @details This class represents a magnetar spin-down model with characteristic luminosity, spin-down
 *          time scale, and power-law index. The magnetar provides continuous energy injection into
 *          the gamma-ray burst ejecta, following the relation: L(t) = L0 * (1 + t/t0)^(-q).
 *          This model is commonly used to explain X-ray plateaus and late-time energy injection
 *          in gamma-ray bursts.
 * <!-- ************************************************************************************** -->
 */
struct PyMagnetar {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct a new PyMagnetar object with spin-down parameters.
     * @details Initializes the magnetar model with the characteristic luminosity at the reference
     *          time, the spin-down time scale, and the power-law index for the energy injection.
     * @param L0 Characteristic luminosity at reference time t = t0 [erg/s]
     * @param t0 Spin-down time scale, time at which luminosity drops to L0/2^q [s]
     * @param q Power-law index for spin-down evolution (typically 1-3, default 2)
     * <!-- ************************************************************************************** -->
     */
    PyMagnetar(Real L0, Real t0, Real q = 2) : L0(L0), t0(t0), q(q) {}

    Real L0; ///< Characteristic luminosity [erg/s]
    Real t0; ///< Spin-down time scale [s]
    Real q;  ///< Power-law index for spin-down
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a top-hat jet model with uniform energy and Lorentz factor distributions.
 * @details Generates an ejecta structure representing a top-hat jet where the energy density and
 *          initial Lorentz factor are constant within the core angle theta_c and zero outside.
 *          This is the simplest jet model and is often used as a baseline for comparison with
 *          more complex structured jets. The model can optionally include lateral spreading
 *          and magnetar energy injection.
 * @param theta_c Core opening angle of the jet [radians]
 * @param E_iso Isotropic-equivalent energy of the entire jet [erg]
 * @param Gamma0 Initial bulk Lorentz factor of the ejecta
 * @param spreading Whether to include lateral jet spreading during expansion (default: false)
 * @param duration Duration of central engine activity [seconds] (default: 1)
 * @param magnetar Optional magnetar model for continuous energy injection (default: none)
 * @return Ejecta Configured ejecta object representing the top-hat jet structure
 * <!-- ************************************************************************************** -->
 */
Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                   const std::optional<PyMagnetar>& magnetar = std::nullopt);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a Gaussian jet model with smooth energy and Lorentz factor profiles.
 * @details Generates an ejecta structure representing a Gaussian-structured jet where both the
 *          energy density and initial Lorentz factor decrease smoothly from the jet axis
 *          following a Gaussian profile. This model provides a more realistic smooth transition
 *          compared to the sharp cutoff of top-hat jets and is often used to model structured
 *          jets in gamma-ray bursts.
 * @param theta_c Characteristic opening angle (1/e width) of the Gaussian jet [radians]
 * @param E_iso Isotropic-equivalent energy at the jet center [erg]
 * @param Gamma0 Initial bulk Lorentz factor at the jet center
 * @param spreading Whether to include lateral jet spreading during expansion (default: false)
 * @param duration Duration of central engine activity [seconds] (default: 1)
 * @param magnetar Optional magnetar model for continuous energy injection (default: none)
 * @return Ejecta Configured ejecta object representing the Gaussian jet structure
 * <!-- ************************************************************************************** -->
 */
Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                     const std::optional<PyMagnetar>& magnetar = std::nullopt);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a power-law structured jet with customizable angular energy distribution.
 * @details Generates an ejecta structure representing a power-law jet where the energy density
 *          and initial Lorentz factor decrease as power laws from the jet center. This model
 *          allows independent control of the energy and Lorentz factor profiles through
 *          separate power-law indices, providing flexibility in modeling diverse jet structures
 *          observed in gamma-ray bursts.
 * @param theta_c Core opening angle defining the transition region [radians]
 * @param E_iso Isotropic-equivalent energy at the jet center [erg]
 * @param Gamma0 Initial bulk Lorentz factor at the jet center
 * @param k_e Power-law index for energy angular distribution (E ∝ θ^(-k_e))
 * @param k_g Power-law index for Lorentz factor angular distribution (Γ ∝ θ^(-k_g))
 * @param spreading Whether to include lateral jet spreading during expansion (default: false)
 * @param duration Duration of central engine activity [seconds] (default: 1)
 * @param magnetar Optional magnetar model for continuous energy injection (default: none)
 * @return Ejecta Configured ejecta object representing the power-law jet structure
 * <!-- ************************************************************************************** -->
 */
Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading = false,
                     Real duration = 1, const std::optional<PyMagnetar>& magnetar = std::nullopt);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a power-law wing jet model with continuous transition from core to wing.
 * @details Generates an ejecta structure representing a jet with a narrow uniform core
 *          and power-law wings that extend smoothly without discontinuities at the core
 *          boundary. This model combines the efficiency of a narrow core with extended
 *          wings, providing a realistic representation of jets with continuous angular
 *          structure transitions.
 * @param theta_c Core opening angle where transition to power-law begins [radians]
 * @param E_iso_w Isotropic-equivalent energy at the theta_c [erg]
 * @param Gamma0_w Initial bulk Lorentz factor at the theta_c
 * @param k_e Power-law index for energy in the wing region (E ∝ θ^(-k_e))
 * @param k_g Power-law index for Lorentz factor in the wing region (Γ ∝ θ^(-k_g))
 * @param spreading Whether to include lateral jet spreading during expansion (default: false)
 * @param duration Duration of central engine activity [seconds] (default: 1)
 * @return Ejecta Configured ejecta object representing the power-law wing jet structure
 * <!-- ************************************************************************************** -->
 */
Ejecta PyPowerLawWing(Real theta_c, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g, bool spreading = false,
                      Real duration = 1);
/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a step power-law jet model: a pencil beam core and powerlaw wing with jump at theta_c.
 * @details Generates an ejecta structure with a narrow uniform core and power-law wings,
 *          with a discontinuous jump in properties at the core boundary. This model represents
 *          jets with distinct core and wing components, useful for modeling jets with clear
 *          transitions between high-energy cores and extended lower-energy wings.
 * @param theta_c Core angle of the jet [radians]
 * @param E_iso Isotropic-equivalent energy of the core region [erg]
 * @param Gamma0 Initial Lorentz factor of the core region
 * @param E_iso_w Isotropic-equivalent energy of the wing region at theta_c [erg]
 * @param Gamma0_w Initial Lorentz factor of the wing region at theta_c
 * @param k_e Power-law index for energy in wing region (E ∝ θ^(-k_e))
 * @param k_g Power-law index for Lorentz factor in wing region (Γ ∝ θ^(-k_g))
 * @param spreading Whether to include jet lateral spreading during expansion
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model for continuous energy injection
 * @return Ejecta Configured jet with step power-law profile
 * <!-- ************************************************************************************** -->
 */
Ejecta PyStepPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g,
                         bool spreading, Real duration, const std::optional<PyMagnetar>& magnetar);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a two-component jet model with different properties for narrow and wide components.
 * @details Generates an ejecta structure representing a jet with two distinct components:
 *          a narrow, high-energy core and a wider, lower-energy wing. This model is useful
 *          for representing jets with complex angular structure where different emission
 *          regions have fundamentally different properties and evolution.
 * @param theta_c Core angle of the narrow component [radians]
 * @param E_iso Isotropic-equivalent energy of the narrow component [erg]
 * @param Gamma0 Initial Lorentz factor of the narrow component
 * @param theta_w Core angle of the wide component [radians]
 * @param E_iso_w Isotropic-equivalent energy of the wide component [erg]
 * @param Gamma0_w Initial Lorentz factor of the wide component
 * @param spreading Whether to include jet lateral spreading during expansion
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model for continuous energy injection
 * @return Ejecta Configured two-component jet with specified properties
 * <!-- ************************************************************************************** -->
 */
Ejecta PyTwoComponentJet(Real theta_c, Real E_iso, Real Gamma0, Real theta_w, Real E_iso_w, Real Gamma0_w,
                         bool spreading = false, Real duration = 1,
                         const std::optional<PyMagnetar>& magnetar = std::nullopt);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a constant density ISM (Interstellar Medium) environment.
 * @details Generates a medium structure representing a uniform interstellar medium
 *          with constant number density. This is the simplest external medium model
 *          and is appropriate for gamma-ray bursts occurring in relatively uniform
 *          environments away from their progenitor stars.
 * @param n_ism Number density of the ISM [cm^-3]
 * @return Medium Configured medium with constant ISM density properties
 * <!-- ************************************************************************************** -->
 */
Medium PyISM(Real n_ism);

/**
 * <!-- ************************************************************************************** -->
 * @brief Creates a stellar wind environment with r^(-2) density profile.
 * @details Generates a medium structure representing a stellar wind from a massive star
 *          progenitor, with density following ρ(r) = A* /(4πr²) where A_* is the wind
 *          parameter. This model is appropriate for gamma-ray bursts occurring in the
 *          vicinity of massive stars such as Wolf-Rayet stars or other evolved massive
 *          stars with strong stellar winds.
 * @param A_star Wind parameter A* = Ṁv_w/(4π) in units of 5×10^11 g/cm [dimensionless]
 * @param n_ism Constant density component for large radii [cm^-3] (default: 0)
 * @param n0 Inner boundary density for stratified wind model [cm^-3] (default: ∞)
 * @param k Power-law index for wind density profile (ρ ∝ r^(-k), default: 2)
 * @return Medium Configured medium object representing stellar wind properties
 * <!-- ************************************************************************************** -->
 */
Medium PyWind(Real A_star, Real n_ism = 0, Real n0 = con::inf, Real k = 2);

/**
 * <!-- ************************************************************************************** -->
 * @class PyObserver
 * @brief Observer configuration for gamma-ray burst afterglow observations.
 * @details This class encapsulates all observational parameters needed to calculate
 *          the observed properties of gamma-ray burst afterglows. It includes the
 *          cosmological distance, redshift, and viewing geometry that determine how
 *          the intrinsic emission appears to the observer.
 * <!-- ************************************************************************************** -->
 */
class PyObserver {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct observer with given cosmological and geometric parameters.
     * @details Initializes the observer configuration with distance, redshift, and viewing
     *          angles. The viewing angle theta_obs is crucial for determining relativistic
     *          beaming effects and the observed light curve evolution.
     * @param lumi_dist Luminosity distance to the gamma-ray burst [cm]
     * @param z Cosmological redshift affecting observed frequencies and times
     * @param theta_obs Viewing angle between jet axis and line of sight [radians]
     * @param phi_obs Azimuthal viewing angle for non-axisymmetric jets [radians]
     * <!-- ************************************************************************************** -->
     */
    PyObserver(Real lumi_dist, Real z, Real theta_obs, Real phi_obs = 0)
        : lumi_dist(lumi_dist * unit::cm), z(z), theta_obs(theta_obs), phi_obs(phi_obs) {}

    Real lumi_dist{1e28}; ///< Luminosity distance [internal units]
    Real z{0};            ///< Redshift
    Real theta_obs{0};    ///< Viewing angle [radians]
    Real phi_obs{0};      ///< Azimuthal angle [radians]
};

/**
 * <!-- ************************************************************************************** -->
 * @class PyRadiation
 * @brief Radiation parameters for synchrotron and inverse Compton emission processes.
 * @details This class encapsulates the microphysical parameters that govern particle
 *          acceleration and radiation processes in relativistic shocks. It controls
 *          the efficiency of converting shock energy into accelerated electrons and
 *          magnetic fields, as well as enabling advanced radiation processes.
 * <!-- ************************************************************************************** -->
 */
class PyRadiation {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct radiation model with given microphysical parameters.
     * @details Initializes the radiation physics model with parameters controlling particle
     *          acceleration efficiency, magnetic field generation, and advanced radiation
     *          processes. These parameters are fundamental to afterglow modeling as they
     *          determine the spectral shape and overall normalization of the emission.
     * @param eps_e Fraction of shock energy transferred to relativistic electrons
     * @param eps_B Fraction of shock energy stored in magnetic field
     * @param p Electron energy spectral index (typically 2.2-2.8)
     * @param xi_e Fraction of shock-heated electrons that are accelerated to relativistic energies
     * @param ssc_cooling Whether to include inverse Compton cooling of electrons
     * @param ssc Whether to include synchrotron self-Compton emission (default: false)
     * @param kn Whether to include Klein-Nishina corrections for IC processes (default: false)
     * <!-- ************************************************************************************** -->
     */
    PyRadiation(Real eps_e, Real eps_B, Real p, Real xi_e = 1, bool ssc_cooling = false, bool ssc = false,
                bool kn = false)
        : rad(RadParams{eps_e, eps_B, p, xi_e}), ssc_cooling(ssc_cooling), ssc(ssc), kn(kn) {}

    RadParams rad;
    bool ssc_cooling{false}; ///< Whether to include IC cooling
    bool ssc{false};         ///< Whether to include SSC
    bool kn{false};          ///< Whether to include KN
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Convert Ejecta and Medium units to internal code units.
 * @details This function performs unit conversion from physical CGS units to the internal
 *          dimensionless units used throughout the afterglow calculations. This conversion
 *          is essential for numerical stability and computational efficiency in the relativistic
 *          hydrodynamics and radiation transfer calculations.
 * @param jet Ejecta object representing jet structure and energy distribution
 * @param medium Medium object representing circumburst environment density profile
 * <!-- ************************************************************************************** -->
 */
void convert_unit(Ejecta& jet, Medium& medium);

using XTArray = xt::xarray<Real>;

/**
 * <!-- ************************************************************************************** -->
 * @struct Flux
 * @brief Container for synchrotron and synchrotron self-Compton flux components.
 * @details This structure organizes the different radiation components calculated
 *          for each shock region, allowing separate tracking of synchrotron emission
 *          and inverse Compton scattering contributions to the total observed flux.
 * <!-- ************************************************************************************** -->
 */
struct Flux {
    XTArray sync; ///< Synchrotron emission flux [mJy]
    XTArray ssc;  ///< Synchrotron self-Compton flux [mJy]
};

/**
 * <!-- ************************************************************************************** -->
 * @struct PyFlux
 * @brief Complete flux structure combining forward and reverse shock contributions.
 * @details This structure aggregates flux from all emission processes and shock regions,
 *          providing both the total observed flux and the individual components for
 *          detailed analysis of the afterglow emission mechanisms.
 * <!-- ************************************************************************************** -->
 */
struct PyFlux {
    XTArray total; ///< Total flux from all components [mJy]
    Flux fwd;      ///< Forward shock emission components
    Flux rvs;      ///< Reverse shock emission components

    ///< Calculate total flux by summing all components
    void calc_total();
};

/**
 * <!-- ************************************************************************************** -->
 * @struct PyShock
 * @brief Comprehensive shock evolution data for detailed afterglow analysis.
 * @details This structure contains all the physical quantities that evolve during
 *          shock propagation, including dynamics, magnetic fields, particle distributions,
 *          and characteristic frequencies. It provides complete information for understanding
 *          the physical processes driving the afterglow emission.
 * <!-- ************************************************************************************** -->
 */
struct PyShock {
    XTArray t_comv;   ///< Comoving time [s]
    XTArray t_obs;    ///< Observer time [s]
    XTArray Gamma;    ///< Bulk Lorentz factor
    XTArray Gamma_th; ///< Thermal Lorentz factor
    XTArray B_comv;   ///< Comoving magnetic field [G]
    XTArray r;        ///< Shock radius [cm]
    XTArray theta;    ///< Polar angle [rad]
    XTArray N_p;      ///< Total proton number in shock
    XTArray N_e;      ///< Total electron number in shock
    XTArray gamma_m;  ///< Minimum electron Lorentz factor
    XTArray gamma_c;  ///< Cooling electron Lorentz factor
    XTArray gamma_M;  ///< Maximum electron Lorentz factor
    XTArray gamma_a;  ///< Absorption electron Lorentz factor
    XTArray gamma_m_hat;
    XTArray gamma_c_hat;
    XTArray nu_m; ///< Synchrotron frequency for γ_m [Hz]
    XTArray nu_c; ///< Synchrotron frequency for γ_c [Hz]
    XTArray nu_M; ///< Synchrotron frequency for γ_M [Hz]
    XTArray nu_a; ///< Synchrotron frequency for γ_a [Hz]
    XTArray nu_m_hat;
    XTArray nu_c_hat;
    XTArray Y_T;
    XTArray I_nu_max; ///< Maximum specific intensity [erg/s/Hz]
    XTArray Doppler;  ///< Doppler factor for beaming
};

/**
 * <!-- ************************************************************************************** -->
 * @struct PyDetails
 * @brief Complete simulation details including coordinate system and shock evolution.
 * @details This structure provides comprehensive information about the simulation grid
 *          and the evolution of both forward and reverse shocks. It enables detailed
 *          analysis of the physical processes and validation of the numerical results.
 * <!-- ************************************************************************************** -->
 */
struct PyDetails {
    XTArray phi;   ///< Azimuthal coordinate array [rad]
    XTArray theta; ///< Polar coordinate array [rad]
    XTArray t_src; ///< Source frame time array [s]

    PyShock fwd; ///< Forward shock evolution details
    PyShock rvs; ///< Reverse shock evolution details
};

/**
 * <!-- ************************************************************************************** -->
 * @class PyModel
 * @brief Main gamma-ray burst afterglow model for comprehensive multi-wavelength predictions.
 * @details This class provides the primary interface for gamma-ray burst afterglow modeling,
 *          combining jet structure, circumburst medium, observational setup, and radiation
 *          physics to compute multi-wavelength light curves, spectra, and internal shock
 *          evolution. It handles both forward and optional reverse shock emission with
 *          advanced radiation processes including synchrotron emission, synchrotron self-
 *          absorption, and inverse Compton scattering.
 * <!-- ************************************************************************************** -->
 */
class PyModel {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct complete afterglow model with physical and numerical parameters.
     * @details Initializes the afterglow model by combining all physical components and
     *          setting numerical parameters for the simulation. The model can handle both
     *          axisymmetric and non-axisymmetric jet structures with customizable resolution.
     * @param jet Ejecta object defining the jet structure and energy distribution
     * @param medium Medium object representing the circumburst environment density profile
     * @param observer Observer configuration including distance, redshift, and viewing angles
     * @param fwd_rad Radiation parameters for forward shock microphysics
     * @param rvs_rad Optional radiation parameters for reverse shock (default: none)
     * @param resolutions Grid resolution tuple (phi_res, theta_res, time_res) in (deg⁻¹, deg⁻¹, decade⁻¹)
     * @param rtol Relative tolerance for numerical integration (default: 1×10⁻⁵)
     * @param axisymmetric Whether to assume axisymmetric jet structure (default: true)
     * <!-- ************************************************************************************** -->
     */
    PyModel(Ejecta jet, Medium medium, const PyObserver& observer, const PyRadiation& fwd_rad,
            const std::optional<PyRadiation>& rvs_rad = std::nullopt,
            const std::tuple<Real, Real, Real>& resolutions = std::make_tuple(0.3, 1, 10), Real rtol = 1e-5,
            bool axisymmetric = true)
        : jet_(std::move(jet)),
          medium_(std::move(medium)),
          obs_setup(observer),
          fwd_rad(fwd_rad),
          rvs_rad_opt(rvs_rad),
          phi_resol(std::get<0>(resolutions)),
          theta_resol(std::get<1>(resolutions)),
          t_resol(std::get<2>(resolutions)),
          rtol(rtol),
          axisymmetric(axisymmetric) {
        convert_unit(this->jet_, this->medium_);
    }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate specific flux at given times and frequencies.
     * @details Computes the observed flux density at all combinations of specified times
     *          and frequencies, creating a complete flux grid. This is the primary method
     *          for generating light curves and spectra from the afterglow model.
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @return PyFlux Structure with synchrotron and inverse Compton flux components
     * <!-- ************************************************************************************** -->
     */
    PyFlux flux_density_grid(PyArray const& t, PyArray const& nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate specific flux at given time and frequency (t_i,nu_i) series.
     * @details Computes the observed flux density at corresponding pairs of times and
     *          frequencies (t[i], nu[i]), useful for modeling observational data where
     *          each measurement has a specific time-frequency combination.
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @return PyFlux Structure with synchrotron and inverse Compton flux components
     * <!-- ************************************************************************************** -->
     */
    PyFlux flux_density(PyArray const& t, PyArray const& nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate bolometric flux integrated over a frequency range.
     * @details Computes the total observed flux by integrating the specific flux
     *          over a specified frequency range [nu_min, nu_max]. This method provides
     *          bolometric light curves that capture the overall energy output of the
     *          afterglow across the electromagnetic spectrum.
     * @param t Observer time array [seconds]
     * @param nu_min Minimum frequency for integration [Hz]
     * @param nu_max Maximum frequency for integration [Hz]
     * @param num_nu Number of frequency sampling points for numerical integration
     * @return PyFlux Structure with synchrotron and inverse Compton flux components
     * <!-- ************************************************************************************** -->
     */
    PyFlux flux(PyArray const& t, double nu_min, double nu_max, size_t num_nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate specific flux with finite exposure time averaging.
     * @details Computes the observed flux density accounting for finite exposure times
     *          by averaging over multiple time points within each exposure window.
     *          This method provides realistic modeling of observational constraints
     *          where the measured flux represents an average over the observation period.
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @param expo_time Exposure time array [seconds]
     * @param num_points Number of sampling points within each exposure time window
     * @return PyFlux Structure with synchrotron and inverse Compton flux components
     * <!-- ************************************************************************************** -->
     */
    PyFlux flux_density_exposures(PyArray const& t, PyArray const& nu, PyArray const& expo_time,
                                  size_t num_points = 10);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Get comprehensive details of the shock evolution and model configuration.
     * @details Provides detailed information about the shock dynamics, particle distributions,
     *          magnetic fields, and characteristic frequencies throughout the simulation.
     *          This method is essential for understanding the physical processes and validating
     *          the model behavior.
     * @param t_min Minimum observer time for detailed output [seconds]
     * @param t_max Maximum observer time for detailed output [seconds]
     * @return PyDetails Structure with comprehensive shock evolution details
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] PyDetails details(Real t_min, Real t_max) const;

    [[nodiscard]] Array medium(Real phi, Real theta, Array const& r) const;

    [[nodiscard]] Array jet_E_iso(Real phi, Array const& theta) const;

    [[nodiscard]] Array jet_Gamma0(Real phi, Array const& theta) const;

  private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Internal emission calculation method using natural units.
     * @details Template method that handles the core emission calculation logic using internal
     *          dimensionless units. This method sets up the coordinate system, generates
     *          shocks, and delegates to the appropriate emission calculation function.
     * @param t_obs Observer time array [internal units]
     * @param nu_obs Observer frequency array [internal units]
     * @param flux_func Function to compute flux (either specific_flux or specific_flux_series)
     * @return PyFlux Structure with flux components
     * <!-- ************************************************************************************** -->
     */
    template <typename Func>
    PyFlux compute_emission(Array const& t_obs, Array const& nu_obs, Func&& flux_func);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Helper method to calculate emission for a given shock region.
     * @details Template method that computes emission from a single shock (forward or reverse)
     *          including synchrotron radiation and optionally synchrotron self-Compton and
     *          inverse Compton cooling effects. This method handles the radiation physics
     *          for a complete shock region.
     * @param shock Forward or reverse shock structure
     * @param coord Coordinate system for the simulation
     * @param t_obs Observer time array [internal units]
     * @param nu_obs Observer frequency array [internal units]
     * @param obs Observer object for flux calculation
     * @param rad Radiation parameters controlling microphysics
     * @param emission Output flux structure to populate
     * @param flux_func Function to compute flux (either specific_flux or specific_flux_series)
     * <!-- ************************************************************************************** -->
     */
    template <typename Func>
    void single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                               Observer& obs, PyRadiation rad, Flux& emission, Func&& flux_func);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Helper method to extract detailed evolution information for a given shock.
     * @details Collects comprehensive physical quantities from a shock region including
     *          dynamics, magnetic fields, particle distributions, and characteristic
     *          frequencies. This information is essential for understanding the physical
     *          processes driving the afterglow emission.
     * @param shock Forward or reverse shock structure
     * @param coord Coordinate system for the simulation
     * @param obs Observer object for coordinate transformations
     * @param rad Radiation parameters for particle calculations
     * @param details Output structure to populate with shock evolution data
     * <!-- ************************************************************************************** -->
     */
    void single_evo_details(Shock const& shock, Coord const& coord, Observer& obs, PyRadiation const& rad,
                            PyShock& details) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Generate time-frequency sampling points with exposure time averaging.
     * @details Creates a sampling grid that accounts for finite exposure times by generating
     *          multiple time points within each exposure window. This method ensures accurate
     *          flux averaging for realistic observational scenarios.
     * <!-- ************************************************************************************** -->
     */
    struct ExposureSampling {
        Array t_obs_sorted;
        Array nu_obs_sorted;
        std::vector<size_t> idx_sorted;
    };

    static ExposureSampling generate_exposure_sampling(PyArray const& t, PyArray const& nu, PyArray const& expo_time,
                                                       size_t num_points);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Average flux results over exposure time samples.
     * @details Processes the detailed flux calculations from multiple time samples within
     *          each exposure window to produce the final averaged flux values that would
     *          be observed with finite exposure times.
     * <!-- ************************************************************************************** -->
     */
    static void average_exposure_flux(PyFlux& result, std::vector<size_t> const& idx_sorted, size_t original_size,
                                      size_t num_points);

    Ejecta jet_;                            ///< Jet model
    Medium medium_;                         ///< Circumburst medium
    PyObserver obs_setup;                   ///< Observer configuration
    PyRadiation fwd_rad;                    ///< Forward shock radiation parameters
    std::optional<PyRadiation> rvs_rad_opt; ///< Optional reverse shock radiation parameters
    Real theta_w{con::pi / 2};              ///< Maximum polar angle to calculate
    Real phi_resol{0.3};                    ///< Azimuthal resolution: number of points per degree
    Real theta_resol{1};                    ///< Polar resolution: number of points per degree
    Real t_resol{10};                       ///< Time resolution: number of points per decade
    Real rtol{1e-5};                        ///< Relative tolerance
    bool axisymmetric{true};                ///< Whether to assume axisymmetric jet
};

//========================================================================================================
//                                  template function implementation
//========================================================================================================

template <typename Func>
void PyModel::single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                    Observer& obs, PyRadiation rad, Flux& emission, Func&& flux_func) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.ssc_cooling) {
        if (rad.kn) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }

    emission.sync = std::invoke(flux_func, obs, t_obs, nu_obs, syn_ph);

    if (rad.ssc) {
        auto IC_ph = generate_IC_photons(syn_e, syn_ph, rad.kn);
        emission.ssc = std::invoke(flux_func, obs, t_obs, nu_obs, IC_ph);
    }
}

template <typename Func>
auto PyModel::compute_emission(Array const& t_obs, Array const& nu_obs, Func&& flux_func) -> PyFlux {
    Coord coord = auto_grid(jet_, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    PyFlux flux;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium_, jet_, fwd_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux.fwd,
                              std::forward<Func>(flux_func));

        return flux;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium_, jet_, fwd_rad.rad, rvs_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux.fwd,
                              std::forward<Func>(flux_func));

        single_shock_emission(rvs_shock, coord, t_obs, nu_obs, observer, rvs_rad, flux.rvs,
                              std::forward<Func>(flux_func));

        return flux;
    }
}
