//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <optional>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "pybind.h"

using ArrayDict = std::unordered_map<std::string, xt::xarray<Real>>;
/**
 * @brief Magnetar model class
 *
 * This class represents a magnetar model with a given luminosity and time scale.
 * The operator() function returns the magnetar luminosity as a function of time.
 */
struct PyMagnetar {
    /**
     * @brief Construct a new PyMagnetar object
     *
     * @param L_0 Luminosity at t = t_0 [erg/s]
     * @param t_0 Time at which luminosity is L_0 [s]
     */
    PyMagnetar(Real L_0, Real t_0, Real q = 2) : L_0(L_0), t_0(t_0), q(q) {}

    Real L_0;
    Real t_0;
    Real q;
};

/**
 * @brief Creates a top-hat jet model where energy and Lorentz factor are constant within theta_c
 *
 * @param theta_c Core angle of the jet [radians]
 * @param E_iso Isotropic-equivalent energy [erg]
 * @param Gamma0 Initial Lorentz factor
 * @param spreading Whether to include jet lateral spreading
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model
 * @return Ejecta Configured jet with top-hat profile
 */
Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                   std::optional<PyMagnetar> magnetar = std::nullopt);

/**
 * @brief Creates a Gaussian jet model where energy and Lorentz factor follow Gaussian distribution
 *
 * @param theta_c Core angle of the jet [radians]
 * @param E_iso Isotropic-equivalent energy at the center [erg]
 * @param Gamma0 Initial Lorentz factor at the center
 * @param spreading Whether to include jet lateral spreading
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model
 * @return Ejecta Configured jet with Gaussian profile
 */
Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                     std::optional<PyMagnetar> magnetar = std::nullopt);

/**
 * @brief Creates a power-law jet model where energy and Lorentz factor follow power-law distribution
 *
 * @param theta_c Core angle of the jet [radians]
 * @param E_iso Isotropic-equivalent energy at the center [erg]
 * @param Gamma0 Initial Lorentz factor at the center
 * @param k Power-law index
 * @param spreading Whether to include jet lateral spreading
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model
 * @return Ejecta Configured jet with power-law profile
 */
Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, bool spreading = false, Real duration = 1,
                     std::optional<PyMagnetar> magnetar = std::nullopt);

/**
 * @brief Creates a general power-law jet model 
 *
 * @param theta_c Core angle of the jet for energy [radians]
 * @param theta_cG Core angle of the jet for LF [radians]
 * @param E_iso Isotropic-equivalent energy at the center [erg]
 * @param Gamma0 Initial Lorentz factor at the center
 * @param k Power-law index for energy
 * @param kG Power-law index for LF
 * @param spreading Whether to include jet lateral spreading
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model
 * @return Ejecta Configured jet with power-law profile
 */
Ejecta PyPowerLawGeneralJet(Real theta_c, Real theta_cG, Real E_iso, Real Gamma0, Real k, Real kG, bool spreading = false, Real duration = 1,
                     std::optional<PyMagnetar> magnetar = std::nullopt);
                     
/**
 * @brief Creates a two-component jet model with different properties for narrow and wide components
 * @param theta_n Core angle of the narrow component [radians]
 * @param E_iso_n Isotropic-equivalent energy of the narrow component [erg]
 * @param Gamma0_n Initial Lorentz factor of the narrow component
 * @param theta_w Core angle of the wide component [radians]
 * @param E_iso_w Isotropic-equivalent energy of the wide component [erg]
 * @param Gamma0_w Initial Lorentz factor of the wide component
 * @param spreading Whether to include jet lateral spreading
 * @param duration Engine activity time [seconds]
 * @param magnetar Optional magnetar model
 * @return Ejecta Configured two-component jet with specified properties
 */
Ejecta PyTwoComponentJet(Real theta_n, Real E_iso_n, Real Gamma0_n, Real theta_w, Real E_iso_w, Real Gamma0_w,
                         bool spreading = false, Real duration = 1, std::optional<PyMagnetar> magnetar = std::nullopt);

/**
 * @brief Creates a constant density ISM (Interstellar Medium) environment
 *
 * @param n_ism Number density of the ISM [cm^-3]
 * @return Medium Configured medium with ISM properties
 */
Medium PyISM(Real n_ism);

/**
 * @brief Creates a wind environment with density profile ρ ∝ r^-2
 *
 * @param A_star Wind parameter in units of 5×10^11 g/cm, typical for Wolf-Rayet stars
 * @return Medium Configured medium with wind properties
 */
Medium PyWind(Real A_star);

/**
 * @brief Class representing the observer configuration
 */
class PyObserver {
   public:
    /**
     * @brief Construct observer with given parameters
     *
     * @param lumi_dist Luminosity distance [cm]
     * @param z Redshift
     * @param theta_obs Viewing angle (between jet axis and line of sight) [radians]
     * @param phi_obs Azimuthal angle [radians]
     */
    PyObserver(Real lumi_dist, Real z, Real theta_obs, Real phi_obs = 0)
        : lumi_dist(lumi_dist * unit::cm), z(z), theta_obs(theta_obs), phi_obs(phi_obs) {}

    Real lumi_dist{1e28};  ///< Luminosity distance [internal units]
    Real z{0};             ///< Redshift
    Real theta_obs{0};     ///< Viewing angle [radians]
    Real phi_obs{0};       ///< Azimuthal angle [radians]
};

/**
 * @brief Class representing radiation parameters for synchrotron and IC emission
 */
class PyRadiation {
   public:
    /**
     * @brief Construct radiation model with given microphysical parameters
     *
     * @param eps_e Fraction of shock energy in electrons
     * @param eps_B Fraction of shock energy in magnetic field
     * @param p Electron energy spectral index
     * @param xi_e Fraction of electrons accelerated
     * @param IC_cooling Whether to include IC cooling
     * @param SSC Whether to include SSC (Synchrotron Self-Compton)
     * @param Klein_Nishina Whether to use Klein-Nishina cross-section for IC scattering
     */
    PyRadiation(Real eps_e, Real eps_B, Real p, Real xi_e = 1, bool IC_cooling = false, bool SSC = false,
                bool KN = false)
        : rad(RadParams{eps_e, eps_B, p, xi_e}), IC_cooling(IC_cooling), SSC(SSC), KN(KN) {}

    RadParams rad;
    bool IC_cooling{false};  ///< Whether to include IC cooling
    bool SSC{false};         ///< Whether to include SSC
    bool KN{false};          ///< Whether to include KN
};

/**
 * @brief Convert Ejecta and Medium units to internal code units
 *
 * This function converts the energy, mass, and other parameters of the Ejecta and Medium
 * objects to the internal code units used in the afterglow calculations.
 *
 * @param jet Ejecta object representing jet structure
 * @param medium Medium object representing circumburst environment
 */
void convert_unit(Ejecta& jet, Medium& medium);

/**
 * @brief Main model class for afterglow calculations
 */
class PyModel {
   public:
    /**
     * @brief Construct afterglow model with given components
     *
     * @param jet Ejecta object representing jet structure
     * @param medium Medium object representing circumburst environment
     * @param observer Observer configuration
     * @param fwd_rad Radiation parameters for forward shock
     * @param rvs_rad Optional radiation parameters for reverse shock
     * @param grid_size Resolution of computational grid (phi, theta, time)
     * @param rtol Relative tolerance for numerical calculations
     */
    PyModel(Ejecta jet, Medium medium, PyObserver observer, PyRadiation fwd_rad,
            std::optional<PyRadiation> rvs_rad = std::nullopt,
            std::tuple<Real, Real, Real> resolutions = std::make_tuple(0.3, 2., 5.), Real rtol = 1e-5,
            bool axisymmetric = true)
        : jet(jet),
          medium(medium),
          obs_setup(observer),
          fwd_rad(fwd_rad),
          rvs_rad_opt(rvs_rad),
          phi_resol(std::get<0>(resolutions)),
          theta_resol(std::get<1>(resolutions)),
          t_resol(std::get<2>(resolutions)),
          rtol(rtol),
          axisymmetric(axisymmetric) {
        convert_unit(this->jet, this->medium);
    }

    /**
     * @brief Calculate specific flux at given times and frequencies
     *
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @return FluxDict Dictionary with synchrotron and IC flux components
     */
    ArrayDict specific_flux(PyArray const& t, PyArray const& nu);

    /**
     * @brief Calculate specific flux at given time and frequency (t_i,nu_i) series.
     *
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @return FluxDict Dictionary with synchrotron and IC flux components
     */
    ArrayDict specific_flux_series(PyArray const& t, PyArray const& nu);

    /**
     * @brief Calculate specific flux at given time and frequency (t_i,nu_i) series, sorted by t_i.
     *
     * @param t Observer time array [seconds]
     * @param nu Observer frequency array [Hz]
     * @return FluxDict Dictionary with synchrotron and IC flux components
     */
    ArrayDict specific_flux_sorted_series(PyArray const& t, PyArray const& nu);

    /**
     * @brief Get details of the model configuration
     *
     * @return ArrayDict Dictionary with model details such as shock, electron and photon grids.
     */
    ArrayDict details(PyArray const& t);

   private:
    /**
     * @brief Internal specific flux calculation method using natural units
     *
     * @param t Observer time array [internal units]
     * @param nu Observer frequency array [internal units]
     * @param trace Whether to return the trace of the flux matrix
     * @return FluxDict Dictionary with flux components
     */
    ArrayDict compute_specific_flux(Array const& t, Array const& nu, bool trace = true);

    /**
     * @brief Helper method to calculate flux for a given shock
     *
     * @param shock Forward or reverse shock structure
     * @param coord Coordinate system
     * @param t Observer time array [internal units]
     * @param nu Observer frequency array [internal units]
     * @param obs Observer object
     * @param rad Radiation parameters
     * @param flux_dict Output flux dictionary
     * @param suffix Key suffix for flux components
     * @param return_trace Whether to return the trace of the flux matrix
     */
    void single_shock_emission(Shock const& shock, Coord const& coord, Array const& t, Array const& nu, Observer& obs,
                               PyRadiation rad, ArrayDict& flux_dict, std::string suffix, bool return_trace);

    /**
     * @brief Helper method to calculate details for a given shock
     * @param shock Forward or reverse shock structure
     * @param coord Coordinate system
     * @param t Observer time array [internal units]
     * @param nu Observer frequency array [internal units]
     * @param obs Observer object
     * @param rad Radiation parameters
     * @param detail_dict Output detail dictionary
     * @param suffix Key suffix for detail components
     */
    void single_shock_details(Shock const& shock, Coord const& coord, Array const& t, Observer& obs, PyRadiation rad,
                              ArrayDict& detail_dict, std::string suffix);

    Ejecta jet;                              ///< Jet model
    Medium medium;                           ///< Circumburst medium
    PyObserver obs_setup;                    ///< Observer configuration
    PyRadiation fwd_rad;                     ///< Forward shock radiation parameters
    std::optional<PyRadiation> rvs_rad_opt;  ///< Optional reverse shock radiation parameters
    Real theta_w{con::pi / 2};               ///< Maximum polar angle to calculate
    Real phi_resol{0.3};                     ///< Azimuthal resolution: number of points per degree
    Real theta_resol{2};                     ///< Polar resolution: number of points per degree
    Real t_resol{5};                         ///< Time resolution: number of points per decade
    Real rtol{1e-5};                         ///< Relative tolerance
    bool axisymmetric{true};                 ///< Whether to assume axisymmetric jet
};