//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <cstdio>
#include <functional>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "error_handling.h"
#include "pybind.h"
#include "shock_dispatch.h"
#include "util/macros.h"
#include "util/profiler.h"

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
    PyMagnetar(Real L0, Real t0, Real q = 2) : L0(L0), t0(t0), q(q) {
        AFTERGLOW_REQUIRE_FINITE_POS(L0);
        AFTERGLOW_REQUIRE_FINITE_POS(t0);
        AFTERGLOW_REQUIRE_FINITE_POS(q);
    }

    Real L0; ///< Characteristic luminosity [erg/s]
    Real t0; ///< Spin-down time scale [s]
    Real q;  ///< Power-law index for spin-down

    [[nodiscard]] std::string repr() const {
        char buf[128];
        snprintf(buf, sizeof(buf), "Magnetar(L0=%.6g, t0=%.6g, q=%.6g)", L0, t0, q);
        return buf;
    }
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
JetVariant PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                       std::optional<PyMagnetar> const& magnetar = std::nullopt);

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
JetVariant PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading = false, Real duration = 1,
                         std::optional<PyMagnetar> const& magnetar = std::nullopt);

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
JetVariant PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading = false,
                         Real duration = 1, std::optional<PyMagnetar> const& magnetar = std::nullopt);

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
JetVariant PyPowerLawWing(Real theta_c, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g, bool spreading = false,
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
JetVariant PyStepPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g,
                             bool spreading, Real duration, std::optional<PyMagnetar> const& magnetar);

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
JetVariant PyTwoComponentJet(Real theta_c, Real E_iso, Real Gamma0, Real theta_w, Real E_iso_w, Real Gamma0_w,
                             bool spreading = false, Real duration = 1,
                             std::optional<PyMagnetar> const& magnetar = std::nullopt);

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
ISM PyISM(Real n_ism);

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
 * @param k_m Power-law index for wind density profile (ρ ∝ r^(-k_m), default: 2)
 * @return Medium Configured medium object representing stellar wind properties
 * <!-- ************************************************************************************** -->
 */
MediumVariant PyWind(Real A_star, std::optional<Real> n_ism = std::nullopt, std::optional<Real> n0 = std::nullopt,
                     Real k_m = 2);

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
        : lumi_dist(lumi_dist * unit::cm), z(z), theta_obs(theta_obs), phi_obs(phi_obs) {
        AFTERGLOW_REQUIRE_FINITE_POS(lumi_dist);
        AFTERGLOW_REQUIRE_FINITE_NONNEG(z);
        // theta_obs in [0, pi] (covers off-axis observers up to behind the jet axis)
        AFTERGLOW_REQUIRE(std::isfinite(theta_obs) && theta_obs >= 0 && theta_obs <= con::pi,
                          std::string("theta_obs must be in [0, pi], got ") + std::to_string(theta_obs));
        AFTERGLOW_REQUIRE(std::isfinite(phi_obs),
                          std::string("phi_obs must be finite, got ") + std::to_string(phi_obs));
    }

    Real lumi_dist{1e28}; ///< Luminosity distance [internal units]
    Real z{0};            ///< Redshift
    Real theta_obs{0};    ///< Viewing angle [radians]
    Real phi_obs{0};      ///< Azimuthal angle [radians]

    [[nodiscard]] Real lumi_dist_cgs() const { return lumi_dist / unit::cm; }

    [[nodiscard]] std::string repr() const {
        char buf[256];
        if (phi_obs != 0) {
            snprintf(buf, sizeof(buf), "Observer(lumi_dist=%.6g, z=%.6g, theta_obs=%.6g, phi_obs=%.6g)",
                     lumi_dist / unit::cm, z, theta_obs, phi_obs);
        } else {
            snprintf(buf, sizeof(buf), "Observer(lumi_dist=%.6g, z=%.6g, theta_obs=%.6g)", lumi_dist / unit::cm, z,
                     theta_obs);
        }
        return buf;
    }
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
     * @param ssc Whether to include synchrotron self-Compton emission and IC cooling (default: false)
     * @param kn Whether to include Klein-Nishina corrections for IC processes (default: false)
     * <!-- ************************************************************************************** -->
     */
    PyRadiation(Real eps_e, Real eps_B, Real p, Real xi_e = 1, bool ssc = false, bool kn = false)
        : rad(RadParams{eps_e, eps_B, p, xi_e}), ssc(ssc), kn(kn) {
        AFTERGLOW_REQUIRE_RANGE_OI(eps_e, 0.0, 1.0);
        AFTERGLOW_REQUIRE_RANGE_OI(eps_B, 0.0, 1.0);
        AFTERGLOW_REQUIRE_RANGE_OI(xi_e, 0.0, 1.0);
        // p > 1 is the physical minimum for a finite electron-energy integral. The model handles
        // both slow-cooling (p > 2) and fast-cooling (1 < p < 2) regimes, so don't impose p > 2.
        AFTERGLOW_REQUIRE_GREATER_THAN(p, 1.0);
    }

    RadParams rad;
    bool ssc{false}; ///< Whether to include SSC emission and IC cooling
    bool kn{false};  ///< Whether to include KN

    [[nodiscard]] std::string repr() const {
        char buf[128];
        snprintf(buf, sizeof(buf), "Radiation(eps_e=%.6g, eps_B=%.6g, p=%.6g", rad.eps_e, rad.eps_B, rad.p);
        std::string s = buf;
        if (rad.xi_e != 1) {
            snprintf(buf, sizeof(buf), ", xi_e=%.6g", rad.xi_e);
            s += buf;
        }
        if (ssc) {
            s += ", ssc=True";
        }
        if (kn) {
            s += ", kn=True";
        }
        s += ")";
        return s;
    }
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
void convert_unit_jet(JetVariant& jet);
void convert_unit_medium(MediumVariant& medium);

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

    [[nodiscard]] std::string repr() const {
        std::string s = "Flux(sync";
        if (ssc.dimension() > 0) {
            s += " + ssc";
        }
        if (sync.dimension() > 0) {
            s += ", shape=(";
            for (size_t i = 0; i < sync.dimension(); ++i) {
                if (i > 0) {
                    s += ", ";
                }
                s += std::to_string(sync.shape()[i]);
            }
            s += ")";
        }
        s += ")";
        return s;
    }
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

    [[nodiscard]] std::string repr() const {
        if (total.size() == 0) {
            return "FluxDict(empty)";
        }
        std::string s = "FluxDict(shape=(";
        for (size_t i = 0; i < total.dimension(); ++i) {
            if (i > 0) {
                s += ", ";
            }
            s += std::to_string(total.shape()[i]);
        }
        s += "), components=[fwd.sync";
        if (fwd.ssc.dimension() > 0) {
            s += ", fwd.ssc";
        }
        if (rvs.sync.dimension() > 0) {
            s += ", rvs.sync";
        }
        if (rvs.ssc.dimension() > 0) {
            s += ", rvs.ssc";
        }
        s += "])";
        return s;
    }
};

/**
 * @struct PySkyImage
 * @brief Resolved sky image sequence of the afterglow at given observer times and frequency.
 */
struct PySkyImage {
    XTArray image;              ///< [n_frames, npixel, npixel] surface brightness (erg/cm²/s/Hz/sr)
    std::array<Real, 4> extent; ///< {x_min, x_max, y_min, y_max} angular extent (rad), shared by all frames
    Real pixel_solid_angle{0};  ///< Pixel solid angle (sr), shared by all frames

    [[nodiscard]] std::string repr() const {
        if (image.size() == 0) {
            return "SkyImage(empty)";
        }
        char buf[256];
        snprintf(buf, sizeof(buf),
                 "SkyImage(n_frames=%zu, npixel=%zu, extent=(%.3e, %.3e, %.3e, %.3e) rad, pixel_solid_angle=%.3e sr)",
                 image.shape()[0], image.shape()[1], extent[0], extent[1], extent[2], extent[3], pixel_solid_angle);
        return buf;
    }
};

// Type aliases for IC photon grid
using SynICPhoton = ICPhoton<SynElectrons, SynPhotons>;
using SynICPhotonGrid = xt::xtensor<SynICPhoton, 3>;

/// Callable evaluator for per-cell spectrum queries (takes comoving frequency in Hz)
struct SpectrumEvaluator {
    std::function<Real(Real)> eval_;
    XTArray operator()(PyArray const& nu_comv) const;
};

/// Callable evaluator for per-cell Y(gamma) queries
struct YEvaluator {
    std::function<Real(Real)> eval_;
    XTArray operator()(PyArray const& gamma) const;
};

/// Grid accessor for synchrotron spectrum: sync_spectrum[i,j,k] → SpectrumEvaluator
struct SynSpectrumGrid {
    SynPhotonGrid const* grid_;
};

/// Grid accessor for IC spectrum: ssc_spectrum[i,j,k] → SpectrumEvaluator
struct ICSpectrumGrid {
    SynICPhotonGrid* grid_; // non-const: ICPhoton::compute_I_nu lazily generates spectrum
};

/// Grid accessor for Y(gamma): Y_spectrum[i,j,k] → YEvaluator
struct YSpectrumGrid {
    SynPhotonGrid const* grid_;
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

    // Stored photon grids for per-cell spectrum evaluation
    // Stored photon grids double as presence flags: size() == 0 until assigned.
    SynPhotonGrid syn_photons_;
    SynICPhotonGrid ic_photons_;

    [[nodiscard]] std::string repr() const {
        if (Gamma.size() == 0) {
            return "ShockDetails(empty)";
        }
        std::string s = "ShockDetails(shape=(";
        for (size_t i = 0; i < Gamma.dimension(); ++i) {
            if (i > 0) {
                s += ", ";
            }
            s += std::to_string(Gamma.shape()[i]);
        }
        s += "))";
        return s;
    }
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

    [[nodiscard]] std::string repr() const {
        char buf[128];
        snprintf(buf, sizeof(buf), "SimulationDetails(phi=%zu, theta=%zu, t_src=%zu)", phi.size(), theta.size(),
                 t_src.size());
        return buf;
    }
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Applies inverse-Compton electron cooling (Klein-Nishina or Thomson) when SSC is enabled.
 * <!-- ************************************************************************************** -->
 */
template <typename Electrons, typename Photons>
void apply_ic_cooling(Electrons& syn_e, Photons& syn_ph, Shock const& shock, Coord const& coord,
                      PyRadiation const& rad) {
    if (rad.ssc) {
        if (rad.kn) {
            KN_cooling(syn_e, syn_ph, shock, coord);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock, coord);
        }
    }
}

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
     * @param resolutions Grid resolution tuple (phi_res, theta_res, time_res) in (deg⁻¹, deg⁻¹, decade⁻¹).
     *        The total grid points in each dimension are computed from these resolutions, but
     *        each has a minimum that cannot be reduced further: phi (min 1), theta (min 36), time (min 24).
     *        When omitted, the default is (0.06, 0.15, 6) for forward-shock-only runs and (0.06, 0.2, 10)
     *        with a reverse shock; explicit values are honored as given.
     * @param rtol Relative tolerance for the shock-dynamics ODE solves (default: 1×10⁻⁶; magnetized shells tighten one decade internally)
     * @param axisymmetric Whether to assume axisymmetric jet structure (default: true)
     * <!-- ************************************************************************************** -->
     */
    PyModel(JetVariant jet, MediumVariant medium, PyObserver const& observer, PyRadiation const& fwd_rad,
            std::optional<PyRadiation> const& rvs_rad = std::nullopt,
            std::optional<std::tuple<Real, Real, Real>> const& resolutions = std::nullopt,
            Real rtol = defaults::solver::dynamics_rtol, bool axisymmetric = true, bool radiative_fireball = true)
        : jet_(std::move(jet)),
          medium_(std::move(medium)),
          obs_setup(observer),
          fwd_rad(fwd_rad),
          rvs_rad_opt(rvs_rad),
          rtol(rtol),
          axisymmetric(axisymmetric) {
        // Radiative losses back-react on the blast-wave dynamics by default; radiative_fireball
        // = false selects the adiabatic approximation used by most afterglow codes.
        this->fwd_rad.rad.radiative = radiative_fireball;
        if (rvs_rad_opt) {
            rvs_rad_opt->rad.radiative = radiative_fireball;
        }
        // Forward-shock-only runs default to the coarser calibrated grid;
        // reverse-shock light curves need denser theta and time (see
        // defaults::grid). An explicit `resolutions` is honored as given.
        const auto res = resolutions.value_or(
            rvs_rad ? std::make_tuple(defaults::grid::phi_resolution, defaults::grid::rvs_theta_resolution,
                                      defaults::grid::rvs_time_resolution)
                    : std::make_tuple(defaults::grid::phi_resolution, defaults::grid::theta_resolution,
                                      defaults::grid::time_resolution));
        phi_resol = std::get<0>(res);
        theta_resol = std::get<1>(res);
        t_resol = std::get<2>(res);
        // rtol = 1.0 means "100% relative error is OK" -- the integrator gives up. Reject both ends.
        AFTERGLOW_REQUIRE(std::isfinite(rtol) && rtol > 0 && rtol < 1,
                          std::string("rtol must be in (0, 1), got ") + std::to_string(rtol));
        AFTERGLOW_REQUIRE_FINITE_POS(phi_resol);
        AFTERGLOW_REQUIRE_FINITE_POS(theta_resol);
        AFTERGLOW_REQUIRE_FINITE_POS(t_resol);
        convert_unit_jet(this->jet_);
        convert_unit_medium(this->medium_);
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
     * @brief Compute resolved sky images at multiple observer times and a single frequency.
     * @param t_obs Array of observer times [seconds]
     * @param nu_obs Observer frequency [Hz]
     * @param fov Total field of view [radians] (pixel scale = fov / npixel)
     * @param npixel Number of pixels per side (default 64)
     * @return PySkyImage with (n_frames, npixel, npixel) surface brightness (erg/cm²/s/Hz/sr), shared extent and pixel_solid_angle
     */
    PySkyImage sky_image(PyArray const& t_obs, double nu_obs, double fov, size_t npixel = 128);

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

#ifdef AFTERGLOW_PROFILE
    static auto profile_data() -> std::unordered_map<std::string, double> { return AFTERGLOW_PROFILE_RESULTS(); }
    static auto profile_counters() -> std::unordered_map<std::string, size_t> { return AFTERGLOW_PROFILE_COUNTERS(); }
    static void profile_reset() { AFTERGLOW_PROFILE_RESET(); }
#endif

    // Read-only accessors for Python properties
    [[nodiscard]] PyObserver const& get_observer() const { return obs_setup; }
    [[nodiscard]] PyRadiation const& get_fwd_rad() const { return fwd_rad; }
    [[nodiscard]] std::optional<PyRadiation> const& get_rvs_rad() const { return rvs_rad_opt; }
    [[nodiscard]] std::tuple<Real, Real, Real> get_resolutions() const { return {phi_resol, theta_resol, t_resol}; }

    [[nodiscard]] bool get_radiative_fireball() const { return fwd_rad.rad.radiative; }
    [[nodiscard]] Real get_rtol() const { return rtol; }
    [[nodiscard]] bool get_axisymmetric() const { return axisymmetric; }

    [[nodiscard]] std::string repr() const {
        // Nested objects format themselves so this can never drift from their own reprs.
        std::string s = "Model(observer=" + obs_setup.repr() + ",\n      fwd_rad=" + fwd_rad.repr();
        if (rvs_rad_opt) {
            s += ",\n      rvs_rad=" + rvs_rad_opt->repr();
        }
        char buf[128];
        snprintf(buf, sizeof(buf), ",\n      resolutions=(%.6g, %.6g, %.6g), rtol=%.6g", phi_resol, theta_resol,
                 t_resol, rtol);
        s += buf;
        if (!fwd_rad.rad.radiative) {
            s += ", radiative_fireball=False";
        }
        return s + ")";
    }

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
     * @param obs Observer object for flux calculation. Precondition: obs.observe() has already run
     *            against kinematically identical shock data — the pair solver evolves both shocks
     *            from one ODE state, so a single EAT grid serves forward and reverse emission.
     * @param rad Radiation parameters controlling microphysics
     * @param emission Output flux structure to populate
     * @param flux_func Function to compute flux (either specific_flux or specific_flux_series)
     * <!-- ************************************************************************************** -->
     */
    template <typename Func>
    void single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                               Observer& obs, PyRadiation const& rad, Flux& emission, Func&& flux_func);

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

    JetVariant jet_;                        ///< Jet model (TophatJet, GaussianJet, PowerLawJet, or Ejecta)
    MediumVariant medium_;                  ///< Circumburst medium (ISM, Wind, or generic Medium)
    PyObserver obs_setup;                   ///< Observer configuration
    PyRadiation fwd_rad;                    ///< Forward shock radiation parameters
    std::optional<PyRadiation> rvs_rad_opt; ///< Optional reverse shock radiation parameters
    Real theta_w{con::pi / 2};              ///< Maximum polar angle to calculate
    Real phi_resol;                         ///< Azimuthal resolution: number of points per degree (set in ctor)
    Real theta_resol;                       ///< Polar resolution: number of points per degree (set in ctor)
    Real t_resol;                           ///< Time resolution: number of points per decade (set in ctor)
    Real rtol;                              ///< Relative tolerance for shock dynamics (set in ctor)
    bool axisymmetric;                      ///< Whether to assume axisymmetric jet (set in ctor)
};

//========================================================================================================
//                                  template function implementation
//========================================================================================================

template <typename Func>
void PyModel::single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                    Observer& obs, PyRadiation const& rad, Flux& emission, Func&& flux_func) {
    auto syn_e = [&] {
        AFTERGLOW_PROFILE_SCOPE(syn_electrons);
        return generate_syn_electrons(shock, coord);
    }();

    auto syn_ph = [&] {
        AFTERGLOW_PROFILE_SCOPE(syn_photons);
        return generate_syn_photons(shock, syn_e, coord);
    }();

    {
        AFTERGLOW_PROFILE_SCOPE(cooling);
        apply_ic_cooling(syn_e, syn_ph, shock, coord, rad);
    }

    {
        AFTERGLOW_PROFILE_SCOPE(sync_flux);
        emission.sync = std::invoke(flux_func, obs, t_obs, nu_obs, syn_ph);
    }

    if (rad.ssc) {
        auto IC_ph = [&] {
            AFTERGLOW_PROFILE_SCOPE(ic_photons);
            return generate_IC_photons(syn_e, syn_ph, rad.kn, coord);
        }();
        {
            AFTERGLOW_PROFILE_SCOPE(ssc_flux);
            emission.ssc = std::invoke(flux_func, obs, t_obs, nu_obs, IC_ph);
        }
    }
}

template <typename Func>
auto PyModel::compute_emission(Array const& t_obs, Array const& nu_obs, Func&& flux_func) -> PyFlux {
    AFTERGLOW_PROFILE_SCOPE(total);

    PyFlux flux;
    Observer observer;

    if (!rvs_rad_opt) {
        auto [coord, fwd_shock] = [&] {
            AFTERGLOW_PROFILE_SCOPE(dynamics);
            return solve_fwd_shock(jet_, medium_, t_obs, theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol,
                                   theta_resol, t_resol, axisymmetric, fwd_rad.rad, rtol);
        }();
        {
            AFTERGLOW_PROFILE_SCOPE(EAT_grid);
            observer.observe(coord, fwd_shock, obs_setup.lumi_dist, obs_setup.z);
        }
        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux.fwd,
                              std::forward<Func>(flux_func));
    } else {
        auto [coord, fwd_shock, rvs_shock] = [&] {
            AFTERGLOW_PROFILE_SCOPE(dynamics);
            return solve_shock_pair(jet_, medium_, t_obs, theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol,
                                    theta_resol, t_resol, axisymmetric, fwd_rad.rad, rvs_rad_opt->rad, rtol);
        }();
        // The pair solver writes both shocks from one ODE state: regions 2 and 3
        // ride the same contact discontinuity, so t_comv, r, theta, and Gamma --
        // everything the observer reads -- are identical, and one set of EAT
        // grids serves both emission passes.
        {
            AFTERGLOW_PROFILE_SCOPE(EAT_grid);
            observer.observe(coord, fwd_shock, obs_setup.lumi_dist, obs_setup.z);
        }
        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux.fwd,
                              std::forward<Func>(flux_func));
        single_shock_emission(rvs_shock, coord, t_obs, nu_obs, observer, *rvs_rad_opt, flux.rvs,
                              std::forward<Func>(flux_func));
    }
    return flux;
}
