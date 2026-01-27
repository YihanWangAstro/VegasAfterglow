//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "../include/afterglow.h"
#include "core/mesh.h"
#include "pybind.h"
#include "util/macros.h"
#include "util/utilities.h"

/**
 * <!-- ************************************************************************************** -->
 * @struct FluxData
 * @brief Container for observational flux data and model predictions for chi-squared analysis.
 * @details This structure organizes time-series flux measurements with their uncertainties,
 *          model predictions, and statistical weights. It provides the fundamental data
 *          structure for parameter estimation and goodness-of-fit evaluation in afterglow
 *          modeling. The structure supports both light curves and spectral measurements.
 * <!-- ************************************************************************************** -->
 */
struct FluxData {
    Array t;        ///< Observation times [seconds]
    Array nu;       ///< Observation frequencies [Hz]
    Array Fv_obs;   ///< Observed flux densities [mJy]
    Array Fv_err;   ///< Observational uncertainties [mJy]
    Array Fv_model; ///< Model-predicted flux densities [mJy]
    Array weights;  ///< Statistical weights for each data point

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate chi-squared goodness-of-fit statistic for this data set.
     * @details Computes the weighted chi-squared statistic comparing observed and model flux
     *          densities, accounting for observational uncertainties and optional weights.
     *          This metric quantifies how well the theoretical model reproduces the observations.
     * @return double Chi-squared value for goodness-of-fit assessment
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] double estimate_chi2() const;
};

/**
 * <!-- ************************************************************************************** -->
 * @struct MultiBandData
 * @brief Comprehensive multi-wavelength observational dataset for gamma-ray burst afterglow fitting.
 * @details This structure manages heterogeneous observational data including light curves
 *          at different frequencies, broadband spectra at specific times, and combined
 *          multi-wavelength datasets. It provides methods for data ingestion, organization,
 *          and statistical analysis, serving as the foundation for MCMC parameter estimation
 *          and model comparison in afterglow studies.
 * <!-- ************************************************************************************** -->
 */
struct MultiBandData {
    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate total chi-squared statistic across all observational datasets.
     * @details Computes the cumulative chi-squared value by summing contributions from all
     *          individual FluxData objects, providing a global measure of model fit quality
     *          across the entire multi-wavelength dataset.
     * @return double Total chi-squared value for the complete observational dataset
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] double estimate_chi2() const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Add a single-frequency light curve to the observational dataset.
     * @details Incorporates time-series flux measurements at a fixed observing frequency,
     *          typically used for adding radio, optical, or X-ray light curves to the
     *          multi-wavelength dataset for comprehensive afterglow modeling.
     * @param nu Observing frequency [Hz]
     * @param t Array of observation times [seconds]
     * @param Fv_obs Array of observed flux densities [mJy]
     * @param Fv_err Array of observational uncertainties [mJy]
     * @param weights Optional statistical weights for each measurement
     * <!-- ************************************************************************************** -->
     */
    void add_flux_density(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err,
                          std::optional<PyArray> const& weights = std::nullopt);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Add frequency-integrated flux measurements over a specified frequency band.
     * @details Incorporates broadband flux measurements integrated over a frequency range,
     *          useful for modeling bolometric or filter-integrated observations where the
     *          effective frequency depends on the spectral shape.
     * @param nu_min Lower frequency bound of the observing band [Hz]
     * @param nu_max Upper frequency bound of the observing band [Hz]
     * @param num_points Number of frequency sampling points for integration
     * @param t Array of observation times [seconds]
     * @param Fv_obs Array of observed integrated flux densities [mJy]
     * @param Fv_err Array of observational uncertainties [mJy]
     * @param weights Optional statistical weights for each measurement
     * <!-- ************************************************************************************** -->
     */
    void add_flux(double nu_min, double nu_max, size_t num_points, PyArray const& t, PyArray const& Fv_obs,
                  PyArray const& Fv_err, const std::optional<PyArray>& weights = std::nullopt);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Add a broadband spectrum at a specific observation time.
     * @details Incorporates simultaneous multi-frequency flux measurements at a fixed time,
     *          providing spectral information crucial for constraining radiation physics
     *          parameters and distinguishing between different emission processes.
     * @param t Observation time [seconds]
     * @param nu Array of observing frequencies [Hz]
     * @param Fv_obs Array of observed flux densities [mJy]
     * @param Fv_err Array of observational uncertainties [mJy]
     * @param weights Optional statistical weights for each measurement
     * <!-- ************************************************************************************** -->
     */
    void add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err,
                      const std::optional<PyArray>& weights = std::nullopt);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Logarithmically subsample data to reduce computational cost while preserving information.
     * @details Selects a subset of data points distributed logarithmically to maintain
     *          representation across different temporal or spectral scales, useful for
     *          reducing computational overhead in MCMC sampling while preserving the
     *          essential characteristics of the observational dataset.
     * @param data Input data array to be subsampled
     * @param num_order Number of logarithmic orders to span in the subsampling
     * @return std::vector<size_t> Indices of selected data points
     * <!-- ************************************************************************************** -->
     */
    static std::vector<size_t> logscale_screen(PyArray const& data, size_t num_order);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Populate consolidated data arrays from individual FluxData objects.
     * @details Converts the collection of FluxData structures into unified arrays suitable
     *          for vectorized operations and MCMC sampling, ensuring efficient access to
     *          observational data during parameter estimation.
     * <!-- ************************************************************************************** -->
     */
    void fill_data_arrays();

    /**
     * <!-- ************************************************************************************** -->
     * @brief Get the total number of individual data points across all datasets.
     * @details Counts all observational measurements in the complete multi-wavelength
     *          dataset, providing the total degrees of freedom for statistical analysis.
     * @return size_t Total number of observational data points
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] size_t data_points_num() const;

    Array times;            ///< Consolidated array of all observation times [seconds]
    Array frequencies;      ///< Consolidated array of all observation frequencies [Hz]
    Array fluxes;           ///< Consolidated array of all observed flux densities [mJy]
    Array errors;           ///< Consolidated array of all observational uncertainties [mJy]
    Array weights;          ///< Consolidated array of all statistical weights
    Array model_fluxes;     ///< Consolidated array of all model-predicted flux densities [mJy]
    double t_min{con::inf}; ///< Earliest observation time [seconds]
    double t_max{0};        ///< Latest observation time [seconds]

    std::vector<FluxData> flux_data; ///< Collection of individual observational datasets

  private:
    std::vector<std::tuple<double, double, double, double, double>> tuple_data; ///< Internal data storage tuples
};

/**
 * <!-- ************************************************************************************** -->
 * @struct Params
 * @brief Complete parameter set for gamma-ray burst afterglow model configuration.
 * @details This structure encapsulates all physical and observational parameters needed
 *          to specify a gamma-ray burst afterglow model, including jet geometry, energy
 *          distribution, external medium properties, and radiation physics. These parameters
 *          serve as the variables for MCMC sampling and likelihood evaluation in parameter
 *          estimation studies.
 * <!-- ************************************************************************************** -->
 */
struct Params {
    double theta_v{0}; ///< Viewing angle between jet axis and line of sight [radians]

    // External medium parameters
    double n_ism{0};     ///< Constant ISM number density [cm^-3]
    double n0{con::inf}; ///< Inner boundary density for stratified wind [cm^-3]
    double A_star{0};    ///< Wind parameter A* in units of 5×10^11 g/cm [dimensionless]
    double k_m{2};       ///< Wind density power-law index (ρ ∝ r^-k_m)

    // Jet core parameters
    double E_iso{1e52};             ///< Isotropic-equivalent energy [erg]
    double Gamma0{300};             ///< Initial bulk Lorentz factor
    double theta_c{0.1};            ///< Core opening angle [radians]
    double k_e{2};                  ///< Energy power-law index for structured jets
    double k_g{2};                  ///< Lorentz factor power-law index for structured jets
    double duration{1 * unit::sec}; ///< Central engine activity duration [seconds]

    // Jet wing/wide component parameters
    double E_iso_w{1e52};        ///< Wing component isotropic-equivalent energy [erg]
    double Gamma0_w{300};        ///< Wing component initial bulk Lorentz factor
    double theta_w{con::pi / 2}; ///< Maximum polar angle for calculations [radians]

    // Magnetar energy injection parameters
    double L0{0}; ///< Characteristic magnetar luminosity [erg/s]
    double t0{1}; ///< Magnetar spin-down time scale [seconds]
    double q{2};  ///< Magnetar power-law index for energy injection

    // Forward shock radiation parameters
    double p{2.3};      ///< Electron energy spectral index
    double eps_e{0.1};  ///< Fraction of shock energy in relativistic electrons
    double eps_B{0.01}; ///< Fraction of shock energy in magnetic field
    double xi_e{1};     ///< Fraction of electrons accelerated to relativistic energies

    // Reverse shock radiation parameters
    double p_r{2.3};      ///< Reverse shock electron energy spectral index
    double eps_e_r{0.1};  ///< Reverse shock fraction of energy in electrons
    double eps_B_r{0.01}; ///< Reverse shock fraction of energy in magnetic field
    double xi_e_r{1};     ///< Reverse shock fraction of accelerated electrons
};

/**
 * <!-- ************************************************************************************** -->
 * @brief Compile-time array mapping parameter indices to Params member pointers.
 * @details Enables O(1) parameter assignment by index without string comparisons.
 *          Index order: theta_v(0), n_ism(1), n0(2), A_star(3), k_m(4), E_iso(5), Gamma0(6),
 *          theta_c(7), k_e(8), k_g(9), duration(10), E_iso_w(11), Gamma0_w(12), theta_w(13),
 *          L0(14), t0(15), q(16), p(17), eps_e(18), eps_B(19), xi_e(20), p_r(21),
 *          eps_e_r(22), eps_B_r(23), xi_e_r(24)
 * <!-- ************************************************************************************** -->
 */
inline constexpr std::array<double Params::*, 25> kParamPtrs = {
    &Params::theta_v,  &Params::n_ism,   &Params::n0,       &Params::A_star,  &Params::k_m,
    &Params::E_iso,    &Params::Gamma0,  &Params::theta_c,  &Params::k_e,     &Params::k_g,
    &Params::duration, &Params::E_iso_w, &Params::Gamma0_w, &Params::theta_w, &Params::L0,
    &Params::t0,       &Params::q,       &Params::p,        &Params::eps_e,   &Params::eps_B,
    &Params::xi_e,     &Params::p_r,     &Params::eps_e_r,  &Params::eps_B_r, &Params::xi_e_r,
};

/// Map parameter name string to index for Python interface setup
int param_name_to_index(const std::string& name);

/**
 * <!-- ************************************************************************************** -->
 * @struct ParamDef
 * @brief Parameter definition for MCMC sampling.
 * @details Specifies parameter name, bounds, scale (linear/log), and optional initial value.
 * <!-- ************************************************************************************** -->
 */
struct ParamDef {
    std::string name;
    double lower;
    double upper;
    int scale; // 0=linear, 1=log, 2=fixed
    double initial;
};

/**
 * <!-- ************************************************************************************** -->
 * @class ParamTransformer
 * @brief Transforms MCMC sampler arrays to Params structs with parameter scaling.
 * @details Handles conversion from sampler coordinates (with log-transformed parameters)
 *          to physical Params structs. Manages both sampled and fixed parameters efficiently.
 *          Designed for both single-sample transformations (flux predictions) and batch
 *          processing (MCMC likelihood evaluations).
 * <!-- ************************************************************************************** -->
 */
class ParamTransformer {
  public:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct transformer from parameter definitions.
     * @param param_defs Vector of parameter definitions with names, bounds, scales
     * <!-- ************************************************************************************** -->
     */
    explicit ParamTransformer(std::vector<ParamDef> const& param_defs);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct transformer from internal state (for unpickling).
     * @param sampled_indices Indices of sampled parameters
     * @param log_mask Log transformation mask for sampled parameters
     * @param fixed_indices Indices of fixed parameters
     * @param fixed_values Values for fixed parameters
     * <!-- ************************************************************************************** -->
     */
    ParamTransformer(std::vector<int> sampled_indices, std::vector<bool> log_mask, std::vector<int> fixed_indices,
                     std::vector<double> fixed_values);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Transform single sampler array to Params struct.
     * @param x Array of sampler coordinates (ndim elements, in log space where specified)
     * @return Params Complete parameter struct with all values set
     * <!-- ************************************************************************************** -->
     */
    Params operator()(PyArray const& x) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Get number of free (sampled) parameters.
     * @return size_t Number of dimensions in sampler space
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] size_t ndim() const { return sampled_indices_.size(); }

    /**
     * <!-- ************************************************************************************** -->
     * @brief Transform single sample to Params (internal helper for batch processing).
     * @param sample Pointer to sampler coordinates
     * @param out Output Params struct to fill
     * <!-- ************************************************************************************** -->
     */
    void transform_sample(double const* sample, Params& out) const;

  public:
    // Public state for pickle support and potential direct access
    std::vector<int> sampled_indices_; ///< Indices of sampled parameters in kParamPtrs
    std::vector<bool> log_mask_;       ///< Whether each sampled param needs 10^x transform
    std::vector<int> fixed_indices_;   ///< Indices of fixed parameters in kParamPtrs
    std::vector<double> fixed_values_; ///< Values for fixed parameters
};

/**
 * <!-- ************************************************************************************** -->
 * @struct ConfigParams
 * @brief Configuration parameters for numerical setup and advanced physics in afterglow modeling.
 * @details This structure controls the numerical resolution, cosmological setup, and optional
 *          advanced physics processes in gamma-ray burst afterglow calculations. These parameters
 *          determine the computational accuracy, efficiency, and the physical processes included
 *          in the simulation without affecting the fundamental model parameters.
 * <!-- ************************************************************************************** -->
 */
struct ConfigParams {
    double lumi_dist{1e26};    ///< Luminosity distance to the gamma-ray burst [cm]
    double z{0};               ///< Cosmological redshift
    std::string medium{"ism"}; ///< External medium type: "ism" or "wind"
    std::string jet{"tophat"}; ///< Jet structure type: "tophat", "gaussian", "powerlaw", etc.
    Real phi_resol{0.3};       ///< Azimuthal resolution: points per degree
    Real theta_resol{1};       ///< Polar resolution: points per degree
    Real t_resol{10};          ///< Temporal resolution: points per decade
    double rtol{1e-6};         ///< Relative tolerance for numerical integration
    bool rvs_shock{false};     ///< Whether to include reverse shock emission
    bool fwd_ssc{false};       ///< Whether to include forward shock synchrotron self-Compton
    bool rvs_ssc{false};       ///< Whether to include reverse shock synchrotron self-Compton
    bool ssc_cooling{false};   ///< Whether to include inverse Compton cooling of electrons
    bool kn{false};            ///< Whether to include Klein-Nishina corrections
    bool magnetar{false};      ///< Whether to include magnetar energy injection
};

/**
 * <!-- ************************************************************************************** -->
 * @struct BatchWorkspace
 * @brief Thread-local workspace for batch chi-squared computation.
 * @details Holds pre-allocated mutable arrays for model outputs, avoiding repeated allocations
 *          in the hot loop. Each thread gets its own workspace while sharing read-only
 *          observation data via const references.
 * <!-- ************************************************************************************** -->
 */
struct BatchWorkspace {
    Array model_fluxes;                  ///< Pre-allocated array for model flux predictions
    std::vector<Array> flux_data_models; ///< Pre-allocated arrays for flux_data Fv_model
    bool initialized{false};             ///< Whether workspace has been sized

    /// Initialize workspace arrays to match observation data dimensions
    void initialize(size_t tuple_size, std::vector<size_t> const& flux_data_sizes) {
        if (tuple_size > 0) {
            model_fluxes = Array::from_shape({tuple_size});
        }
        flux_data_models.resize(flux_data_sizes.size());
        for (size_t i = 0; i < flux_data_sizes.size(); ++i) {
            flux_data_models[i] = Array::from_shape({flux_data_sizes[i]});
        }
        initialized = true;
    }
};

/**
 * <!-- ************************************************************************************** -->
 * @struct MultiBandModel
 * @brief Primary interface for MCMC-based parameter estimation of gamma-ray burst afterglows.
 * @details This class provides a high-level interface for fitting gamma-ray burst afterglow
 *          models to multi-wavelength observational data using Markov Chain Monte Carlo methods.
 *          It combines theoretical model calculations with observational datasets to enable
 *          parameter estimation, uncertainty quantification, and model comparison for afterglow
 *          physics studies.
 * <!-- ************************************************************************************** -->
 */
struct MultiBandModel {
    MultiBandModel() = delete; ///< Default construction disabled - requires observational data

    /**
     * <!-- ************************************************************************************** -->
     * @brief Construct MCMC model interface with observational dataset.
     * @details Initializes the parameter estimation framework with a complete multi-wavelength
     *          observational dataset, preparing the model for likelihood evaluation and MCMC
     *          sampling across the parameter space.
     * @param data Complete multi-wavelength observational dataset for fitting
     * <!-- ************************************************************************************** -->
     */
    explicit MultiBandModel(MultiBandData data);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Lightweight constructor for batch processing mode.
     * @details Creates a model configured for use with estimate_chi2_with_workspace().
     *          Does not copy observation data - only stores configuration.
     * @param cfg Configuration parameters for physics calculations
     * <!-- ************************************************************************************** -->
     */
    explicit MultiBandModel(ConfigParams const& cfg);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Configure numerical and physics settings for the afterglow model.
     * @details Sets up the computational parameters including grid resolution, tolerance
     *          settings, and optional advanced physics processes that will be used during
     *          MCMC sampling and likelihood evaluations.
     * @param param Configuration parameters for numerical setup and physics options
     * <!-- ************************************************************************************** -->
     */
    void configure(ConfigParams const& param);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Batch chi-squared computation with clean transformer-based interface.
     * @details Computes chi-squared values for multiple parameter sets simultaneously.
     *          - Single sample (nwalkers=1): Uses optimized serial execution
     *          - Small batches (nwalkers<4): Uses serial execution to avoid OpenMP overhead
     *          - Large batches (nwalkers>=4): Uses OpenMP parallelization for speedup
     *
     *          Parameter transformation is handled by the ParamTransformer object.
     * @param samples 2D array of shape (nwalkers, ndim) containing sampler coordinates
     * @param transformer ParamTransformer to convert sampler coordinates to Params
     * @return PyArray 1D array of chi-squared values, one per walker
     * <!-- ************************************************************************************** -->
     */
    PyArray batch_estimate_chi2(PyGrid const& samples, ParamTransformer const& transformer);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate theoretical flux predictions for given parameters and observation grid.
     * @details Computes the complete theoretical afterglow flux predictions at specified times
     *          and frequencies using the given parameter values, providing detailed model
     *          predictions for comparison with observations or posterior predictive analysis.
     * @param param Complete set of physical model parameters
     * @param t Array of observation times [seconds]
     * @param nu Array of observation frequencies [Hz]
     * @return PyGrid Flux prediction grid with time-frequency dimensions
     * <!-- ************************************************************************************** -->
     */
    PyGrid flux_density_grid(Params const& param, PyArray const& t, PyArray const& nu);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Calculate frequency-integrated flux over specified bands and times.
     * @details Computes broadband flux predictions integrated over specified frequency ranges,
     *          useful for modeling bolometric observations or filter-integrated measurements
     *          where the effective frequency depends on the spectral shape.
     * @param param Complete set of physical model parameters
     * @param t Array of observation times [seconds]
     * @param nu_min Lower frequency bound for integration [Hz]
     * @param nu_max Upper frequency bound for integration [Hz]
     * @param num_points Number of frequency sampling points for integration
     * @return PyArray Time series of integrated flux values
     * <!-- ************************************************************************************** -->
     */
    PyArray flux(Params const& param, PyArray const& t, double nu_min, double nu_max, size_t num_points);

  private:
    /**
     * <!-- ************************************************************************************** -->
     * @brief Generate synchrotron and inverse Compton photon grids for emission calculations.
     * @details Template method that computes the fundamental radiation components (synchrotron
     *          and optionally inverse Compton) from both forward and reverse shocks, providing
     *          the basis for all flux calculations and likelihood evaluations.
     * @param param Complete set of physical model parameters
     * @param t_min Minimum time for photon grid generation [seconds]
     * @param t_max Maximum time for photon grid generation [seconds]
     * @param obs Observer object for flux calculations
     * @param fwd_photons Forward shock synchrotron photon grid
     * @param rvs_photons Reverse shock synchrotron photon grid
     * @param fwd_IC_photons Forward shock inverse Compton photon grid
     * @param rvs_IC_photons Reverse shock inverse Compton photon grid
     * <!-- ************************************************************************************** -->
     */
    template <typename ICPhotonGrid>
    void generate_photons(Params const& param, double t_min, double t_max, Observer& obs, SynPhotonGrid& fwd_photons,
                          SynPhotonGrid& rvs_photons, ICPhotonGrid& fwd_IC_photons, ICPhotonGrid& rvs_IC_photons);

    /**
     * <!-- ************************************************************************************** -->
     * @brief Select appropriate jet model based on configuration parameters.
     * @details Creates the appropriate jet structure (top-hat, Gaussian, power-law, etc.)
     *          using the specified parameters, enabling flexible jet geometry modeling
     *          within the MCMC framework.
     * @param param Complete set of physical model parameters
     * @return Ejecta Configured jet structure for afterglow calculations
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Ejecta select_jet(Params const& param) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Select appropriate external medium model based on configuration parameters.
     * @details Creates the appropriate circumburst medium (ISM or wind) using the specified
     *          parameters, enabling flexible environmental modeling within the MCMC framework.
     * @param param Complete set of physical model parameters
     * @return Medium Configured external medium for afterglow calculations
     * <!-- ************************************************************************************** -->
     */
    [[nodiscard]] Medium select_medium(Params const& param) const;

    /**
     * <!-- ************************************************************************************** -->
     * @brief Internal chi-squared estimation using pre-allocated workspace.
     * @details Helper method for batch processing that writes model outputs to a pre-allocated
     *          workspace, avoiding repeated memory allocations. Used internally by both
     *          estimate_chi2() and batch_estimate_chi2().
     * @param param Complete set of physical model parameters
     * @param obs_data_ref Const reference to observation data (shared across threads)
     * @param workspace Pre-allocated workspace for model outputs (thread-local)
     * @return double Chi-squared value for goodness-of-fit assessment
     * <!-- ************************************************************************************** -->
     */
    double estimate_chi2_with_workspace(Params const& param, MultiBandData const& obs_data_ref,
                                        BatchWorkspace& workspace);

    MultiBandData obs_data; ///< Multi-wavelength observational dataset for fitting
    ConfigParams config;    ///< Numerical and physics configuration parameters
};

template <typename ICPhotonGrid>
void MultiBandModel::generate_photons(Params const& param, double t_min, double t_max, Observer& obs,
                                      SynPhotonGrid& fwd_photons, SynPhotonGrid& rvs_photons,
                                      ICPhotonGrid& fwd_IC_photons, ICPhotonGrid& rvs_IC_photons) {
    Real theta_v = param.theta_v;
    Real theta_w = param.theta_w;
    RadParams rad;
    rad.p = param.p;
    rad.eps_e = param.eps_e;
    rad.eps_B = param.eps_B;
    rad.xi_e = param.xi_e;

    Real lumi_dist = config.lumi_dist * unit::cm;
    Real z = config.z;

    Medium medium = select_medium(param);
    Ejecta jet = select_jet(param);

    Real t_resol = config.t_resol;
    Real theta_resol = config.theta_resol;
    Real phi_resol = config.phi_resol;

    Array t_eval = xt::linspace<Real>(t_min, t_max, 5);

    if (config.rvs_shock == false) {
        auto coord = auto_grid(jet, t_eval, theta_w, theta_v, z, phi_resol, theta_resol, t_resol, true, 0, 48, 0.7);

        auto shock = generate_fwd_shock(coord, medium, jet, rad, config.rtol);

        obs.observe(coord, shock, lumi_dist, z);

        auto electrons = generate_syn_electrons(shock);
        fwd_photons = generate_syn_photons(shock, electrons);

        if (config.ssc_cooling) {
            if (config.kn) {
                KN_cooling(electrons, fwd_photons, shock);
            } else {
                Thomson_cooling(electrons, fwd_photons, shock);
            }
        }

        if (config.fwd_ssc) {
            fwd_IC_photons = generate_IC_photons(electrons, fwd_photons, config.kn);
        }

    } else {
        auto coord = auto_grid(jet, t_eval, theta_w, theta_v, z, phi_resol, theta_resol, t_resol, true, 0, 56, 0.3);

        RadParams rad_rvs;

        rad_rvs.p = param.p_r;
        rad_rvs.eps_e = param.eps_e_r;
        rad_rvs.eps_B = param.eps_B_r;
        rad_rvs.xi_e = param.xi_e_r;

        auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad, rad_rvs, config.rtol);

        obs.observe(coord, f_shock, lumi_dist, z);

        auto f_electrons = generate_syn_electrons(f_shock);
        fwd_photons = generate_syn_photons(f_shock, f_electrons);

        auto r_electrons = generate_syn_electrons(r_shock);
        rvs_photons = generate_syn_photons(r_shock, r_electrons);

        if (config.ssc_cooling) {
            if (config.kn) {
                KN_cooling(f_electrons, fwd_photons, f_shock);
                KN_cooling(r_electrons, rvs_photons, r_shock);
            } else {
                Thomson_cooling(f_electrons, fwd_photons, f_shock);
                Thomson_cooling(r_electrons, rvs_photons, r_shock);
            }
        }

        if (config.fwd_ssc) {
            fwd_IC_photons = generate_IC_photons(f_electrons, fwd_photons, config.kn);
        }

        if (config.rvs_ssc) {
            rvs_IC_photons = generate_IC_photons(r_electrons, rvs_photons, config.kn);
        }
    }
}
