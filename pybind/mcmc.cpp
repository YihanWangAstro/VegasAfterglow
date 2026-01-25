//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mcmc.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef HAVE_OPENMP
    #include <omp.h>
#endif

#include "error_handling.h"
#include "pybind.h"

// Parameter name to index mapping (computed once, used for Python interface setup)
int param_name_to_index(const std::string& name) {
    static const std::unordered_map<std::string, int> name_map = {
        {"theta_v", 0},  {"n_ism", 1},     {"n0", 2},       {"A_star", 3}, {"k_m", 4},       {"E_iso", 5},
        {"Gamma0", 6},   {"theta_c", 7},   {"k_e", 8},      {"k_g", 9},    {"duration", 10}, {"tau", 10},
        {"E_iso_w", 11}, {"Gamma0_w", 12}, {"theta_w", 13}, {"L0", 14},    {"t0", 15},       {"q", 16},
        {"p", 17},       {"eps_e", 18},    {"eps_B", 19},   {"xi_e", 20},  {"p_r", 21},      {"eps_e_r", 22},
        {"eps_B_r", 23}, {"xi_e_r", 24},
    };
    auto it = name_map.find(name);
    if (it == name_map.end()) {
        return -1; // Invalid parameter name
    }
    return it->second;
}

std::vector<size_t> MultiBandData::logscale_screen(PyArray const& data, size_t num_order) {
    const size_t total_size = data.size();

    // Handle edge cases: empty or single-element arrays
    if (total_size == 0) {
        return {};
    }
    if (total_size == 1) {
        return {0};
    }

    if (num_order == 0) {
        std::vector<size_t> indices(total_size);
        std::iota(indices.begin(), indices.end(), 0);
        return indices;
    }

    const double log_start = std::log10(static_cast<double>(data(0)));
    const double log_end = std::log10(static_cast<double>(data(total_size - 1)));
    const double log_range = log_end - log_start;
    const size_t total_points = static_cast<size_t>(std::ceil(log_range * static_cast<double>(num_order))) + 1;

    std::vector<size_t> indices;
    indices.reserve(total_points);

    // Always include the first point
    indices.push_back(0);

    // Handle case where total_points <= 1 (avoid division by zero)
    if (total_points <= 1) {
        if (total_size > 1) {
            indices.push_back(total_size - 1);
        }
        return indices;
    }

    const double step = log_range / static_cast<double>(total_points - 1);

    for (size_t i = 1; i < total_points - 1; ++i) {
        const double log_target = log_start + static_cast<double>(i) * step;
        const double target_value = std::pow(10.0, log_target);

        size_t best_idx = 1;
        double min_diff = std::abs(static_cast<double>(data(1)) - target_value);

        for (size_t j = 2; j < total_size - 1; ++j) {
            const double diff = std::abs(static_cast<double>(data(j)) - target_value);
            if (diff < min_diff) {
                min_diff = diff;
                best_idx = j;
            }
        }

        if (std::ranges::find(indices, best_idx) == indices.end()) {
            indices.push_back(best_idx);
        }
    }

    // Always include the last point
    if (total_size > 1) {
        indices.push_back(total_size - 1);
    }

    std::ranges::sort(indices);
    indices.erase(std::ranges::unique(indices).begin(), indices.end());

    return indices;
}

double FluxData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        const double error = Fv_err(i);
        //if (error == 0)
        //    continue;
        const double diff = Fv_obs(i) - Fv_model(i);
        chi_square += weights(i) * (diff * diff) / (error * error);
    }
    return chi_square;
}

double MultiBandData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < times.size(); ++i) {
        const double error = errors(i);
        //if (error == 0)
        //    continue;
        const double diff = fluxes(i) - model_fluxes(i);
        chi_square += weights(i) * (diff * diff) / (error * error);
    }
    for (auto& d : flux_data) {
        chi_square += d.estimate_chi2();
    }

    return chi_square;
}

Ejecta MultiBandModel::select_jet(Params const& param) const {
    const Real eps_iso = param.E_iso * unit::erg / (4 * con::pi);
    const Real Gamma0 = param.Gamma0;
    const Real theta_c = param.theta_c;
    const Real theta_w = param.theta_w;
    const Real eps_iso_w = param.E_iso_w * unit::erg / (4 * con::pi);
    const Real Gamma0_w = param.Gamma0_w;
    Ejecta jet;
    jet.T0 = param.duration * unit::sec;
    if (config.jet == "tophat") {
        jet.eps_k = math::tophat(theta_c, eps_iso);
        jet.Gamma0 = math::tophat_plus_one(theta_c, Gamma0 - 1);
    } else if (config.jet == "gaussian") {
        jet.eps_k = math::gaussian(theta_c, eps_iso);
        jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
    } else if (config.jet == "powerlaw") {
        jet.eps_k = math::powerlaw(theta_c, eps_iso, param.k_e);
        jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, param.k_g);
    } else if (config.jet == "powerlaw_wing") {
        jet.eps_k = math::powerlaw_wing(theta_c, eps_iso_w, param.k_e);
        jet.Gamma0 = math::powerlaw_wing_plus_one(theta_c, Gamma0_w - 1, param.k_g);
    } else if (config.jet == "uniform") {
        jet.eps_k = math::tophat(con::pi / 2, eps_iso);
        jet.Gamma0 = math::tophat_plus_one(con::pi / 2, Gamma0 - 1);
    } else if (config.jet == "two_component") {
        jet.eps_k = math::two_component(theta_c, theta_w, eps_iso, eps_iso_w);
        jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0 - 1, Gamma0_w - 1);
    } else if (config.jet == "step_powerlaw") {
        jet.eps_k = math::step_powerlaw(theta_c, eps_iso, eps_iso_w, param.k_e);
        jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma0 - 1, Gamma0_w - 1, param.k_g);
    } else {
        AFTERGLOW_ENSURE(false, "Unknown jet type");
    }

    if (config.magnetar == true) {
        jet.deps_dt =
            math::magnetar_injection(param.t0 * unit::sec, param.q, param.L0 * unit::erg / unit::sec, theta_c);
    }
    return jet;
}

Medium MultiBandModel::select_medium(Params const& param) const {
    Medium medium;
    if (config.medium == "ism") {
        medium.rho = evn::ISM(param.n_ism / unit::cm3);
    } else if (config.medium == "wind") {
        medium.rho = evn::wind(param.A_star, param.n_ism / unit::cm3, param.n0 / unit::cm3, param.k_m);
    } else {
        AFTERGLOW_ENSURE(false, "Unknown medium type");
    }
    return medium;
}

void MultiBandData::add_flux_density(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err,
                                     std::optional<PyArray> const& weights) {
    AFTERGLOW_REQUIRE(t.size() == Fv_obs.size() && t.size() == Fv_err.size(), "light curve array inconsistent length!");

    Array w = xt::ones<Real>({t.size()});

    if (weights) {
        w = *weights;
        AFTERGLOW_REQUIRE(t.size() == w.size(), "weights array inconsistent length!");
    }

    for (size_t i = 0; i < t.size(); ++i) {
        tuple_data.emplace_back(t(i) * unit::sec, nu * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                Fv_err(i) * unit::flux_den_cgs, w(i));
    }
}

void MultiBandData::add_flux(double nu_min, double nu_max, size_t num_points, PyArray const& t, PyArray const& Fv_obs,
                             PyArray const& Fv_err, const std::optional<PyArray>& weights) {
    AFTERGLOW_REQUIRE(t.size() == Fv_obs.size() && t.size() == Fv_err.size(), "light curve array inconsistent length!");
    AFTERGLOW_REQUIRE(is_ascending(t), "Time array must be in ascending order!");
    AFTERGLOW_REQUIRE(nu_min < nu_max, "nu_min must be less than nu_max!");

    Array w = xt::ones<Real>({t.size()});

    if (weights) {
        w = *weights;
        AFTERGLOW_REQUIRE(t.size() == w.size(), "weights array inconsistent length!");

        const size_t len = w.size();
        Real weight_sum = 0;
        for (size_t i = 0; i < len; ++i) {
            weight_sum += w(i);
        }
        if (weight_sum > 0) {
            w /= (weight_sum / static_cast<double>(len));
        }
    }

    const Array nu = xt::logspace(std::log10(nu_min * unit::Hz), std::log10(nu_max * unit::Hz), num_points);

    flux_data.emplace_back(
        FluxData{t * unit::sec, nu, Fv_obs * unit::flux_cgs, Fv_err * unit::flux_cgs, xt::zeros<Real>({t.size()}), w});
}

void MultiBandData::add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err,
                                 const std::optional<PyArray>& weights) {
    AFTERGLOW_REQUIRE(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size(), "spectrum array inconsistent length!");

    Array w = xt::ones<Real>({nu.size()});

    if (weights) {
        w = *weights;
        AFTERGLOW_REQUIRE(nu.size() == w.size(), "weights array inconsistent length!");
    }

    for (size_t i = 0; i < nu.size(); ++i) {
        tuple_data.emplace_back(t * unit::sec, nu(i) * unit::Hz, Fv_obs(i) * unit::flux_den_cgs,
                                Fv_err(i) * unit::flux_den_cgs, w(i));
    }
}

size_t MultiBandData::data_points_num() const {
    size_t num = tuple_data.size();
    for (auto& d : flux_data) {
        num += d.t.size();
    }
    return num;
}

void MultiBandData::fill_data_arrays() {
    // Skip if arrays are already filled (e.g., from pickle deserialization)
    if (times.size() > 0 || (tuple_data.empty() && !flux_data.empty())) {
        return;
    }

    const size_t len = tuple_data.size();
    std::ranges::sort(tuple_data, [](auto const& a, auto const& b) { return std::get<0>(a) < std::get<0>(b); });
    times = Array::from_shape({len});
    frequencies = Array::from_shape({len});
    fluxes = Array::from_shape({len});
    errors = Array::from_shape({len});
    model_fluxes = Array::from_shape({len});
    weights = Array::from_shape({len});

    Real weight_sum = 0;
    for (size_t i = 0; i < len; ++i) {
        times(i) = std::get<0>(tuple_data[i]);
        frequencies(i) = std::get<1>(tuple_data[i]);
        fluxes(i) = std::get<2>(tuple_data[i]);
        errors(i) = std::get<3>(tuple_data[i]);
        weights(i) = std::get<4>(tuple_data[i]);
        model_fluxes(i) = 0; // Placeholder for model fluxes
        weight_sum += weights(i);
    }
    if (weight_sum > 0) {
        weights /= (weight_sum / static_cast<double>(len));
    }

    if (len > 0) {
        this->t_min = times.front();
        this->t_max = times.back();
    }

    for (auto& d : flux_data) {
        if (d.t.size() == 0)
            continue;
        if (d.t.front() < t_min)
            t_min = d.t.front();
        if (d.t.back() > t_max)
            t_max = d.t.back();
    }
}

MultiBandModel::MultiBandModel(MultiBandData data) : obs_data(std::move(data)) {
    obs_data.fill_data_arrays();

    AFTERGLOW_REQUIRE((obs_data.times.size() > 0 || !obs_data.flux_data.empty()), "No observation time data provided!");
}

MultiBandModel::MultiBandModel(ConfigParams const& cfg) : config(cfg) {
    // Lightweight constructor for batch processing - obs_data remains empty
    // Use with estimate_chi2_with_workspace() which takes obs_data as parameter
}

void MultiBandModel::configure(ConfigParams const& param) {
    this->config = param;
}

double MultiBandModel::estimate_chi2(Params const& param) {
    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, obs_data.t_min, obs_data.t_max, obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    obs_data.model_fluxes = obs.specific_flux_series(obs_data.times, obs_data.frequencies, f_photons);
    for (auto& d : obs_data.flux_data) {
        d.Fv_model = obs.flux(d.t, d.nu, f_photons);
    }

    if (r_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, r_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, r_photons);
        }
    }

    if (f_IC_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, f_IC_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, f_IC_photons);
        }
    }

    if (r_IC_photons.size() > 0) {
        obs_data.model_fluxes += obs.specific_flux_series(obs_data.times, obs_data.frequencies, r_IC_photons);
        for (auto& d : obs_data.flux_data) {
            d.Fv_model += obs.flux(d.t, d.nu, r_IC_photons);
        }
    }

    return obs_data.estimate_chi2();
}

double MultiBandModel::estimate_chi2_with_workspace(Params const& param, MultiBandData const& obs_data_ref,
                                                    BatchWorkspace& workspace) {
    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, obs_data_ref.t_min, obs_data_ref.t_max, obs, f_photons, r_photons, f_IC_photons,
                     r_IC_photons);

    // Compute model fluxes into workspace arrays (not obs_data)
    if (obs_data_ref.times.size() > 0) {
        workspace.model_fluxes = obs.specific_flux_series(obs_data_ref.times, obs_data_ref.frequencies, f_photons);
    }
    for (size_t i = 0; i < obs_data_ref.flux_data.size(); ++i) {
        workspace.flux_data_models[i] = obs.flux(obs_data_ref.flux_data[i].t, obs_data_ref.flux_data[i].nu, f_photons);
    }

    if (r_photons.size() > 0) {
        if (obs_data_ref.times.size() > 0) {
            workspace.model_fluxes += obs.specific_flux_series(obs_data_ref.times, obs_data_ref.frequencies, r_photons);
        }
        for (size_t i = 0; i < obs_data_ref.flux_data.size(); ++i) {
            workspace.flux_data_models[i] +=
                obs.flux(obs_data_ref.flux_data[i].t, obs_data_ref.flux_data[i].nu, r_photons);
        }
    }

    if (f_IC_photons.size() > 0) {
        if (obs_data_ref.times.size() > 0) {
            workspace.model_fluxes +=
                obs.specific_flux_series(obs_data_ref.times, obs_data_ref.frequencies, f_IC_photons);
        }
        for (size_t i = 0; i < obs_data_ref.flux_data.size(); ++i) {
            workspace.flux_data_models[i] +=
                obs.flux(obs_data_ref.flux_data[i].t, obs_data_ref.flux_data[i].nu, f_IC_photons);
        }
    }

    if (r_IC_photons.size() > 0) {
        if (obs_data_ref.times.size() > 0) {
            workspace.model_fluxes +=
                obs.specific_flux_series(obs_data_ref.times, obs_data_ref.frequencies, r_IC_photons);
        }
        for (size_t i = 0; i < obs_data_ref.flux_data.size(); ++i) {
            workspace.flux_data_models[i] +=
                obs.flux(obs_data_ref.flux_data[i].t, obs_data_ref.flux_data[i].nu, r_IC_photons);
        }
    }

    // Compute chi-squared using workspace arrays and const observation data
    double chi_square = 0;

    // Chi2 from tuple data (times/frequencies/fluxes arrays)
    for (size_t i = 0; i < obs_data_ref.times.size(); ++i) {
        const double error = obs_data_ref.errors(i);
        const double diff = obs_data_ref.fluxes(i) - workspace.model_fluxes(i);
        chi_square += obs_data_ref.weights(i) * (diff * diff) / (error * error);
    }

    // Chi2 from flux_data entries
    for (size_t d = 0; d < obs_data_ref.flux_data.size(); ++d) {
        const auto& fd = obs_data_ref.flux_data[d];
        const auto& model = workspace.flux_data_models[d];
        for (size_t i = 0; i < fd.t.size(); ++i) {
            const double error = fd.Fv_err(i);
            const double diff = fd.Fv_obs(i) - model(i);
            chi_square += fd.weights(i) * (diff * diff) / (error * error);
        }
    }

    return chi_square;
}

auto MultiBandModel::flux_density_grid(Params const& param, PyArray const& t, PyArray const& nu) -> PyGrid {
    Array t_bins = t * unit::sec;
    Array nu_bins = nu * unit::Hz;
    MeshGrid F_nu = MeshGrid::from_shape({nu.size(), t.size()});

    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, t_bins.front(), t_bins.back(), obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    F_nu = obs.specific_flux(t_bins, nu_bins, f_photons);

    if (r_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, r_photons);
    }

    if (f_IC_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, f_IC_photons);
    }

    if (r_IC_photons.size() > 0) {
        F_nu += obs.specific_flux(t_bins, nu_bins, r_IC_photons);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_den_cgs;
}

auto MultiBandModel::flux(Params const& param, PyArray const& t, double nu_min, double nu_max, size_t num_points)
    -> PyArray {
    Array t_bins = t * unit::sec;
    Array nu_bins = xt::logspace(std::log10(nu_min * unit::Hz), std::log10(nu_max * unit::Hz), num_points);
    Array F_nu = Array::from_shape({t.size()});

    Observer obs;
    SynPhotonGrid f_photons;
    SynPhotonGrid r_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> f_IC_photons;
    ICPhotonGrid<SynElectrons, SynPhotons> r_IC_photons;

    generate_photons(param, t_bins.front(), t_bins.back(), obs, f_photons, r_photons, f_IC_photons, r_IC_photons);

    F_nu = obs.flux(t_bins, nu_bins, f_photons);

    if (r_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, r_photons);
    }

    if (f_IC_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, f_IC_photons);
    }

    if (r_IC_photons.size() > 0) {
        F_nu += obs.flux(t_bins, nu_bins, r_IC_photons);
    }

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_cgs;
}

PyArray MultiBandModel::batch_estimate_chi2(PyGrid const& samples, std::vector<int> const& param_indices,
                                            std::vector<bool> const& log_mask, std::vector<int> const& fixed_indices,
                                            std::vector<double> const& fixed_values) {
    const size_t nwalkers = samples.shape(0);
    const size_t ndim = samples.shape(1);

    AFTERGLOW_REQUIRE(ndim == param_indices.size(), "param_indices size must match ndim");
    AFTERGLOW_REQUIRE(ndim == log_mask.size(), "log_mask size must match ndim");
    AFTERGLOW_REQUIRE(fixed_indices.size() == fixed_values.size(),
                      "fixed_indices and fixed_values must have same size");

    // Output array
    Array chi2_results = Array::from_shape({nwalkers});

    // Shared read-only observation data (no copying needed)
    const MultiBandData& data_ref = obs_data;
    const ConfigParams& config_ref = config;

    // Pre-compute workspace dimensions (same for all threads)
    const size_t tuple_size = data_ref.times.size();
    std::vector<size_t> flux_data_sizes;
    flux_data_sizes.reserve(data_ref.flux_data.size());
    for (const auto& fd : data_ref.flux_data) {
        flux_data_sizes.push_back(fd.t.size());
    }

#ifdef HAVE_OPENMP
    #pragma omp parallel
    {
        // Thread-local workspace (only allocates model output arrays, not observation data)
        BatchWorkspace workspace;
        workspace.initialize(tuple_size, flux_data_sizes);

        // Lightweight model - only config, no obs_data copy
        MultiBandModel local_model(config_ref);

    #pragma omp for schedule(dynamic)
        for (size_t i = 0; i < nwalkers; ++i) {
#else
    {
        // Serial fallback - single workspace instance
        BatchWorkspace workspace;
        workspace.initialize(tuple_size, flux_data_sizes);

        // Lightweight model - only config, no obs_data copy
        MultiBandModel local_model(config_ref);

        for (size_t i = 0; i < nwalkers; ++i) {
#endif
            // Build Params struct using pointer-to-member for O(1) access
            Params p;

            // Set sampled parameters with optional log transformation
            for (size_t j = 0; j < ndim; ++j) {
                const double raw_val = samples(i, j);
                const double val = log_mask[j] ? std::pow(10.0, raw_val) : raw_val;
                p.*(kParamPtrs[static_cast<size_t>(param_indices[j])]) = val;
            }

            // Set fixed parameters
            for (size_t j = 0; j < fixed_indices.size(); ++j) {
                p.*(kParamPtrs[static_cast<size_t>(fixed_indices[j])]) = fixed_values[j];
            }

            // Compute chi2 using workspace (avoids copying obs_data)
            double chi2;
            try {
                chi2 = local_model.estimate_chi2_with_workspace(p, data_ref, workspace);
                // Check for NaN
                if (chi2 != chi2) {
                    chi2 = std::numeric_limits<double>::infinity();
                }
            } catch (...) {
                chi2 = std::numeric_limits<double>::infinity();
            }

            chi2_results(i) = chi2;
        }
    }

    // Acquire GIL before creating Python object
    pybind11::gil_scoped_acquire acquire;
    return chi2_results;
}
