//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "mcmc.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <utility>
#include <vector>

#include "error_handling.h"
#include "pybind.h"

std::vector<size_t> MultiBandData::logscale_screen(PyArray const& data, size_t num_order) {
    const size_t total_size = data.size();

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
        w /= (weight_sum / static_cast<double>(len));
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
    weights /= (weight_sum / static_cast<double>(len));

    if (len > 0) {
        this->t_min = times.front();
        this->t_max = times.back();
    }

    for (auto& d : flux_data) {
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
