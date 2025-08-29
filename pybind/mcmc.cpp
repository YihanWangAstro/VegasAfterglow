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

#include "pybind.h"

template <typename Array>
void sort_synchronized(Array& a, Array& b, Array& c) {
    std::size_t N = a.size();

    std::vector<std::size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);

    std::sort(idx.begin(), idx.end(), [&](std::size_t i, std::size_t j) { return a(i) < a(j); });

    Array a_sorted = xt::empty<double>({N});
    Array b_sorted = xt::empty<double>({N});
    Array c_sorted = xt::empty<double>({N});

    for (std::size_t i = 0; i < N; ++i) {
        a_sorted(i) = a(idx[i]);
        b_sorted(i) = b(idx[i]);
        c_sorted(i) = c(idx[i]);
    }

    a = std::move(a_sorted);
    b = std::move(b_sorted);
    c = std::move(c_sorted);
}

double LightCurveData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < t.size(); ++i) {
        if (Fv_err[i] == 0) continue;
        double diff = Fv_obs[i] - Fv_model[i];
        chi_square += (diff * diff) / (Fv_err[i] * Fv_err[i]);
    }
    return chi_square;
}

double SpectrumData::estimate_chi2() const {
    double chi_square = 0;
    for (size_t i = 0; i < nu.size(); ++i) {
        if (Fv_err[i] == 0) continue;
        double diff = Fv_obs[i] - Fv_model[i];
        chi_square += (diff * diff) / (Fv_err[i] * Fv_err[i]);
    }
    return chi_square;
}

double MultiBandData::estimate_chi2() const {
    double chi_square = 0;
    for (auto const& data : light_curve) {
        chi_square += data.estimate_chi2();
    }
    for (auto const& data : spectrum) {
        chi_square += data.estimate_chi2();
    }
    return chi_square;
}

void MultiBandData::add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err) {
    assert(t.size() == Fv_obs.size() && t.size() == Fv_err.size() && "light curve array inconsistent length!");
    LightCurveData data;

    data.nu = nu * unit::Hz;
    data.t = t;
    data.Fv_obs = Fv_obs;
    data.Fv_err = Fv_err;
    data.Fv_model = xt::zeros_like(data.Fv_obs);

    sort_synchronized(data.t, data.Fv_obs, data.Fv_err);

    for (auto& t : data.t) {
        t *= unit::sec;
    }
    for (auto& Fv_obs : data.Fv_obs) {
        Fv_obs *= unit::flux_den_cgs;
    }
    for (auto& Fv_err : data.Fv_err) {
        Fv_err *= unit::flux_den_cgs;
    }

    light_curve.push_back(std::move(data));
}

void MultiBandData::add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err) {
    assert(nu.size() == Fv_obs.size() && nu.size() == Fv_err.size() && "spectrum array inconsistent length!");
    SpectrumData data;

    data.t = t * unit::sec;
    data.nu = nu;
    data.Fv_obs = Fv_obs;
    data.Fv_err = Fv_err;

    sort_synchronized(data.nu, data.Fv_obs, data.Fv_err);

    for (auto& nu : data.nu) {
        nu *= unit::Hz;
    }
    for (auto& Fv_obs : data.Fv_obs) {
        Fv_obs *= unit::flux_den_cgs;
    }
    for (auto& Fv_err : data.Fv_err) {
        Fv_err *= unit::flux_den_cgs;
    }
    spectrum.push_back(std::move(data));
}

MultiBandModel::MultiBandModel(MultiBandData const& obs_data) : obs_data(obs_data) {
    std::vector<Real> evaluate_t;
    evaluate_t.reserve(100);

    for (auto& data : obs_data.light_curve) {
        for (auto t : data.t) {
            evaluate_t.push_back(t);
        }
    }

    for (auto& data : obs_data.spectrum) {
        evaluate_t.push_back(data.t);
    }
    std::sort(evaluate_t.begin(), evaluate_t.end());
    this->t_eval = xt::eval(xt::adapt(evaluate_t));

    if (t_eval.size() == 0) {
        std::cerr << "Error: No observation time data provided!" << std::endl;
    }
}

void MultiBandModel::configure(ConfigParams const& param) { this->config = param; }

/*
std::vector<double> MultiBandModel::chiSquareBatch(std::vector<Params> const& param_batch) {
    const size_t N = param_batch.size();
    std::vector<double> results(N);

#pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        results[i] = chiSquare(param_batch[i]);
    }

    return results;
}
*/

void MultiBandModel::build_system(Params const& param, Array const& t_eval, Observer& obs, SynElectronGrid& electrons,
                                  SynPhotonGrid& photons) {
    Real E_iso = param.E_iso * unit::erg;
    Real Gamma0 = param.Gamma0;
    Real theta_c = param.theta_c;
    Real theta_v = param.theta_v;
    Real theta_w = param.theta_w;
    RadParams rad;

    rad.p = param.p;
    rad.eps_e = param.eps_e;
    rad.eps_B = param.eps_B;
    rad.xi_e = param.xi_e;

    Real lumi_dist = config.lumi_dist * unit::cm;
    Real z = config.z;

    // create model
    Medium medium;
    if (config.medium == "ism") {
        Real n_ism = param.n_ism / unit::cm3;
        std::tie(medium.rho, medium.mass) = evn::ISM(n_ism);
    } else if (config.medium == "wind") {
        std::tie(medium.rho, medium.mass) = evn::wind(param.A_star);
    } else {
        std::cerr << "Error: Unknown medium type" << std::endl;
    }

    Ejecta jet;
    if (config.jet == "tophat") {
        jet.eps_k = math::tophat(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::tophat(theta_c, Gamma0);
    } else if (config.jet == "gaussian") {
        jet.eps_k = math::gaussian(theta_c, E_iso / (4 * con::pi));
        jet.Gamma0 = math::gaussian(theta_c, Gamma0);
    } else if (config.jet == "powerlaw") {
        jet.eps_k = math::powerlaw(theta_c, E_iso / (4 * con::pi), param.k_jet);
        jet.Gamma0 = math::powerlaw(theta_c, Gamma0, param.k_jet);
    } else if (config.jet == "powerlawgeneral") {
        jet.eps_k = math::powerlaw(param.theta_c, E_iso / (4 * con::pi), param.k_jet);
        jet.Gamma0 = math::powerlaw(param.theta_cG, Gamma0, param.kG_jet);
    } else {
        std::cerr << "Error: Unknown jet type" << std::endl;
    }

    Real t_resol = config.t_resol;
    Real theta_resol = config.theta_resol;
    Real phi_resol = config.phi_resol;

    auto coord = auto_grid(jet, t_eval, theta_w, theta_v, z, phi_resol, theta_resol, t_resol);

    auto shock = generate_fwd_shock(coord, medium, jet, rad, config.rtol);

    // obs.observe_at(t_eval, coord, shock, lumi_dist, z);
    obs.observe(coord, shock, lumi_dist, z);

    electrons = generate_syn_electrons(shock);

    photons = generate_syn_photons(shock, electrons);
}

double MultiBandModel::estimate_chi2(Params const& param) {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    build_system(param, t_eval, obs, electrons, photons);

    for (auto& data : obs_data.light_curve) {
        data.Fv_model = obs.specific_flux(data.t, data.nu, photons);
    }

    for (auto& data : obs_data.spectrum) {
        data.Fv_model = obs.spectrum(data.nu, data.t, photons);
    }

    return obs_data.estimate_chi2();
}

auto MultiBandModel::light_curves(Params const& param, PyArray const& t, PyArray const& nu) -> PyGrid {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    Array t_bins = t * unit::sec;

    build_system(param, t_bins, obs, electrons, photons);
    Array nu_bins = nu * unit::Hz;
    auto F_nu = obs.specific_flux(t_bins, nu_bins, photons);

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_den_cgs;
}

auto MultiBandModel::spectra(Params const& param, PyArray const& nu, PyArray const& t) -> PyGrid {
    Observer obs;
    SynElectronGrid electrons;
    SynPhotonGrid photons;

    Array t_bins = t * unit::sec;
    build_system(param, t_bins, obs, electrons, photons);
    Array nu_bins = nu * unit::Hz;
    auto F_nu = obs.spectra(nu_bins, t_bins, photons);

    // we bind this function for GIL free. As the return will create a pyobject, we need to get the GIL.
    pybind11::gil_scoped_acquire acquire;
    return F_nu / unit::flux_den_cgs;
}