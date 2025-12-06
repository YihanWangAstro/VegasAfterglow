//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include <algorithm>
#include <numeric>

#include "../include/afterglow.h"
#include "error_handling.h"
#include "xtensor/misc/xsort.hpp"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                   const std::optional<PyMagnetar>& magnetar) {
    Ejecta jet;
    jet.eps_k = math::tophat(theta_c, E_iso);
    jet.Gamma0 = math::tophat_plus_one(theta_c, Gamma0 - 1);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                     const std::optional<PyMagnetar>& magnetar) {
    Ejecta jet;
    jet.eps_k = math::gaussian(theta_c, E_iso);
    jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading, Real duration,
                     const std::optional<PyMagnetar>& magnetar) {
    Ejecta jet;
    jet.eps_k = math::powerlaw(theta_c, E_iso, k_e);
    jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, k_g);
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyPowerLawWing(Real theta_c, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g, bool spreading, Real duration) {
    Ejecta jet;
    jet.eps_k = math::powerlaw_wing(theta_c, E_iso_w, k_e);
    jet.Gamma0 = math::powerlaw_wing_plus_one(theta_c, Gamma0_w - 1, k_g);
    jet.spreading = spreading;
    jet.T0 = duration;

    return jet;
}

Ejecta PyStepPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g,
                         bool spreading, Real duration, const std::optional<PyMagnetar>& magnetar) {
    Ejecta jet;
    jet.eps_k = math::step_powerlaw(theta_c, E_iso, E_iso_w, k_e);
    jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma0 - 1, Gamma0_w - 1, k_g);

    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Ejecta PyTwoComponentJet(Real theta_c, Real E_iso, Real Gamma0, Real theta_w, Real E_iso_w, Real Gamma0_w,
                         bool spreading, Real duration, const std::optional<PyMagnetar>& magnetar) {
    Ejecta jet;
    jet.eps_k = math::two_component(theta_c, theta_w, E_iso, E_iso_w);

    jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0 - 1, Gamma0_w - 1);

    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }

    return jet;
}

Medium PyISM(Real n_ism) {
    Medium medium;
    medium.rho = [=](Real phi, Real theta, Real r) { return n_ism * 1.67e-24; };

    return medium;
}

Medium PyWind(Real A_star, Real n_ism, Real n0, Real k) {
    Medium medium;

    constexpr Real r0 = 1e17; // reference radius
    const Real A = A_star * 5e11 * std::pow(r0, k - 2);

    if (k >= 0) {
        const Real rho_ism = n_ism * 1.67e-24;
        const Real r0k = A / (n0 * 1.67e-24);
        medium.rho = [=](Real phi, Real theta, Real r) { return A / (r0k + std::pow(r, k)) + rho_ism; };
    } else {
        const Real rho0 = n0 * 1.67e-24;
        const Real r_ism = A / (n_ism * 1.67e-24);
        medium.rho = [=](Real phi, Real theta, Real r) { return A / (r_ism + std::pow(r, k)) + rho0; };
    }
    return medium;
}

void convert_unit(Ejecta& jet, Medium& medium) {
    const auto eps_k_cgs = jet.eps_k;
    jet.eps_k = [=](Real phi, Real theta) { return eps_k_cgs(phi, theta) * (unit::erg / (4 * con::pi)); };

    const auto deps_dt_cgs = jet.deps_dt;
    jet.deps_dt = [=](Real phi, Real theta, Real t) {
        return deps_dt_cgs(phi, theta, t / unit::sec) * (unit::erg / (4 * con::pi * unit::sec));
    };

    const auto dm_dt_cgs = jet.dm_dt;
    jet.dm_dt = [=](Real phi, Real theta, Real t) {
        return dm_dt_cgs(phi, theta, t / unit::sec) * (unit::g / (4 * con::pi * unit::sec));
    };

    jet.T0 *= unit::sec;

    const auto rho_cgs = medium.rho; // number density from python side
    medium.rho = [=](Real phi, Real theta, Real r) {
        return rho_cgs(phi, theta, r / unit::cm) * (unit::g / unit::cm3); // convert to density
    };
}

void save_shock_details(Shock const& shock, PyShock& details) {
    details.Gamma = shock.Gamma;
    details.Gamma_th = shock.Gamma_th;
    details.r = shock.r / unit::cm;
    details.t_comv = shock.t_comv / unit::sec;
    details.B_comv = shock.B / unit::Gauss;
    details.N_p = shock.N_p;
    details.theta = shock.theta;
}

template <typename ElectronGrid>
void save_electron_details(ElectronGrid const& electrons, PyShock& details) {
    auto shape = electrons.shape();

    details.gamma_m = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.gamma_c = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.gamma_a = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.gamma_M = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.gamma_m_hat = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.gamma_c_hat = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.N_e = xt::zeros<Real>({shape[0], shape[1], shape[2]});

    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                details.gamma_a(i, j, k) = electrons(i, j, k).gamma_a;
                details.gamma_m(i, j, k) = electrons(i, j, k).gamma_m;
                details.gamma_c(i, j, k) = electrons(i, j, k).gamma_c;
                details.gamma_M(i, j, k) = electrons(i, j, k).gamma_M;
                details.gamma_m_hat(i, j, k) = electrons(i, j, k).Ys.gamma_m_hat;
                details.gamma_c_hat(i, j, k) = electrons(i, j, k).Ys.gamma_c_hat;
                details.N_e(i, j, k) = electrons(i, j, k).N_e;
            }
        }
    }
}
template <typename PhotonGrid>
void save_photon_details(PhotonGrid const& photons, PyShock& details) {
    auto shape = photons.shape();

    details.nu_m = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.nu_c = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.nu_a = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.nu_M = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.nu_m_hat = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.nu_c_hat = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.I_nu_max = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    details.Y_T = xt::zeros<Real>({shape[0], shape[1], shape[2]});

    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                details.nu_a(i, j, k) = photons(i, j, k).nu_a / unit::Hz;
                details.nu_m(i, j, k) = photons(i, j, k).nu_m / unit::Hz;
                details.nu_c(i, j, k) = photons(i, j, k).nu_c / unit::Hz;
                details.nu_M(i, j, k) = photons(i, j, k).nu_M / unit::Hz;
                details.nu_m_hat(i, j, k) = photons(i, j, k).Ys.nu_m_hat / unit::Hz;
                details.nu_c_hat(i, j, k) = photons(i, j, k).Ys.nu_c_hat / unit::Hz;
                details.I_nu_max(i, j, k) =
                    photons(i, j, k).I_nu_max / (unit::erg / (unit::Hz * unit::sec * unit::cm2));
                details.Y_T(i, j, k) = photons(i, j, k).Ys.Y_T;
            }
        }
    }
}

void PyModel::single_evo_details(Shock const& shock, Coord const& coord, Observer& obs, PyRadiation const& rad,
                                 PyShock& details) const {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    details.t_obs = obs.time / unit::sec;
    details.Doppler = xt::exp2(obs.lg2_doppler);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.ssc_cooling) {
        if (rad.kn) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }
    save_electron_details(syn_e, details);
    save_photon_details(syn_ph, details);
}

auto PyModel::details(Real t_min, Real t_max) const -> PyDetails {
    const Array t_obs = xt::logspace(std::log10(t_min * unit::sec), std::log10(t_max * unit::sec), 10);
    Coord coord = auto_grid(jet_, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    PyDetails details;

    details.phi = coord.phi;
    details.theta = coord.theta;
    details.t_src = coord.t / unit::sec;

    Observer observer;

    if (!rvs_rad_opt) {
        const auto fwd_shock = generate_fwd_shock(coord, medium_, jet_, fwd_rad.rad, rtol);

        save_shock_details(fwd_shock, details.fwd);

        single_evo_details(fwd_shock, coord, observer, fwd_rad, details.fwd);

        return details;
    } else {
        const auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium_, jet_, fwd_rad.rad, rvs_rad.rad, rtol);

        save_shock_details(fwd_shock, details.fwd);

        save_shock_details(rvs_shock, details.rvs);

        single_evo_details(fwd_shock, coord, observer, fwd_rad, details.fwd);

        single_evo_details(rvs_shock, coord, observer, rvs_rad, details.rvs);

        return details;
    }
}

void PyFlux::calc_total() {
    total = xt::zeros<Real>(fwd.sync.shape());
    if (fwd.sync.size() > 0) {
        total += fwd.sync;
    }
    if (fwd.ssc.size() > 0) {
        total += fwd.ssc;
    }
    if (rvs.sync.size() > 0) {
        total += rvs.sync;
    }
    if (rvs.ssc.size() > 0) {
        total += rvs.ssc;
    }
}

auto PyModel::flux_density(PyArray const& t, PyArray const& nu) -> PyFlux {
    AFTERGLOW_REQUIRE(
        t.size() == nu.size(),
        "time and frequency arrays must have the same size\nIf you intend to get grid-like output, use the "
        "generic `flux_density_grid` instead");
    AFTERGLOW_REQUIRE(is_ascending(t), "time array must be in ascending order");

    const Array t_obs = t * unit::sec;
    const Array nu_obs = nu * unit::Hz;

    auto flux_func = [](Observer& obs, Array const& time, Array const& freq, auto& photons) -> XTArray {
        return obs.specific_flux_series(time, freq, photons) / unit::flux_den_cgs;
    };

    auto result = compute_emission(t_obs, nu_obs, flux_func);
    result.calc_total();
    return result;
}

auto PyModel::flux(PyArray const& t, double nu_min, double nu_max, size_t num_nu) -> PyFlux {
    AFTERGLOW_REQUIRE(is_ascending(t), "time array must be in ascending order");

    // Generate frequency array
    const Array nu_obs = xt::logspace(std::log10(nu_min * unit::Hz), std::log10(nu_max * unit::Hz), num_nu);
    const Array t_obs = t * unit::sec;

    auto flux_func = [](Observer& obs, Array const& time, Array const& freq, auto& photons) -> XTArray {
        return obs.flux(time, freq, photons) / unit::flux_cgs;
    };

    auto result = compute_emission(t_obs, nu_obs, flux_func);
    result.calc_total();
    return result;
}

auto PyModel::generate_exposure_sampling(PyArray const& t, PyArray const& nu, PyArray const& expo_time,
                                         size_t num_points) -> ExposureSampling {
    const size_t total_points = t.size() * num_points;
    Array t_obs = Array::from_shape({total_points});
    Array nu_obs = Array::from_shape({total_points});
    std::vector<size_t> idx(total_points);

    // Generate time-frequency samples within each exposure window
    for (size_t i = 0, j = 0; i < t.size() && j < total_points; ++i) {
        const Real t_start = t(i);
        const Real dt = expo_time(i) / static_cast<Real>(num_points - 1);

        for (size_t k = 0; k < num_points && j < total_points; ++k, ++j) {
            t_obs(j) = t_start + k * dt;
            nu_obs(j) = nu(i);
            idx[j] = i;
        }
    }

    std::vector<size_t> sort_indices(total_points);
    std::iota(sort_indices.begin(), sort_indices.end(), 0);
    std::ranges::sort(sort_indices, [&t_obs](size_t i, size_t j) { return t_obs(i) < t_obs(j); });

    Array t_obs_sorted = Array::from_shape({total_points});
    Array nu_obs_sorted = Array::from_shape({total_points});
    std::vector<size_t> idx_sorted(idx.size());

    for (size_t i = 0; i < sort_indices.size(); ++i) {
        const size_t orig_idx = sort_indices[i];
        t_obs_sorted(i) = t_obs(orig_idx);
        nu_obs_sorted(i) = nu_obs(orig_idx);
        idx_sorted[i] = idx[orig_idx];
    }

    t_obs_sorted *= unit::sec;
    nu_obs_sorted *= unit::Hz;

    return {std::move(t_obs_sorted), std::move(nu_obs_sorted), std::move(idx_sorted)};
}

void PyModel::average_exposure_flux(PyFlux& result, std::vector<size_t> const& idx_sorted, size_t original_size,
                                    size_t num_points) {
    auto average_component = [&](XTArray& component) {
        if (component.size() > 0) {
            Array summed = xt::zeros<Real>({original_size});
            for (size_t j = 0; j < component.size(); j++) {
                const size_t orig_time_idx = idx_sorted[j];
                summed(orig_time_idx) += component(j);
            }
            summed /= static_cast<Real>(num_points);
            component = std::move(summed);
        }
    };

    average_component(result.fwd.sync);
    average_component(result.fwd.ssc);
    average_component(result.rvs.sync);
    average_component(result.rvs.ssc);
}

auto PyModel::flux_density_exposures(PyArray const& t, PyArray const& nu, PyArray const& expo_time, size_t num_points)
    -> PyFlux {
    AFTERGLOW_REQUIRE(t.size() == nu.size() && t.size() == expo_time.size(),
                      "time, frequency, and exposure time arrays must have the same size");
    AFTERGLOW_REQUIRE(num_points >= 2, "num_points must be at least 2 to sample within each exposure time");

    const auto [t_obs_sorted, nu_obs_sorted, idx_sorted] = generate_exposure_sampling(t, nu, expo_time, num_points);

    auto flux_func = [](Observer& obs, Array const& time, Array const& freq, auto& photons) -> XTArray {
        return obs.specific_flux_series(time, freq, photons) / unit::flux_den_cgs;
    };

    auto result = compute_emission(t_obs_sorted, nu_obs_sorted, flux_func);

    average_exposure_flux(result, idx_sorted, t.size(), num_points);

    result.calc_total();
    return result;
}

auto PyModel::flux_density_grid(PyArray const& t, PyArray const& nu) -> PyFlux {
    AFTERGLOW_REQUIRE(is_ascending(t), "time array must be in ascending order");

    const Array t_obs = t * unit::sec;
    const Array nu_obs = nu * unit::Hz;

    auto flux_func = [](Observer& obs, Array const& time, Array const& freq, auto& photons) -> XTArray {
        return obs.specific_flux(time, freq, photons) / unit::flux_den_cgs;
    };

    auto result = compute_emission(t_obs, nu_obs, flux_func);
    result.calc_total();
    return result;
}

Array PyModel::jet_E_iso(Real phi, Array const& theta) const {
    Array E_iso = xt::zeros<Real>(theta.shape());
    for (size_t i = 0; i < theta.size(); ++i) {
        E_iso(i) = jet_.eps_k(phi, theta(i)) / (unit::erg / (4 * con::pi));
    }
    return E_iso;
}

Array PyModel::jet_Gamma0(Real phi, Array const& theta) const {
    Array Gamma0 = xt::zeros<Real>(theta.shape());
    for (size_t i = 0; i < theta.size(); ++i) {
        Gamma0(i) = jet_.Gamma0(phi, theta(i));
    }
    return Gamma0;
}

Array PyModel::medium(Real phi, Real theta, Array const& r) const {
    Array rho = xt::zeros<Real>(r.shape());
    for (size_t i = 0; i < r.size(); ++i) {
        rho(i) = medium_.rho(phi, theta, r(i) * unit::cm) / (unit::g / unit::cm3);
    }
    return rho;
}
/*
Array PyModel::medium(Real phi, Real theta, Array const& r) {
    Array rho = xt::zeros<Real>(r.shape());
    for (size_t i = 0; i < r.size(); ++i) {
        rho(i) = medium.rho(phi, theta, r(i) * unit::cm) / (unit::g / unit::cm3);
    }
    return rho;
}

Array PyModel::jet_E_iso(Real phi, Array const& theta) {
    Array E_iso = xt::zeros<Real>(theta.shape());
    for (size_t i = 0; i < theta.size(); ++i) {
        E_iso(i) = jet.eps_k(phi, theta(i)) / (unit::erg / (4 * con::pi));
    }
    return E_iso;
}

Array PyModel::jet_Gamma0(Real phi, Array const& theta) {
    Array Gamma0 = xt::zeros<Real>(theta.shape());
    for (size_t i = 0; i < theta.size(); ++i) {
        Gamma0(i) = jet.Gamma0(phi, theta(i));
    }
    return Gamma0;
}*/
