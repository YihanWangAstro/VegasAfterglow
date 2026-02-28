//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include <algorithm>
#include <numeric>

#include "error_handling.h"

//========================================================================================================
//                                  SpectrumEvaluator / YEvaluator
//========================================================================================================

XTArray SpectrumEvaluator::operator()(PyArray const& nu_comv) const {
    const size_t n = nu_comv.size();
    XTArray result = xt::zeros<Real>({n});
    for (size_t i = 0; i < n; ++i) {
        result(i) = eval_(nu_comv(i));
    }
    return result;
}

XTArray YEvaluator::operator()(PyArray const& gamma) const {
    const size_t n = gamma.size();
    XTArray result = xt::zeros<Real>({n});
    for (size_t i = 0; i < n; ++i) {
        result(i) = eval_(gamma(i));
    }
    return result;
}

// Overload for magnetar
void initialize_ejecta(Ejecta& jet, bool spreading, Real duration, std::optional<PyMagnetar> const& magnetar,
                       Real theta_c) {
    jet.spreading = spreading;
    jet.T0 = duration;
    if (magnetar) {
        jet.deps_dt = math::magnetar_injection(magnetar->t0, magnetar->q, magnetar->L0, theta_c);
    }
}

JetVariant PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                       std::optional<PyMagnetar> const& magnetar) {
    if (magnetar) {
        Ejecta jet;
        jet.eps_k = math::tophat(theta_c, E_iso);
        jet.Gamma0 = math::tophat_plus_one(theta_c, Gamma0 - 1);
        initialize_ejecta(jet, spreading, duration, magnetar, theta_c);
        return jet;
    }
    return TophatJet(theta_c, E_iso * unit::erg, Gamma0, spreading, duration * unit::sec);
}

JetVariant PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                         std::optional<PyMagnetar> const& magnetar) {
    if (magnetar) {
        Ejecta jet;
        jet.eps_k = math::gaussian(theta_c, E_iso);
        jet.Gamma0 = math::gaussian_plus_one(theta_c, Gamma0 - 1);
        initialize_ejecta(jet, spreading, duration, magnetar, theta_c);
        return jet;
    }
    return GaussianJet(theta_c, E_iso * unit::erg, Gamma0, spreading, duration * unit::sec);
}

JetVariant PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k_e, Real k_g, bool spreading, Real duration,
                         std::optional<PyMagnetar> const& magnetar) {
    if (magnetar) {
        Ejecta jet;
        jet.eps_k = math::powerlaw(theta_c, E_iso, k_e);
        jet.Gamma0 = math::powerlaw_plus_one(theta_c, Gamma0 - 1, k_g);
        initialize_ejecta(jet, spreading, duration, magnetar, theta_c);
        return jet;
    }
    return PowerLawJet(theta_c, E_iso * unit::erg, Gamma0, k_e, k_g, spreading, duration * unit::sec);
}

JetVariant PyPowerLawWing(Real theta_c, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g, bool spreading,
                          Real duration) {
    Ejecta jet;
    jet.eps_k = math::powerlaw_wing(theta_c, E_iso_w, k_e);
    jet.Gamma0 = math::powerlaw_wing_plus_one(theta_c, Gamma0_w - 1, k_g);
    jet.spreading = spreading;
    jet.T0 = duration;
    return jet;
}

JetVariant PyStepPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real E_iso_w, Real Gamma0_w, Real k_e, Real k_g,
                             bool spreading, Real duration, std::optional<PyMagnetar> const& magnetar) {
    Ejecta jet;
    jet.eps_k = math::step_powerlaw(theta_c, E_iso, E_iso_w, k_e);
    jet.Gamma0 = math::step_powerlaw_plus_one(theta_c, Gamma0 - 1, Gamma0_w - 1, k_g);
    initialize_ejecta(jet, spreading, duration, magnetar, theta_c);
    return jet;
}

JetVariant PyTwoComponentJet(Real theta_c, Real E_iso, Real Gamma0, Real theta_w, Real E_iso_w, Real Gamma0_w,
                             bool spreading, Real duration, std::optional<PyMagnetar> const& magnetar) {
    Ejecta jet;
    jet.eps_k = math::two_component(theta_c, theta_w, E_iso, E_iso_w);
    jet.Gamma0 = math::two_component_plus_one(theta_c, theta_w, Gamma0 - 1, Gamma0_w - 1);
    initialize_ejecta(jet, spreading, duration, magnetar, theta_c);
    return jet;
}

ISM PyISM(Real n_ism) {
    return ISM(n_ism / unit::cm3);
}

MediumVariant PyWind(Real A_star, std::optional<Real> n_ism_opt, std::optional<Real> n0_opt, Real k_m) {
    const Real n_ism = n_ism_opt.value_or(0);
    const Real n0 = n0_opt.value_or(con::inf);

    if (k_m == 2) {
        return Wind(A_star, n_ism / unit::cm3, n0 / unit::cm3);
    }
    // General k_m: fall back to Medium with std::function (evn::wind handles unit conversion)
    Medium medium;
    medium.rho = evn::wind(A_star, n_ism / unit::cm3, n0 / unit::cm3, k_m);
    medium.isotropic = true;
    return medium;
}

void convert_unit_jet(JetVariant& jet) {
    std::visit(
        [](auto& j) {
            if constexpr (std::is_same_v<std::decay_t<decltype(j)>, Ejecta>) {
                const auto eps_k_cgs = j.eps_k;
                j.eps_k = [=](Real phi, Real theta) { return eps_k_cgs(phi, theta) * (unit::erg / (4 * con::pi)); };

                const auto deps_dt_cgs = j.deps_dt;
                j.deps_dt = [=](Real phi, Real theta, Real t) {
                    return deps_dt_cgs(phi, theta, t / unit::sec) * (unit::erg / (4 * con::pi * unit::sec));
                };

                const auto dm_dt_cgs = j.dm_dt;
                j.dm_dt = [=](Real phi, Real theta, Real t) {
                    return dm_dt_cgs(phi, theta, t / unit::sec) * (unit::g / (4 * con::pi * unit::sec));
                };

                j.T0 *= unit::sec;
            }
            // TophatJet, GaussianJet, PowerLawJet are already constructed in internal units
        },
        jet);
}

void convert_unit_medium(MediumVariant& medium) {
    std::visit(
        [](auto& m) {
            if constexpr (std::is_same_v<std::decay_t<decltype(m)>, Medium>) {
                const auto rho_cgs = m.rho;
                m.rho = [=](Real phi, Real theta, Real r) {
                    return rho_cgs(phi, theta, r / unit::cm) * (unit::g / unit::cm3);
                };
            }
            // ISM and Wind are already constructed in internal units — no conversion needed
        },
        medium);
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
void save_photon_details(PhotonGrid const& photons, PyShock& details, Shock const& shock) {
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
                details.nu_m_hat(i, j, k) =
                    compute_syn_freq(photons(i, j, k).Ys.gamma_m_hat, shock.B(i, j, k)) / unit::Hz;
                details.nu_c_hat(i, j, k) =
                    compute_syn_freq(photons(i, j, k).Ys.gamma_c_hat, shock.B(i, j, k)) / unit::Hz;
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

    auto syn_e = generate_syn_electrons(shock, coord);

    auto syn_ph = generate_syn_photons(shock, syn_e, coord);

    if (rad.ssc) {
        if (rad.kn) {
            KN_cooling(syn_e, syn_ph, shock, coord, obs_setup.z);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock, coord, obs_setup.z);
        }
    } else if (rad.rad.cmb_cooling) {
        CMB_cooling(syn_e, syn_ph, shock, coord, obs_setup.z);
    }
    save_electron_details(syn_e, details);
    save_photon_details(syn_ph, details, shock);

    // Store photon grids for per-cell spectrum evaluation
    details.syn_photons_ = syn_ph;
    details.has_syn_spectrum_ = true;

    if (rad.ssc) {
        details.ic_photons_ = generate_IC_photons(syn_e, syn_ph, rad.kn, coord);
        details.has_ssc_spectrum_ = true;
    }
}

auto PyModel::details(Real t_min, Real t_max) const -> PyDetails {
    const Array t_obs = xt::logspace(std::log10(t_min * unit::sec), std::log10(t_max * unit::sec), 10);

    PyDetails details;
    Observer observer;

    if (!rvs_rad_opt) {
        auto [coord, fwd_shock] = solve_fwd_shock(jet_, medium_, t_obs, grid_config(), fwd_rad.rad, rtol);

        details.phi = coord.phi;
        details.theta = coord.theta;
        details.t_src = coord.t / unit::sec;
        save_shock_details(fwd_shock, details.fwd);
        single_evo_details(fwd_shock, coord, observer, fwd_rad, details.fwd);
    } else {
        auto [coord, fwd_shock, rvs_shock] =
            solve_shock_pair(jet_, medium_, t_obs, grid_config(), fwd_rad.rad, rvs_rad_opt->rad, rtol);

        details.phi = coord.phi;
        details.theta = coord.theta;
        details.t_src = coord.t / unit::sec;
        save_shock_details(fwd_shock, details.fwd);
        save_shock_details(rvs_shock, details.rvs);
        single_evo_details(fwd_shock, coord, observer, fwd_rad, details.fwd);
        single_evo_details(rvs_shock, coord, observer, *rvs_rad_opt, details.rvs);
    }
    return details;
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
        E_iso(i) = jet_eps_k(jet_, phi, theta(i)) / (unit::erg / (4 * con::pi));
    }
    return E_iso;
}

Array PyModel::jet_Gamma0(Real phi, Array const& theta) const {
    Array Gamma0 = xt::zeros<Real>(theta.shape());
    for (size_t i = 0; i < theta.size(); ++i) {
        Gamma0(i) = ::jet_Gamma0(jet_, phi, theta(i));
    }
    return Gamma0;
}

Array PyModel::medium(Real phi, Real theta, Array const& r) const {
    Array rho = xt::zeros<Real>(r.shape());
    for (size_t i = 0; i < r.size(); ++i) {
        rho(i) = medium_rho(medium_, phi, theta, r(i) * unit::cm) / (unit::g / unit::cm3);
    }
    return rho;
}
