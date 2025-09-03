//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "pymodel.h"

#include <algorithm>
#include <numeric>

#include "afterglow.h"

Ejecta PyTophatJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                   std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = [=](Real phi, Real theta) { return theta < theta_c ? E_iso : 0.; };
    jet.Gamma0 = [=](Real phi, Real theta) { return theta < theta_c ? Gamma0 : 1.; };
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 * std::pow(tt, -magnetar->q);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Ejecta PyGaussianJet(Real theta_c, Real E_iso, Real Gamma0, bool spreading, Real duration,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    Real norm = -1 / (2 * theta_c * theta_c);
    jet.eps_k = [=](Real phi, Real theta) { return E_iso * fast_exp(norm * theta * theta); };
    jet.Gamma0 = [=](Real phi, Real theta) { return (Gamma0 - 1) * fast_exp(norm * theta * theta) + 1; };
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 * std::pow(tt, -magnetar->q);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Ejecta PyPowerLawJet(Real theta_c, Real E_iso, Real Gamma0, Real k, bool spreading, Real duration,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = [=](Real phi, Real theta) { return E_iso / (1 + fast_pow(theta / theta_c, k)); };
    jet.Gamma0 = [=](Real phi, Real theta) { return (Gamma0 - 1) / (1 + fast_pow(theta / theta_c, k)) + 1; };
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 * std::pow(tt, -magnetar->q);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Ejecta PyPowerLawGeneralJet(Real theta_c, Real theta_cG, Real E_iso, Real Gamma0, Real k, Real kG, bool spreading, Real duration,
                     std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = [=](Real phi, Real theta) { return E_iso / (1 + fast_pow(theta / theta_c, k)); };
    jet.Gamma0 = [=](Real phi, Real theta) { return (Gamma0 - 1) / (1 + fast_pow(theta / theta_cG, kG)) + 1; };
    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_c) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 * std::pow(tt, -magnetar->q);
            } else {
                return 0.;
            }
        };
    } 

    return jet;
}

Ejecta PyTwoComponentJet(Real theta_n, Real E_iso_n, Real Gamma0_n, Real theta_w, Real E_iso_w, Real Gamma0_w,
                         bool spreading, Real duration, std::optional<PyMagnetar> magnetar) {
    Ejecta jet;
    jet.eps_k = [=](Real phi, Real theta) {
        if (theta <= theta_n) {
            return E_iso_n;
        } else if (theta <= theta_w) {
            return E_iso_w;
        } else {
            return 0.;
        }
    };

    jet.Gamma0 = [=](Real phi, Real theta) {
        if (theta <= theta_n) {
            return Gamma0_n;
        } else if (theta <= theta_w) {
            return Gamma0_w;
        } else {
            return 1.0;  // Default Lorentz factor outside the jet
        }
    };

    jet.spreading = spreading;
    jet.T0 = duration;

    if (magnetar) {
        jet.deps_dt = [=](Real phi, Real theta, Real t) {
            if (theta <= theta_w) {
                Real tt = 1 + t / magnetar->t_0;
                return magnetar->L_0 * std::pow(tt, -magnetar->q);
            } else {
                return 0.;
            }
        };
    }

    return jet;
}

Medium PyISM(Real n_ism) {
    Medium medium;
    medium.rho = [=](Real phi, Real theta, Real r) { return n_ism * 1.67e-24; };
    medium.mass = [=](Real phi, Real theta, Real r) { return n_ism * (1.67e-24 / 3) * r * r * r; };

    return medium;
}

Medium PyWind(Real A_star) {
    Medium medium;
    medium.rho = [=](Real phi, Real theta, Real r) { return A_star * 5e11 / (r * r); };
    medium.mass = [=](Real phi, Real theta, Real r) {
        return A_star * 5e11 * r;  // Integrated mass per solid angle
    };
    return medium;
}

void convert_unit(Ejecta& jet, Medium& medium) {
    auto eps_k_cgs = jet.eps_k;
    jet.eps_k = [=](Real phi, Real theta) { return eps_k_cgs(phi, theta) * (unit::erg / (4 * con::pi)); };

    auto deps_dt_cgs = jet.deps_dt;
    jet.deps_dt = [=](Real phi, Real theta, Real t) {
        return deps_dt_cgs(phi, theta, t / unit::sec) * (unit::erg / (4 * con::pi * unit::sec));
    };

    auto dm_dt_cgs = jet.dm_dt;
    jet.dm_dt = [=](Real phi, Real theta, Real t) {
        return dm_dt_cgs(phi, theta, t / unit::sec) * (unit::g / (4 * con::pi * unit::sec));
    };

    jet.T0 *= unit::sec;

    auto rho_cgs = medium.rho;
    medium.rho = [=](Real phi, Real theta, Real r) {
        return rho_cgs(phi, theta, r / unit::cm) * (unit::g / unit::cm3);
    };

    auto mass_cgs = medium.mass;
    medium.mass = [=](Real phi, Real theta, Real r) { return mass_cgs(phi, theta, r / unit::cm) * unit::g; };
}

void PyModel::single_shock_emission(Shock const& shock, Coord const& coord, Array const& t_obs, Array const& nu_obs,
                                    Observer& obs, PyRadiation rad, ArrayDict& flux_dict, std::string suffix,
                                    bool serilized) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.IC_cooling) {
        if (rad.KN) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }

    if (rad.SSC) {
        auto IC_ph = generate_IC_photons(syn_e, syn_ph, rad.KN);

        if (serilized) {
            flux_dict["IC" + suffix] = obs.specific_flux_series(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        } else {
            flux_dict["IC" + suffix] = obs.specific_flux(t_obs, nu_obs, IC_ph) / unit::flux_den_cgs;
        }
    }

    if (serilized) {
        flux_dict["syn" + suffix] = obs.specific_flux_series(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    } else {
        flux_dict["syn" + suffix] = obs.specific_flux(t_obs, nu_obs, syn_ph) / unit::flux_den_cgs;
    }
}

auto PyModel::compute_specific_flux(Array const& t_obs, Array const& nu_obs, bool serilized) -> ArrayDict {
    Coord coord = auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    ArrayDict flux_dict;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", serilized);

        return flux_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium, jet, fwd_rad.rad, rvs_rad.rad, rtol);

        single_shock_emission(fwd_shock, coord, t_obs, nu_obs, observer, fwd_rad, flux_dict, "", serilized);

        single_shock_emission(rvs_shock, coord, t_obs, nu_obs, observer, rvs_rad, flux_dict, "_rvs", serilized);

        return flux_dict;
    }
}

void save_shock_details(Shock const& shock, ArrayDict& detail_dict, std::string suffix) {
    detail_dict["Gamma_downstr" + suffix] = shock.Gamma;
    detail_dict["Gamma_rel" + suffix] = shock.Gamma_rel;
    detail_dict["r" + suffix] = shock.r / unit::cm;
    detail_dict["t_comv" + suffix] = shock.t_comv / unit::sec;
    detail_dict["B" + suffix] = shock.B / unit::Gauss;
    detail_dict["N_p" + suffix] = shock.proton_num;
    detail_dict["theta" + suffix] = shock.theta;
}

template <typename ElectronGrid>
void save_electron_details(ElectronGrid const& electrons, ArrayDict& detail_dict, std::string suffix) {
    auto shape = electrons.shape();
    detail_dict["gamma_m" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_c" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_a" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["gamma_M" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["N_e" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                detail_dict["gamma_a" + suffix](i, j, k) = electrons(i, j, k).gamma_a;
                detail_dict["gamma_m" + suffix](i, j, k) = electrons(i, j, k).gamma_m;
                detail_dict["gamma_c" + suffix](i, j, k) = electrons(i, j, k).gamma_c;
                detail_dict["gamma_M" + suffix](i, j, k) = electrons(i, j, k).gamma_M;
                detail_dict["N_e" + suffix](i, j, k) = electrons(i, j, k).N_e;
            }
        }
    }
}
template <typename PhotonGrid>
void save_photon_details(PhotonGrid const& photons, ArrayDict& detail_dict, std::string suffix) {
    auto shape = photons.shape();
    detail_dict["nu_m" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_c" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_a" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["nu_M" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    detail_dict["I_nu_max" + suffix] = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                detail_dict["nu_a" + suffix](i, j, k) = photons(i, j, k).nu_a;
                detail_dict["nu_m" + suffix](i, j, k) = photons(i, j, k).nu_m;
                detail_dict["nu_c" + suffix](i, j, k) = photons(i, j, k).nu_c;
                detail_dict["nu_M" + suffix](i, j, k) = photons(i, j, k).nu_M;
                detail_dict["I_nu_max" + suffix](i, j, k) = photons(i, j, k).I_nu_max;
            }
        }
    }
    detail_dict["nu_m" + suffix] /= unit::Hz;
    detail_dict["nu_c" + suffix] /= unit::Hz;
    detail_dict["nu_a" + suffix] /= unit::Hz;
    detail_dict["nu_M" + suffix] /= unit::Hz;
    detail_dict["I_nu_max" + suffix] /= (unit::erg / (unit::Hz * unit::sec * unit::cm2));
}

void PyModel::single_shock_details(Shock const& shock, Coord const& coord, Array const& t_obs, Observer& obs,
                                   PyRadiation rad, ArrayDict& detail_dict, std::string suffix) {
    obs.observe(coord, shock, obs_setup.lumi_dist, obs_setup.z);

    detail_dict["EAT" + suffix] = obs.time / unit::sec;
    detail_dict["Doppler" + suffix] = xt::exp2(obs.lg2_doppler);
    // detail_dict["Omega" + suffix] = xt::exp2(obs.lg2_emission_area);

    auto syn_e = generate_syn_electrons(shock);

    auto syn_ph = generate_syn_photons(shock, syn_e);

    if (rad.IC_cooling) {
        if (rad.KN) {
            KN_cooling(syn_e, syn_ph, shock);
        } else {
            Thomson_cooling(syn_e, syn_ph, shock);
        }
    }
    save_electron_details(syn_e, detail_dict, suffix);
    save_photon_details(syn_ph, detail_dict, suffix);
}

auto PyModel::details(PyArray const& t_obs) -> ArrayDict {
    Coord coord = auto_grid(jet, t_obs, this->theta_w, obs_setup.theta_obs, obs_setup.z, phi_resol, theta_resol,
                            t_resol, axisymmetric);

    ArrayDict details_dict;

    details_dict["phi"] = coord.phi;
    details_dict["theta"] = coord.theta;
    details_dict["t_src"] = coord.t / unit::sec;

    Observer observer;

    if (!rvs_rad_opt) {
        auto fwd_shock = generate_fwd_shock(coord, medium, jet, fwd_rad.rad, rtol);

        save_shock_details(fwd_shock, details_dict, "_fwd");

        single_shock_details(fwd_shock, coord, t_obs, observer, fwd_rad, details_dict, "_fwd");

        return details_dict;
    } else {
        auto rvs_rad = *rvs_rad_opt;
        auto [fwd_shock, rvs_shock] = generate_shock_pair(coord, medium, jet, fwd_rad.rad, rvs_rad.rad, rtol);

        save_shock_details(fwd_shock, details_dict, "_fwd");

        save_shock_details(rvs_shock, details_dict, "_rvs");

        single_shock_details(fwd_shock, coord, t_obs, observer, fwd_rad, details_dict, "_fwd");

        single_shock_details(rvs_shock, coord, t_obs, observer, rvs_rad, details_dict, "_rvs");

        return details_dict;
    }
}

auto PyModel::specific_flux_sorted_series(PyArray const& t, PyArray const& nu) -> ArrayDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool serilized = true;

    if (t_obs.size() != nu_obs.size()) {
        throw std::invalid_argument(
            "time and frequency arrays must have the same size\n"
            "If you intend to get matrix-like output, use the generic `specific_flux` instead");
    }

    return compute_specific_flux(t_obs, nu_obs, serilized);
}

auto PyModel::specific_flux_series(PyArray const& t, PyArray const& nu) -> ArrayDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool serilized = true;

    if (t_obs.size() != nu_obs.size()) {
        throw std::invalid_argument(
            "time and frequency arrays must have the same size\n"
            "If you intend to get matrix-like output, use the generic `specific_flux` instead");
    }

    // Create sorted indices to handle random order
    std::vector<size_t> sorted_indices(t_obs.size());
    std::iota(sorted_indices.begin(), sorted_indices.end(), 0);
    std::sort(sorted_indices.begin(), sorted_indices.end(),
              [&t_obs](size_t i, size_t j) { return t_obs(i) < t_obs(j); });

    // Create sorted arrays
    Array t_sorted = xt::zeros<Real>({t_obs.size()});
    Array nu_sorted = xt::zeros<Real>({nu_obs.size()});
    for (size_t i = 0; i < sorted_indices.size(); ++i) {
        t_sorted(i) = t_obs(sorted_indices[i]);
        nu_sorted(i) = nu_obs(sorted_indices[i]);
    }

    // Compute flux with sorted arrays
    ArrayDict sorted_flux_dict = compute_specific_flux(t_sorted, nu_sorted, serilized);

    // Reorder results back to original order
    ArrayDict flux_dict;
    for (auto const& [key, sorted_flux] : sorted_flux_dict) {
        Array reordered_flux = xt::zeros<Real>({sorted_flux.size()});
        for (size_t i = 0; i < sorted_indices.size(); ++i) {
            reordered_flux(sorted_indices[i]) = sorted_flux(i);
        }
        flux_dict[key] = reordered_flux;
    }

    return flux_dict;
}

auto PyModel::specific_flux(PyArray const& t, PyArray const& nu) -> ArrayDict {
    Array t_obs = t * unit::sec;
    Array nu_obs = nu * unit::Hz;
    bool return_trace = false;

    return compute_specific_flux(t_obs, nu_obs, return_trace);
}
