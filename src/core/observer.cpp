//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "observer.h"

#include "../util/macros.h"

//========================================================================================================
//                                  File-Local Helper Functions
//========================================================================================================

static Array compute_dphi(Coord const& coord, size_t eff_phi_grid) {
    Array dphi({eff_phi_grid}, 0);
    if (eff_phi_grid == 1) {
        dphi(0) = 2 * con::pi;
    } else {
        const size_t last = eff_phi_grid - 1;
        for (size_t i = 0; i < eff_phi_grid; ++i) {
            dphi(i) = 0.5 * (coord.phi(std::min(i + 1, last)) - coord.phi(i > 0 ? i - 1 : size_t(0)));
        }
    }
    return dphi;
}

//========================================================================================================
//                                  Static Member Functions
//========================================================================================================

Real Observer::loglog_interpolate(InterpState const& state, Real lg2_t_obs, Real lg2_t_lo) noexcept {
    const Real dlg2_t = lg2_t_obs - lg2_t_lo;
    return fast_exp2(state.lg2_L_nu_lo + dlg2_t * state.slope);
}

//========================================================================================================
//                                  Private Methods - Basic Calculations
//========================================================================================================

void Observer::calc_t_obs(Coord const& coord, Shock const& shock) {
    const Real cos_obs = std::cos(coord.theta_view);
    const Real sin_obs = std::sin(coord.theta_view);

    // Spreading jets: theta varies along k, need full 3D trig arrays
    static thread_local MeshGrid3d cos_theta, sin_theta;
    if (cos_theta.shape() != shock.theta.shape()) {
        cos_theta = MeshGrid3d::from_shape(shock.theta.shape());
        sin_theta = MeshGrid3d::from_shape(shock.theta.shape());
    }

    const size_t shock_phi_size = shock.theta.shape(0);
    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_grid; ++j) {
            for (size_t k = 0; k < t_grid; ++k) {
                cos_theta(i, j, k) = std::cos(shock.theta(i, j, k));
                sin_theta(i, j, k) = std::sin(shock.theta(i, j, k));
            }
        }
    }

    for (size_t i = 0; i < eff_phi_grid; ++i) {
        const Real cos_phi = std::cos(coord.phi[i] - coord.phi_view);
        const size_t i_eff = i * jet_3d;

        for (size_t j = 0; j < theta_grid; ++j) {
            for (size_t k = 0; k < t_grid; ++k) {
                const Real gamma_ = shock.Gamma(i_eff, j, k);
                const Real r = shock.r(i_eff, j, k);
                const Real t_eng_ = coord.t(i_eff, j, k);
                const Real cos_v = sin_theta(i_eff, j, k) * cos_phi * sin_obs + cos_theta(i_eff, j, k) * cos_obs;
                const Real t_val = (t_eng_ + (1 - cos_v) * r / con::c) * one_plus_z;

                lg2_doppler(i, j, k) = -fast_log2(gamma_ - std::sqrt(gamma_ * gamma_ - 1) * cos_v);
                time(i, j, k) = t_val;
                lg2_t(i, j, k) = fast_log2(t_val);
            }
        }
    }
}

void Observer::calc_solid_angle(Coord const& coord, Shock const& shock) {
    const Array dphi = compute_dphi(coord, eff_phi_grid);
    const size_t last = theta_grid - 1;
    const size_t shock_phi_size = shock.theta.shape(0);

    // Spreading jets: dcos varies along k, need full 3D storage
    static thread_local MeshGrid3d dcos;
    if (dcos.shape() != shock.theta.shape()) {
        dcos = MeshGrid3d::from_shape(shock.theta.shape());
    }

    // Interpolate theta of neighbor cell (i, j_nb) at engine time t_target
    auto interp_theta = [&](size_t i, size_t j_nb, Real t_target, size_t& k_hint) -> Real {
        while (k_hint + 1 < t_grid && coord.t(i, j_nb, k_hint + 1) < t_target)
            k_hint++;
        if (k_hint + 1 >= t_grid)
            return shock.theta(i, j_nb, t_grid - 1);
        const Real w =
            (t_target - coord.t(i, j_nb, k_hint)) / (coord.t(i, j_nb, k_hint + 1) - coord.t(i, j_nb, k_hint));
        return shock.theta(i, j_nb, k_hint) + w * (shock.theta(i, j_nb, k_hint + 1) - shock.theta(i, j_nb, k_hint));
    };

    for (size_t i = 0; i < shock_phi_size; ++i) {
        for (size_t j = 0; j < theta_grid; ++j) {
            const size_t j_p1 = (j == last) ? last : (j + 1);
            size_t k_hint_lo = 0, k_hint_hi = 0;
            for (size_t k = 0; k < t_grid; ++k) {
                const Real t_target = coord.t(i, j, k);
                const Real th_lo =
                    (j == 0) ? 0.0 : 0.5 * (shock.theta(i, j, k) + interp_theta(i, j - 1, t_target, k_hint_lo));
                const Real th_hi = (j == last)
                                       ? shock.theta(i, j, k)
                                       : 0.5 * (shock.theta(i, j, k) + interp_theta(i, j_p1, t_target, k_hint_hi));
                dcos(i, j, k) = std::cos(th_hi) - std::cos(th_lo);
            }
        }
    }

    for (size_t i = 0; i < eff_phi_grid; ++i) {
        const size_t i_eff = i * jet_3d;
        for (size_t j = 0; j < theta_grid; ++j) {
            for (size_t k = 0; k < t_grid; ++k) {
                const Real dOmega = std::fabs(dcos(i_eff, j, k) * dphi(i));
                const Real r = shock.r(i_eff, j, k);
                lg2_geom_factor(i, j, k) = fast_log2(dOmega * r * r) + 3 * lg2_doppler(i, j, k);
            }
        }
    }
}

void Observer::calc_eat_non_spreading(Coord const& coord, Shock const& shock) {
    const Real cos_obs = std::cos(coord.theta_view);
    const Real sin_obs = std::sin(coord.theta_view);
    const Array dphi = compute_dphi(coord, eff_phi_grid);
    const size_t last = theta_grid - 1;

    for (size_t i = 0; i < eff_phi_grid; ++i) {
        const Real cos_phi = std::cos(coord.phi[i] - coord.phi_view);
        const size_t i_eff = i * jet_3d;
        const Real dphi_i = dphi(i);

        for (size_t j = 0; j < theta_grid; ++j) {
            const Real ct = std::cos(shock.theta(i_eff, j, 0));
            const Real st = std::sin(shock.theta(i_eff, j, 0));
            const Real cos_v = st * cos_phi * sin_obs + ct * cos_obs;
            const Real t_coeff = (1 - cos_v) / con::c * one_plus_z;

            // Solid angle: constant along k for non-spreading jets
            const size_t j_p1 = (j == last) ? last : (j + 1);
            const Real th_lo = (j == 0) ? 0.0 : 0.5 * (shock.theta(i_eff, j, 0) + shock.theta(i_eff, j - 1, 0));
            const Real th_hi = 0.5 * (shock.theta(i_eff, j, 0) + shock.theta(i_eff, j_p1, 0));
            const Real dOmega = std::fabs((std::cos(th_hi) - std::cos(th_lo)) * dphi_i);
            const Real lg2_dOmega = fast_log2(dOmega);

            for (size_t k = 0; k < t_grid; ++k) {
                const Real gamma_ = shock.Gamma(i_eff, j, k);
                const Real r = shock.r(i_eff, j, k);
                const Real t_val = coord.t(i_eff, j, k) * one_plus_z + t_coeff * r;
                const Real lg2_dop = -fast_log2(gamma_ - std::sqrt(gamma_ * gamma_ - 1) * cos_v);

                lg2_doppler(i, j, k) = lg2_dop;
                time(i, j, k) = t_val;
                lg2_t(i, j, k) = fast_log2(t_val);
                lg2_geom_factor(i, j, k) = lg2_dOmega + fast_log2(r * r) + 3 * lg2_dop;
            }
        }
    }
}

//========================================================================================================
//                                  Private Methods - Initialization
//========================================================================================================

void Observer::build_time_grid(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    auto [phi_size, theta_size, t_size] = shock.shape();

    // Determine if the jet is 3D (more than one phi value)
    jet_3d = static_cast<size_t>((phi_size > 1));

    // Set the effective phi grid size based on the observation angle and jet dimensionality.
    if (coord.theta_view == 0 && jet_3d == 0) {
        eff_phi_grid = 1; // optimize for on-axis observer
    } else {
        eff_phi_grid = coord.phi.size();
    }

    theta_grid = theta_size;
    t_grid = t_size;
    lumi_dist = luminosity_dist;
    one_plus_z = 1 + redshift;

    // Detect jet spreading: theta varies along the k (time) axis
    jet_spreading_ = false;
    if (t_size > 1) {
        for (size_t i = 0; i < phi_size && !jet_spreading_; ++i) {
            for (size_t j = 0; j < theta_size && !jet_spreading_; ++j) {
                jet_spreading_ = (shock.theta(i, j, 0) != shock.theta(i, j, t_size - 1));
            }
        }
    }

    time = MeshGrid3d::from_shape({eff_phi_grid, theta_size, t_size});
    lg2_t = MeshGrid3d::from_shape({eff_phi_grid, theta_size, t_size});
    lg2_doppler = MeshGrid3d::from_shape({eff_phi_grid, theta_size, t_size});
    lg2_geom_factor = MeshGrid3d::from_shape({eff_phi_grid, theta_size, t_size});
}

//========================================================================================================
//                                  Public Interface Methods
//========================================================================================================

void Observer::observe(Coord const& coord, Shock const& shock, Real luminosity_dist, Real redshift) {
    build_time_grid(coord, shock, luminosity_dist, redshift);
    if (jet_spreading_) {
        calc_t_obs(coord, shock);
        calc_solid_angle(coord, shock);
    } else {
        calc_eat_non_spreading(coord, shock);
    }
}
