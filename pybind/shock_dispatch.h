//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include "../include/afterglow.h"

/// Solve forward shock dynamics with variant dispatch.
/// Handles std::visit over JetVariant × MediumVariant to invoke auto_grid + generate_fwd_shock
/// with concrete types for optimized template instantiation.
inline auto solve_fwd_shock(JetVariant const& jet, MediumVariant const& medium, Array const& t_obs, Real theta_w,
                            Real theta_obs, Real z, Real phi_resol, Real theta_resol, Real t_resol, bool axisymmetric,
                            size_t min_theta_num, RadParams const& rad, Real rtol) -> std::pair<Coord, Shock> {
    return std::visit(
        [&](auto const& j, auto const& med) {
            auto coord = auto_grid(j, med, t_obs, theta_w, theta_obs, z, phi_resol, theta_resol, t_resol, axisymmetric,
                                   0, min_theta_num);
            auto shock = generate_fwd_shock(coord, med, j, rad, rtol);
            return std::pair{std::move(coord), std::move(shock)};
        },
        jet, medium);
}

/// Solve forward + reverse shock pair with variant dispatch.
inline auto solve_shock_pair(JetVariant const& jet, MediumVariant const& medium, Array const& t_obs, Real theta_w,
                             Real theta_obs, Real z, Real phi_resol, Real theta_resol, Real t_resol, bool axisymmetric,
                             size_t min_theta_num, RadParams const& fwd_rad, RadParams const& rvs_rad, Real rtol)
    -> std::tuple<Coord, Shock, Shock> {
    return std::visit(
        [&](auto const& j, auto const& med) {
            auto coord = auto_grid(j, med, t_obs, theta_w, theta_obs, z, phi_resol, theta_resol, t_resol, axisymmetric,
                                   0, min_theta_num);
            auto [fwd, rvs] = generate_shock_pair(coord, med, j, fwd_rad, rvs_rad, rtol);
            return std::tuple{std::move(coord), std::move(fwd), std::move(rvs)};
        },
        jet, medium);
}
