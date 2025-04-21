//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "observer.h"

#include <cmath>

#include "macros.h"
#include "physics.h"
#include "utilities.h"

/********************************************************************************************************************
 * METHOD: LogScaleInterp::interpLuminosity
 * DESCRIPTION: Interpolates the luminosity at a given observation time (t) using linear
 *              interpolation in log-space. The result is returned in linear space.
 ********************************************************************************************************************/

Real LogScaleInterp::interpLuminosity(Real t_obs) const {
    Real log_t = fastLog(t_obs / t_obs_lo);
    return L_lo * fastExp(log_t * log_L_ratio / log_t_ratio);
}

/********************************************************************************************************************
 * METHOD: Observer::calcSolidAngle
 * DESCRIPTION: Calculates the solid angle (dOmega) for each effective phi and theta grid point.
 *              The solid angle is computed as the product of the differential cosine of theta and either 2π (if
 *              the effective phi size is 1) or the differential phi value.
 ********************************************************************************************************************/
void Observer::calcSolidAngle(Coord const& coord, MeshGrid3d const& theta_grid) {
    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real phi_lo = 0;
        Real phi_hi = 0;
        if (eff_phi_size == 1) {
            phi_lo = 0;
            phi_hi = 4 * con::pi;
        } else if (i == 0) {  // note this also implys phi.size() > 1
            phi_lo = coord.phi(0);
            phi_hi = coord.phi(1);
        } else if (i == coord.phi.size() - 1) {
            phi_lo = coord.phi(i - 1);
            phi_hi = coord.phi(i);
        } else {
            phi_lo = coord.phi(i - 1);
            phi_hi = coord.phi(i + 1);
        }
        Real dphi = std::abs(phi_hi - phi_lo) / 2;
        size_t i_eff = i * interp.jet_3d;
        for (size_t j = 0; j < theta_size; ++j) {
            for (size_t k = 0; k < t_size; ++k) {
                Real theta_lo = 0;
                Real theta_hi = 0;
                if (j == 0) {
                    theta_lo = theta_grid(i_eff, 0, k);
                    theta_hi = 0.5 * (theta_grid(i_eff, 1, k) + theta_grid(i_eff, 0, k));
                } else if (j == coord.theta.size() - 1) {
                    theta_lo = 0.5 * (theta_grid(i_eff, j - 1, k) + theta_grid(i_eff, j, k));
                    theta_hi = theta_grid(i_eff, j, k);
                } else {
                    theta_lo = 0.5 * (theta_grid(i_eff, j, k) + theta_grid(i_eff, j - 1, k));
                    theta_hi = 0.5 * (theta_grid(i_eff, j, k) + theta_grid(i_eff, j + 1, k));
                }
                Real dcos = std::abs(std::cos(theta_hi) - std::cos(theta_lo));
                dOmega(i, j, k) = dcos * dphi;
            }
        }
    }
}

/********************************************************************************************************************
 * METHOD: Observer::calcObsTimeGrid
 * DESCRIPTION: Calculates the observation time grid (t_obs_grid) and updates the doppler factor grid based on the
 *              Gamma (Lorentz factor) and engine time (t) array.
 *              For each grid point, the Doppler factor is computed and the observed time is calculated taking
 *              redshift into account.
 ********************************************************************************************************************/
void Observer::calcObsTimeGrid(Coord const& coord, MeshGrid3d const& Gamma, Real theta_obs) {
    Real cos_obs = std::cos(theta_obs);
    Real sin_obs = std::sin(theta_obs);
    for (size_t i = 0; i < eff_phi_size; ++i) {
        Real cos_phi = std::cos(coord.phi[i]);
        for (size_t j = 0; j < theta_size; ++j) {
            // Compute the cosine of the angle between the local velocity vector and the observer's line of sight.
            Real cos_v = std::sin(coord.theta[j]) * cos_phi * sin_obs + std::cos(coord.theta[j]) * cos_obs;
            for (size_t k = 0; k < t_size; ++k) {
                Real gamma_ = Gamma(i * interp.jet_3d, j, k);  // Get Gamma at the grid point.
                Real r = r_grid(i * interp.jet_3d, j, k);
                Real t_eng_ = coord.t[k];         // Get engine time at the grid point.
                Real beta = gammaTobeta(gamma_);  // Convert Gamma to beta.
                // Compute the Doppler factor: D = 1 / [Gamma * (1 - beta * cos_v)]
                doppler(i, j, k) = 1 / (gamma_ * (1 - beta * cos_v));
                // Compute the observed time: t_obs = [t_eng + (1 - cos_v) * r / c] * (1 + z)
                t_obs_grid(i, j, k) = (t_eng_ + (1 - cos_v) * r / con::c) * (1 + z);
            }
        }
    }
}

void Observer::updateRequired(MaskGrid& required, Array const& t_obs) {
    size_t t_obs_size = t_obs.size();

    // Loop over effective phi and theta grid points.
    for (size_t i = 0; i < eff_phi_size; i++) {
        size_t i_eff = i * interp.jet_3d;
        for (size_t j = 0; j < theta_size; j++) {
            // Skip observation times that are below the grid's start time
            size_t t_idx = 0;
            while (t_idx < t_obs_size && t_obs(t_idx) < t_obs_grid(i, j, 0)) {
                t_idx++;
            }

            // find the grid points that are required for the interpolation.
            for (size_t k = 0; k < t_size - 1 && t_idx < t_obs_size; k++) {
                Real const t_lo = t_obs_grid(i, j, k);
                Real const t_hi = t_obs_grid(i, j, k + 1);

                if (t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi) {
                    required(i_eff, j, k) = true;
                    required(i_eff, j, k + 1) = true;
                }

                for (; t_idx < t_obs_size && t_lo <= t_obs(t_idx) && t_obs(t_idx) < t_hi; t_idx++) {
                }
            }
        }
    }
}