#pragma once

#include "../environment/medium.h"
#include "mesh.h"

inline Real structure_weight(Real Gamma) {
    return Gamma * std::sqrt(std::fabs((Gamma - 1) * Gamma));
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Builds an adaptive log-spaced grid with explicit break points.
 *
 * Places grid points at endpoints and break frequencies, then fills segments
 * between with approximately `pts_per_decade` log-spaced points. Nearby points
 * are merged automatically based on the local step size.
 *
 * @param lg2_min  Log2 of minimum grid value
 * @param lg2_max  Log2 of maximum grid value
 * @param breaks   Break frequencies (linear scale) to include
 * @param break_weights  Importance weights for each break (top values get refined regions)
 * @param pts_per_decade  Base density of fill points per decade
 * @param grid     Output array of grid point values (resized internally)
 * @param max_refined_breaks  Maximum number of breaks to refine (sorted by weight)
 * @param refine_radius_decades  Half-width of the refined region around each break (in decades)
 * @param refine_factor  Grid density multiplier inside refined regions
 * @return         Number of grid points
 * <!-- ************************************************************************************** -->
 */
size_t adaptive_grid_with_breaks(Real lg2_min, Real lg2_max, std::span<const Real> breaks,
                                 std::span<const Real> break_weights, Real pts_per_decade, Array& grid,
                                 size_t max_refined_breaks = 3, Real refine_radius_decades = 0.5,
                                 Real refine_factor = 3.0);

template <typename Ejecta>
Array find_jet_jumps(Ejecta const& jet, Real gamma_cut, [[maybe_unused]] bool is_axisymmetric) {
    // Detect sharp discontinuities in the Gamma0 profile via binary-search refinement.
    constexpr size_t n_scan = 512;
    constexpr Real eps = defaults::solver::binary_search_eps;
    const Real theta_lo = defaults::grid::theta_min;
    const Real theta_hi = con::pi / 2;
    const Real dtheta = (theta_hi - theta_lo) / (n_scan - 1);

    if (jet.Gamma0(0, theta_hi) >= gamma_cut) {
        return {theta_hi};
    }

    std::vector<Real> jumps;
    jumps.reserve(16);

    Real prev_th = theta_lo;
    Real prev_G = jet.Gamma0(0, theta_lo);

    for (size_t j = 1; j < n_scan; ++j) {
        const Real cur_th = theta_lo + dtheta * static_cast<Real>(j);
        const Real cur_G = jet.Gamma0(0, cur_th);

        if (prev_G >= gamma_cut || cur_G >= gamma_cut) {
            const Real dG = std::abs(cur_G - prev_G);
            const Real scale = std::max(prev_G - 1, cur_G - 1);

            if (scale > 0 && dG > 0.5 * scale) {
                Real lo = prev_th, hi = cur_th;
                while (hi - lo > eps) {
                    const Real mid = 0.5 * (lo + hi);
                    const Real G_mid = jet.Gamma0(0, mid);
                    if (std::abs(G_mid - prev_G) < std::abs(G_mid - cur_G)) {
                        lo = mid;
                    } else {
                        hi = mid;
                    }
                }
                jumps.push_back(prev_G > cur_G ? lo : hi);
            }
        }
        prev_th = cur_th;
        prev_G = cur_G;
    }

    return xt::adapt(jumps, {jumps.size()});
}

template <typename Ejecta>
Real find_theta_max(Ejecta const& jet, Real gamma_cut) {
    constexpr size_t n_scan = 512;
    const Real theta_lo = defaults::grid::theta_min;
    const Real theta_hi = con::pi / 2;

    const Real step = (theta_hi - theta_lo) / n_scan;
    for (Real th = theta_hi; th >= theta_lo; th -= step) {
        if (jet.Gamma0(0, th) >= gamma_cut) {
            return th;
        }
    }
    return theta_lo;
}

template <typename Ejecta, typename Medium>
Real jet_spreading_edge(Ejecta const& jet, Medium const& /*medium*/, Real phi, Real theta_min, Real theta_max,
                        Real /*t0*/) {
    const Real step = (theta_max - theta_min) / 256;
    Real theta_s = theta_min;
    Real dp_min = 0;

    for (Real theta = theta_min; theta <= theta_max; theta += step) {
        const Real th_lo = std::max(theta - step, theta_min);
        const Real th_hi = std::min(theta + step, theta_max);
        const Real dp = (jet.Gamma0(phi, th_hi) - jet.Gamma0(phi, th_lo)) / (th_hi - th_lo);

        if (dp < dp_min) {
            dp_min = dp;
            theta_s = theta;
        }
    }
    if (dp_min == 0) {
        theta_s = theta_max;
    }

    return theta_s;
}

template <typename Func>
Array inverse_CFD_sampling(Func&& pdf, Real min, Real max, size_t num,
                           size_t sample_num = defaults::sampling::theta_samples, bool log_sample = false) {
    using namespace boost::numeric::odeint;
    constexpr Real rtol = defaults::solver::ode_rtol;
    Array x_i = log_sample ? xt::eval(xt::logspace(std::log10(min), std::log10(max), sample_num))
                           : xt::eval(xt::linspace(min, max, sample_num));
    Array CDF_i = xt::zeros<Real>({sample_num});

    auto stepper = make_dense_output(rtol, rtol, runge_kutta_dopri5<Real>());
    stepper.initialize(0, min, (max - min) / 1e3);

    for (size_t k = 1, steps = 0; stepper.current_time() <= max;) {
        stepper.do_step(pdf);
        if (++steps > defaults::solver::max_ode_steps) {
            break;
        }
        while (k < x_i.size() && stepper.current_time() > x_i(k)) {
            stepper.calc_state(x_i(k), CDF_i(k));
            ++k;
        }
    }

    Array CDF_out = xt::linspace(CDF_i.front(), CDF_i.back(), num);
    Array x_out = xt::zeros<Real>({num});

    for (size_t k = 0; k < num; ++k) {
        for (size_t j = 0; j < sample_num; ++j) {
            if (CDF_out(k) <= CDF_i(j)) {
                if (j == 0) {
                    x_out(k) = x_i(j);
                } else {
                    Real slope = (x_i(j) - x_i(j - 1)) / (CDF_i(j) - CDF_i(j - 1));
                    x_out(k) = x_i(j - 1) + slope * (CDF_out(k) - CDF_i(j - 1));
                }
                break;
            }
        }
    }
    return x_out;
}

template <typename Ejecta>
Array adaptive_theta_grid(Ejecta const& jet, Real theta_min, Real theta_max, size_t base_pts, Real theta_v,
                          Real theta_resol, ThetaGridParams const& tgp = {}) {
    // Pre-scan: find peak structure weight, Gamma_peak, CDF estimate, and bright edge.
    // Bright edge = outermost angle where structure > 1% of peak; limits core beam extent
    // so structured jets (Gaussian, power-law) don't waste beam points in dim wings.
    constexpr size_t scan_pts = 100;
    const Real theta_extent = theta_max - theta_min;
    Real peak_weight = 0;
    Real Gamma_peak = 1.0;
    Real struct_sum = 0;
    size_t last_bright = 0;
    for (size_t i = 0; i <= scan_pts; ++i) {
        const Real Gamma = jet.Gamma0(0, theta_min + theta_extent * i / scan_pts);
        const Real w = structure_weight(Gamma);
        struct_sum += w;
        if (w > peak_weight) {
            peak_weight = w;
            Gamma_peak = Gamma;
            last_bright = i;
        } else if (w > 0.01 * peak_weight) {
            last_bright = i;
        }
    }
    const Real floor_weight = tgp.floor_fraction * peak_weight;
    const Real CDF_est = (struct_sum / scan_pts + floor_weight) * theta_extent;
    const Real theta_bright = theta_min + theta_extent * last_bright / scan_pts;

    const Real Gamma_v = jet.Gamma0(0, std::clamp(theta_v, theta_min, theta_max));
    Gamma_peak = std::max(Gamma_peak, Gamma_v);
    const Real doppler_alpha = tgp.doppler_alpha * std::sqrt(peak_weight / std::max(structure_weight(Gamma_v), 1.0));

    // Beaming density: b(θ) = Γ²θ/(1+Γ²θ²), derived from |d/dθ ln(δ³)|.
    // Core term (centered at θ=0): resolves jet core, extent limited to theta_bright.
    // View term (centered at θ_v): resolves observer's Doppler transition, only when θ_v
    // is outside the core beaming cone (θ_v·Γ > 3).
    constexpr Real beam_offset = 1.0;

    auto compute_beam_pts = [&](Real log_decades, Real beam_coeff) -> size_t {
        return static_cast<size_t>(std::max(0.0, log_decades - beam_offset) * theta_resol * beam_coeff);
    };

    const Real Gamma_peak_sq = Gamma_peak * Gamma_peak;
    const Real Gamma_v_sq = Gamma_v * Gamma_v;
    const size_t core_beam_pts =
        compute_beam_pts(std::log10(std::max(1.0, Gamma_peak * (theta_bright - theta_min))), tgp.core_beam_coeff);
    const size_t view_beam_pts =
        (theta_v * Gamma_peak > 3.0)
            ? compute_beam_pts(std::log10(std::max(1.0, Gamma_v * std::max(theta_v - theta_min, theta_max - theta_v))),
                               tgp.view_beam_coeff)
            : 0;
    const size_t total_pts = base_pts + core_beam_pts + view_beam_pts;

    // Calibrate beam weights so each term attracts its target fraction of points.
    auto calibrate_beam_weight = [&](size_t n_pts, Real beam_cdf) -> Real {
        return (n_pts > 0 && beam_cdf > 0) ? static_cast<Real>(n_pts) / base_pts * CDF_est / beam_cdf : 0.0;
    };

    const Real core_weight =
        calibrate_beam_weight(core_beam_pts, 0.5 * std::log((1.0 + Gamma_peak_sq * theta_max * theta_max) /
                                                            (1.0 + Gamma_peak_sq * theta_min * theta_min)));
    const Real theta_v_left = theta_v - theta_min;
    const Real theta_v_right = theta_max - theta_v;
    const Real view_weight =
        calibrate_beam_weight(view_beam_pts, 0.5 * (std::log(1.0 + Gamma_v_sq * theta_v_left * theta_v_left) +
                                                    std::log(1.0 + Gamma_v_sq * theta_v_right * theta_v_right)));

    auto pdf_eqn = [=, &jet](Real const& /*cdf*/, Real& pdf, Real theta) {
        const Real Gamma = jet.Gamma0(0, theta);
        const Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
        const Real doppler = (1 - beta) / (1 - beta * std::cos(theta - theta_v));
        const Real structure = structure_weight(Gamma);
        const Real dtheta = theta - theta_v;
        pdf = core_weight * Gamma_peak_sq * theta / (1.0 + Gamma_peak_sq * theta * theta) +
              view_weight * Gamma_v_sq * std::fabs(dtheta) / (1.0 + Gamma_v_sq * dtheta * dtheta) +
              (1 + doppler_alpha * doppler) * structure + floor_weight;
    };

    return inverse_CFD_sampling(pdf_eqn, theta_min, theta_max, total_pts, defaults::sampling::theta_samples,
                                /*log_sample=*/true);
}

Array jump_refinement_grid(Array const& jumps, Real theta_min, Real theta_max, Real avg_spacing);

template <typename Ejecta>
Array adaptive_phi_grid(Ejecta const& jet, size_t phi_num, Real theta_v, Array const& theta_grid,
                        bool is_axisymmetric) {
    if (theta_v == 0 && is_axisymmetric) {
        return xt::linspace(0., 2 * con::pi, phi_num);
    }

    const Real cos_tv = std::cos(theta_v);
    const Real sin_tv = std::sin(theta_v);
    const size_t n_theta = theta_grid.size();

    // Precompute dcos_theta bin widths from theta grid
    std::vector<Real> dcos(n_theta);
    for (size_t it = 0; it < n_theta; ++it) {
        const Real left = (it == 0) ? 0.0 : 0.5 * (theta_grid(it - 1) + theta_grid(it));
        const Real right = (it == n_theta - 1) ? theta_grid(it) : 0.5 * (theta_grid(it) + theta_grid(it + 1));
        dcos[it] = std::fabs(std::cos(left) - std::cos(right));
    }

    // Weight function: Doppler-weighted structure summed across theta bins
    auto phi_weight = [=, &jet](Real phi) {
        const Real cos_phi = std::cos(phi);
        Real w = 0;
        for (size_t it = 0; it < n_theta; ++it) {
            const Real theta = theta_grid(it);
            const Real Gamma = jet.Gamma0(phi, theta);
            const Real beta = std::sqrt(std::fabs(Gamma * Gamma - 1)) / Gamma;
            const Real cos_alpha = std::cos(theta) * cos_tv + std::sin(theta) * sin_tv * cos_phi;
            const Real a = (1 - beta) / (1 - beta * cos_alpha);
            w += a * structure_weight(Gamma) * dcos[it];
        }
        return w;
    };

    // Pre-scan to find peak weight for floor calculation
    constexpr size_t scan_pts = 100;
    Real peak_weight = 0;
    for (size_t s = 0; s <= scan_pts; ++s) {
        const Real phi = 2 * con::pi * static_cast<Real>(s) / scan_pts;
        peak_weight = std::max(peak_weight, phi_weight(phi));
    }
    const Real floor_weight = 0.05 * peak_weight;

    auto eqn = [=](Real const& /*cdf*/, Real& pdf, Real phi) { pdf = phi_weight(phi) + floor_weight; };

    return inverse_CFD_sampling(eqn, 0, 2 * con::pi, phi_num);
}

template <typename Arr>
Array merge_grids(Arr const& arr1, Arr const& arr2) {
    std::vector<Real> result;
    result.reserve(arr1.size() + arr2.size());

    size_t i = 0;
    size_t j = 0;
    auto add_unique = [&](Real val) {
        if (result.empty() || result.back() != val) {
            result.push_back(val);
        }
    };

    while (i < arr1.size() && j < arr2.size()) {
        if (arr1[i] <= arr2[j]) {
            add_unique(arr1[i++]);
            if (arr1[i - 1] == arr2[j]) {
                j++;
            }
        } else {
            add_unique(arr2[j++]);
        }
    }
    while (i < arr1.size()) {
        add_unique(arr1[i++]);
    }
    while (j < arr2.size()) {
        add_unique(arr2[j++]);
    }

    return xt::adapt(result);
}

/**
 * <!-- ************************************************************************************** -->
 * @brief Estimates deceleration time for a given angle from ejecta and medium properties.
 * @details Integrates swept mass outward until m_swept = m_jet/Gamma0(phi,theta). Returns lab-frame time.
 * <!-- ************************************************************************************** -->
 */
template <typename Ejecta, typename Medium>
Real estimate_t_dec(Ejecta const& jet, Medium const& medium, Real phi, Real theta) {
    /*const Real gamma = jet.Gamma0(phi, theta);
    if (gamma <= 1) {
        return con::inf;
    }
    const Real beta = physics::relativistic::gamma_to_beta(gamma);
    Real m_jet = jet.eps_k(phi, theta) / (gamma * con::c2);
    if constexpr (HasSigma<Ejecta>) {
        m_jet /= (1.0 + jet.sigma0(phi, theta));
    }
    const Real target = m_jet / gamma;

    constexpr size_t N = 256;
    const Real u_min = std::log(1e-3);
    const Real u_max = u_min + 40 * std::log(10.0);
    const Real du = (u_max - u_min) / N;

    Real mass = 0;
    Real r_prev = std::exp(u_min);
    Real f_prev = medium.rho(phi, theta, r_prev) * r_prev * r_prev;

    for (size_t i = 1; i <= N; ++i) {
        const Real r_i = std::exp(u_min + i * du);
        const Real f_i = medium.rho(phi, theta, r_i) * r_i * r_i;
        const Real dr = r_i - r_prev;
        mass += 0.5 * (f_prev + f_i) * dr;

        if (mass >= target) {
            const Real r_dec = r_prev + (target - (mass - 0.5 * (f_prev + f_i) * dr)) / f_i;
            return r_dec * (1 - beta) / (beta * con::c);
        }
        f_prev = f_i;
        r_prev = r_i;
    }
    return std::exp(u_max) * (1 - beta) / (beta * con::c);*/

    const Real gamma = jet.Gamma0(phi, theta);
    const Real beta = physics::relativistic::gamma_to_beta(gamma);
    Real m_jet = jet.eps_k(phi, theta) / (gamma * con::c2);
    if constexpr (HasSigma<decltype(jet)>) {
        m_jet /= (1.0 + jet.sigma0(phi, theta));
    }
    const Real target = m_jet / gamma;

    constexpr Real r_min = 1e-3;
    const Real r_max = r_min * std::pow(10.0, 40.0);
    if (target <= 0) {
        return r_min * (1 - beta) / (beta * con::c);
    }

    if constexpr (std::is_same_v<Medium, ISM>) {
        const Real rho = medium.rho(phi, theta, r_min);
        if (rho > 0) {
            const Real r3_dec = r_min * r_min * r_min + 3 * target / rho;
            const Real r_dec = std::cbrt(std::max(r3_dec, 0.0));
            return std::min(r_dec, r_max) * (1 - beta) / (beta * con::c);
        }
        return r_max * (1 - beta) / (beta * con::c);
    }

    auto rho = [&](Real r) { return medium.rho(phi, theta, r); };

    // Trapezoidal integration in log-space with early exit at deceleration radius
    constexpr size_t N = 256;
    const Real u_min = std::log(1e-3);
    const Real u_max = u_min + 40 * std::log(10.0);
    const Real du = (u_max - u_min) / N;

    Real mass = 0;
    Real r_prev = std::exp(u_min);
    Real f_prev = rho(r_prev) * r_prev * r_prev;

    for (size_t i = 1; i <= N; ++i) {
        const Real r_i = std::exp(u_min + i * du);
        const Real f_i = rho(r_i) * r_i * r_i;
        const Real dr = r_i - r_prev;
        mass += 0.5 * (f_prev + f_i) * dr;

        if (mass >= target) {
            const Real r_dec = r_prev + (target - (mass - 0.5 * (f_prev + f_i) * dr)) / f_i;
            return r_dec * (1 - beta) / (beta * con::c);
        }
        f_prev = f_i;
        r_prev = r_i;
    }
    return std::exp(u_max) * (1 - beta) / (beta * con::c);
}

/// Create a two-segment logspace grid with extra resolution before t_cross.
/// Post-crossing density matches base logspace (base_t_num); extra points all go to pre-crossing.
Array shock_crossing_refinement_logspace(Real t_start, Real t_end, Real t_cross, size_t t_num, size_t base_t_num);

template <typename Ejecta, typename Medium>
void Coord::detect_symmetry(Ejecta const& jet, Medium const& medium, Real t_min, Real t_max, Real z, Real t_resol) {
    const size_t theta_size = theta.size();

    if (jet.spreading || !medium.isotropic) {
        symmetry = Symmetry::structured;
        theta_reps.resize(theta_size);
        std::iota(theta_reps.begin(), theta_reps.end(), size_t(0));
        return;
    }

    const Real phi0 = phi(0);
    theta_reps.clear();
    theta_reps.reserve(theta_size);
    theta_reps.push_back(0);

    // Logspaced engine time for time-dependent injection checks
    const size_t t_check = std::max<size_t>(static_cast<size_t>(std::log10(t_max / t_min) * t_resol), 24);
    const Array temp_t = xt::logspace(std::log10(t_min / (1 + z)), std::log10(t_max / (1 + z)), t_check);

    auto jet_ic_differs = [&](size_t ja, size_t jb) {
        const Real theta_a = theta(ja);
        const Real theta_b = theta(jb);
        if (jet.eps_k(phi0, theta_a) != jet.eps_k(phi0, theta_b)) {
            return true;
        }
        if (jet.Gamma0(phi0, theta_a) != jet.Gamma0(phi0, theta_b)) {
            return true;
        }
        if constexpr (HasSigma<Ejecta>) {
            if (jet.sigma0(phi0, theta_a) != jet.sigma0(phi0, theta_b)) {
                return true;
            }
        }
        if constexpr (HasDedt<Ejecta>) {
            for (size_t k = 0; k < t_check; ++k) {
                if (jet.deps_dt(phi0, theta_a, temp_t(k)) != jet.deps_dt(phi0, theta_b, temp_t(k))) {
                    return true;
                }
            }
        }
        if constexpr (HasDmdt<Ejecta>) {
            for (size_t k = 0; k < t_check; ++k) {
                if (jet.dm_dt(phi0, theta_a, temp_t(k)) != jet.dm_dt(phi0, theta_b, temp_t(k))) {
                    return true;
                }
            }
        }
        return false;
    };

    for (size_t j = 1; j < theta_size; ++j) {
        if (jet_ic_differs(j - 1, j)) {
            theta_reps.push_back(j);
        }
    }

    if (theta_reps.size() == 1) {
        symmetry = Symmetry::isotropic;
    } else if (theta_reps.size() < theta_size) {
        symmetry = Symmetry::piecewise;
    } else {
        symmetry = Symmetry::phi_symmetric;
    }
}

/// Result of scanning cell time bounds for grid sizing.
struct TimeScanResult {
    Real min_t_start_raw;
    Real min_t_start;
    bool has_early_point;
    Real max_t_refine;
    xt::xtensor<Real, 2> t_dec; // cached deceleration times [phi_size x theta_size]
};

/// Scan all cells to determine global time bounds for grid construction.
template <typename Ejecta, typename Medium>
TimeScanResult scan_time_bounds(Coord const& coord, Ejecta const& jet, Medium const& medium, Real t_min, Real t_end,
                                Real z, bool is_rvs, size_t phi_size) {
    const size_t theta_size = coord.theta.size();
    const Real cos_tv = std::cos(coord.theta_view);
    const Real sin_tv = std::sin(coord.theta_view);

    Real min_raw = t_end, min_guarded = t_end, min_cut = t_end, max_ref = 0;
    auto t_dec_cache = xt::xtensor<Real, 2>::from_shape({phi_size, theta_size});

    for (size_t i = 0; i < phi_size; ++i) {
        for (size_t j = 0; j < theta_size; ++j) {
            const Real b = physics::relativistic::gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
            const Real cos_a =
                std::cos(coord.theta(j)) * cos_tv + std::sin(coord.theta(j)) * sin_tv * std::cos(coord.phi(i));
            const Real ts = 0.99 * t_min * (1 - b) / (1 - cos_a * b) / (1 + z);

            const Real td = estimate_t_dec(jet, medium, coord.phi(i), coord.theta(j));
            t_dec_cache(i, j) = td;
            Real cut = std::min(0.01 * td, 1e-2 * unit::sec);
            if (is_rvs) {
                cut = std::min(cut, 0.01 * jet.T0);
                max_ref = std::max(max_ref, 10.0 * std::max(td, jet.T0));
            }

            min_raw = std::min(min_raw, ts);
            min_guarded = std::min(min_guarded, std::max(ts, cut));
            min_cut = std::min(min_cut, cut);
        }
    }

    return {min_raw, min_guarded, min_raw < min_cut, max_ref, std::move(t_dec_cache)};
}

template <typename Ejecta, typename Medium>
void build_time_grid(Coord& coord, Ejecta const& jet, Medium const& medium, Real t_min, Real t_max, Real z,
                     Real t_resol, bool is_axisymmetric, bool is_rvs) {
    const size_t theta_size = coord.theta.size();
    const size_t phi_size = is_axisymmetric ? 1 : coord.phi.size();
    const Real t_end = 1.01 * t_max / (1 + z);

    auto [min_t_start_raw, min_t_start, has_early_point, max_t_refine, t_dec] =
        scan_time_bounds(coord, jet, medium, t_min, t_end, z, is_rvs, phi_size);

    size_t n_pre_extra = 0;
    if (is_rvs && max_t_refine > min_t_start) {
        n_pre_extra = static_cast<size_t>(std::log10(std::min(max_t_refine, t_end) / min_t_start) * t_resol);
    }
    const size_t t_num_base_fwd = static_cast<size_t>(std::max(std::log10(t_end / min_t_start), 1.0) * t_resol);
    const size_t t_num_base = t_num_base_fwd + n_pre_extra;
    const size_t t_num = t_num_base + (has_early_point ? 1 : 0);

    coord.t = xt::zeros<Real>({phi_size, theta_size, t_num});

    auto make_grid = [&](Real t_start, Real t_dec) -> Array {
        if (is_rvs) {
            return shock_crossing_refinement_logspace(t_start, t_end, 10.0 * std::max(t_dec, jet.T0), t_num_base,
                                                      t_num_base_fwd);
        }
        return xt::logspace(std::log10(t_start), std::log10(t_end), t_num_base);
    };

    auto store_time_grid = [&](size_t i, size_t j, Real early_t, auto const& grid) {
        if (has_early_point) {
            xt::view(coord.t, i, j, 0) = early_t;
            xt::view(coord.t, i, j, xt::range(1, t_num)) = grid;
        } else {
            xt::view(coord.t, i, j, xt::all()) = grid;
        }
    };

    if (coord.symmetry >= Symmetry::phi_symmetric) {
        for (size_t r = 0; r < coord.theta_reps.size(); ++r) {
            const size_t j_rep = coord.theta_reps[r];
            const size_t j_end = (r + 1 < coord.theta_reps.size()) ? coord.theta_reps[r + 1] : theta_size;
            auto grid = make_grid(min_t_start, t_dec(0, j_rep));
            for (size_t i = 0; i < phi_size; ++i) {
                for (size_t j = j_rep; j < j_end; ++j) {
                    store_time_grid(i, j, min_t_start_raw, grid);
                }
            }
        }
    } else {
        const Real cos_tv = std::cos(coord.theta_view);
        const Real sin_tv = std::sin(coord.theta_view);
        for (size_t i = 0; i < phi_size; ++i) {
            for (size_t j = 0; j < theta_size; ++j) {
                const Real b = physics::relativistic::gamma_to_beta(jet.Gamma0(coord.phi(i), coord.theta(j)));
                const Real cos_a =
                    std::cos(coord.theta(j)) * cos_tv + std::sin(coord.theta(j)) * sin_tv * std::cos(coord.phi(i));
                const Real t_raw = 0.99 * t_min * (1 - b) / (1 - cos_a * b) / (1 + z);

                const Real td = t_dec(i, j);
                Real cut = std::min(0.01 * td, 1e-2 * unit::sec);
                if (is_rvs) {
                    cut = std::min(cut, 0.01 * jet.T0);
                }
                const Real t_start = std::max(t_raw, cut);
                auto grid = make_grid(t_start, td);
                store_time_grid(i, j, 0.99 * std::min(t_raw, cut), grid);
            }
        }
    }
}

template <typename Ejecta, typename Medium>
Coord auto_grid(Ejecta const& jet, Medium const& medium, Array const& t_obs, Real theta_cut, Real theta_view, Real z,
                bool is_rvs, Real phi_resol, Real theta_resol, Real t_resol, bool is_axisymmetric) {
    Coord coord;
    coord.theta_view = theta_view;

    size_t min_theta_num = defaults::grid::min_theta_points;

    const Array jet_jumps = find_jet_jumps(jet, con::Gamma_cut, is_axisymmetric);
    Real jet_edge = find_theta_max(jet, con::Gamma_cut);
    for (size_t i = 0; i < jet_jumps.size(); ++i) {
        jet_edge = std::max(jet_edge, jet_jumps(i));
    }
    const Real theta_min = defaults::grid::theta_min;
    const Real theta_max = std::min(jet_edge, theta_cut);

    const size_t theta_num = min_theta_num + static_cast<size_t>((theta_max - theta_min) * 180 / con::pi * theta_resol);

    constexpr ThetaGridParams tgp{};
    const Array base_theta = adaptive_theta_grid(jet, theta_min, theta_max, theta_num, theta_view, theta_resol, tgp);

    const Real avg_spacing = (theta_max - theta_min) / base_theta.size();
    const Array feature_theta = jump_refinement_grid(jet_jumps, theta_min, theta_max, avg_spacing);
    coord.theta = merge_grids(base_theta, feature_theta);

    const size_t phi_base = std::max<size_t>(static_cast<size_t>(360 * phi_resol), 1);
    const Real doppler_sharpness = jet.Gamma0(0, theta_view) * std::sin(theta_view);
    const Real phi_boost = std::sqrt(std::max(doppler_sharpness / (2 * con::pi), 1.0));
    const size_t phi_num = std::clamp(static_cast<size_t>(phi_base * phi_boost), size_t(1), phi_base * 5);

    if (phi_num <= 2) {
        coord.phi = xt::linspace(0., 2 * con::pi, phi_num);
    } else {
        coord.phi = adaptive_phi_grid(jet, phi_num, theta_view, coord.theta, is_axisymmetric);
    }
    // Shift phi grid by half the first spacing so no point lands on the
    // Doppler peak at phi=0
    if (phi_num >= 2) {
        const Real shift = 0.5 * (coord.phi(1) - coord.phi(0));
        coord.phi += shift;
    }

    if (t_obs.size() == 0) {
        assert(false && "auto_grid: t_obs is empty");
        return coord;
    }

    const Real t_max = *std::ranges::max_element(t_obs);
    const Real t_min = *std::ranges::min_element(t_obs);

    coord.detect_symmetry(jet, medium, t_min, t_max, z, t_resol);
    build_time_grid(coord, jet, medium, t_min, t_max, z, t_resol, is_axisymmetric, is_rvs);

    return coord;
}
