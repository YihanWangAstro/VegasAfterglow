#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <filesystem>
#include <format>
#include <fstream>

#include "afterglow.h"

void tests(Real phi_resol, Real theta_resol, Real t_resol, Real n_ism, Real eps_e, Real eps_B, Real p, Real E_iso,
           Real Gamma0, Real theta_c, Real theta_v, bool verbose = false) {
    Real z = 0.009;

    Real lumi_dist = 1.23e26 * unit::cm;

    size_t t_num = 100;

    if (verbose) {
        t_num = 100;
    }

    Array t_obs = xt::logspace(std::log10(1e3 * unit::sec), std::log10(1e7 * unit::sec), t_num);

    ISM medium(n_ism);

    GaussianJet jet(theta_c, E_iso, Gamma0);

    Coord coord = auto_grid(jet, t_obs, 0.6, theta_v, z, phi_resol, theta_resol, t_resol);

    RadParams rad_fwd;
    rad_fwd.eps_e = eps_e;
    rad_fwd.eps_B = eps_B;
    rad_fwd.p = p;

    RadParams rad_rev = rad_fwd;

    Shock f_shock = generate_fwd_shock(coord, medium, jet, rad_fwd);
    // auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad_fwd, rad_rev);

    Observer obs;

    //obs.observe_at(t_obs, coord, f_shock, lumi_dist, z);

    obs.observe(coord, f_shock, lumi_dist, z);

    auto syn_e = generate_syn_electrons(f_shock);

    auto syn_ph = generate_syn_photons(f_shock, syn_e);

    // auto ic = generate_IC_photons(syn_e, syn_ph);

    Real nu_obs = eVtoHz(1 * unit::keV);

    Array F_nu = obs.specific_flux(t_obs, nu_obs, syn_ph);

    if (verbose) {
        //write_npz(std::format("F_nu{:.1f}-{:.1f}-{:.1f}", phi_resol, theta_resol, t_resol), "F_nu",
        //          xt::eval(F_nu / unit::Jy), "t_obs", xt::eval(t_obs / unit::sec));
    }

    return;
}

int main() {
    Real n_ism = 2 / unit::cm3;
    Real eps_e = 1e-2;
    Real eps_B = 1e-4;
    Real p = 2.1;
    Real Gamma0 = 300;

    Array E_iso = xt::logspace(std::log10(1e48 * unit::erg), std::log10(1e52 * unit::erg), 100);
    Array theta_v = xt::linspace(0.01, 0.5, 10);

    Real phi_resolu[] = {0.1, 0.2, 0.4, 0.8, 1.6};
    Real theta_resolu[] = {1, 2, 4, 8, 16};
    Real t_resolu[] = {5, 10, 20, 40, 80};

    for (size_t i = 0; i < 5; ++i) {
        tests(phi_resolu[i], theta_resolu[i], t_resolu[i], n_ism, eps_e, eps_B, p, 1e52 * unit::erg, Gamma0, 0.1, 0.3,
              true);
    }

    Real benchmark_phi_resolu[] = {0.1, 0.2, 0.4};
    Real benchmark_theta_resolu[] = {1, 2, 4};
    Real benchmark_t_resolu[] = {5, 10, 20};

    for (size_t l = 0; l < 3; ++l) {
        std::ofstream file(std::format("benchmark{:.1f}-{:.1f}-{:.1f}.txt", benchmark_phi_resolu[l],
                                       benchmark_theta_resolu[l], benchmark_t_resolu[l]));

        for (size_t i = 0; i < E_iso.size(); ++i) {
            auto start = std::chrono::high_resolution_clock::now();
            for (size_t k = 0; k < theta_v.size(); ++k) {
                tests(benchmark_phi_resolu[l], benchmark_theta_resolu[l], benchmark_t_resolu[l], n_ism, eps_e, eps_B, p,
                      E_iso[i], Gamma0, 0.1, theta_v[k]);
                // tests(r, r, r, n_ism, eps_e, eps_B, p, E_iso[i], Gamma0, theta_c[j], 0);
            }
            auto end = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
            file << duration.count() / 1000000. / theta_v.size() << std::endl;
        }
    }
    return 0;
}
