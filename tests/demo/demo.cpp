
#include "afterglow.h"
auto test_reverse_shock(double xi, double sigma, bool output = true) {
    Real E_iso = 1e49 * unit::erg;
    Real theta_c = 0.1;
    Real theta_v = 0;

    Real n_ism = 1 / unit::cm3;
    Real Gamma0 = 200;
    Real z = 0;

    RadParams rad_fwd;
    rad_fwd.eps_e = 1e-2;
    rad_fwd.eps_B = 1e-4;

    RadParams rad_rvs = rad_fwd;

    Array t_obs = xt::logspace(std::log10(1e-9 * unit::sec), std::log10(1e9 * unit::sec), 230);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.sigma0 = math::tophat(theta_c, sigma);

    jet.T0 = calc_engine_duration(E_iso, n_ism, Gamma0, xi);

    std::cout << "T0: " << jet.T0 / unit::sec << ' ' << xi << ' ' << sigma << ' '
              << thin_shell_dec_radius(E_iso, n_ism, Gamma0) / (2 * Gamma0 * Gamma0) / unit::sec << std::endl;

    Coord coord = auto_grid(jet, t_obs, 0.6, theta_v, z, 0.5, 10, 50);

    auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad_fwd, rad_rvs);
    // auto f_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);
    // auto r_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    if (output) {
        /* write_npz("rshock-data/coord" + std::to_string(xi) + "-" + std::to_string(sigma), coord);
        write_npz("rshock-data/f_shock" + std::to_string(xi) + "-" + std::to_string(sigma), f_shock);
        write_npz("rshock-data/r_shock" + std::to_string(xi) + "-" + std::to_string(sigma), r_shock);
        */
    }

    auto t_cross = coord.t(0, 0, r_shock.injection_idx(0, 0)) / unit::sec;
    auto duration = jet.T0 / unit::sec;
    auto t_dec = sedov_length(E_iso, n_ism) / std::pow(Gamma0, 8. / 3) / 2 / con::c / unit::sec;

    return std::make_tuple(t_cross, duration, t_dec);
}

void test_spreading() {
    Real E_iso = 1e51 * unit::erg;
    Real theta_c = 10 * unit::deg;
    Real theta_v = 0;

    Real n_ism = 1 / unit::cm3;

    Real Gamma0 = 300;
    Real z = 0;

    RadParams rad_fwd;
    rad_fwd.eps_e = 1e-2;
    rad_fwd.eps_B = 1e-3;

    Array t_obs = xt::logspace(std::log10(0.1 * unit::sec), std::log10(1e8 * unit::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);

    jet.spreading = true;

    Coord coord = auto_grid(jet, t_obs, con::pi / 2, theta_v, z);

    auto shock = generate_fwd_shock(coord, medium, jet, rad_fwd);

    //write_npz("spreading-data/shock", shock);
}

void test_ic(Real theta_c_) {
    Real E_iso = 1e51 * unit::erg;
    Real theta_c = theta_c_ * unit::deg;
    Real theta_v = 0;

    Real n_ism = 1 / unit::cm3;

    Real Gamma0 = 300;
    Real z = 0;

    RadParams rad_fwd;
    rad_fwd.eps_e = 1e-2;
    rad_fwd.eps_B = 1e-3;

    Array t_obs = xt::logspace(std::log10(0.1 * unit::sec), std::log10(1e8 * unit::sec), 5);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);

    Coord coord = auto_grid(jet, t_obs, con::pi / 2, theta_v, z);

    auto shock = generate_fwd_shock(coord, medium, jet, rad_fwd);

    auto elec = generate_syn_electrons(shock);

    auto photons = generate_syn_photons(shock, elec);

    auto ic = generate_IC_photons(elec, photons, false);

    Observer obs;

    Real lumi_dist = 1;

    obs.observe(coord, shock, lumi_dist, z);

    auto flux = obs.specific_flux(t_obs, 1e22 * unit::Hz, ic);

    // write_npz("spreading-data/shock", shock);
}

void test_FRS() {
    Real E_iso = std::pow(10, 54.43) * unit::erg;
    Real theta_c = 0.1;
    Real theta_v = 0;

    Real n_ism = std::pow(10, -0.32) / unit::cm3;

    RadParams rad_fwd;
    rad_fwd.eps_e = std::pow(10, -1.81);
    rad_fwd.eps_B = std::pow(10, -3.43);
    rad_fwd.p = 2.77;

    RadParams rad_rvs;
    rad_rvs.eps_e = std::pow(10, -0.14);
    rad_rvs.eps_B = std::pow(10, -5.10);
    rad_rvs.p = 2.71;

    Real Gamma0 = std::pow(10, 2.14);
    Real z = 1.88;

    Array t_obs = xt::logspace(std::log10(1e2 * unit::sec), std::log10(1e8 * unit::sec), 130);

    ISM medium(n_ism);

    Ejecta jet;

    jet.eps_k = math::tophat(theta_c, E_iso / (4 * con::pi));
    jet.Gamma0 = math::tophat(theta_c, Gamma0);
    jet.T0 = 1 * unit::sec;

    std::cout << "FRS:" << sedov_length(E_iso, n_ism) / (2 * con::c * std::pow(Gamma0, 8. / 3)) / unit::sec << ' '
              << shell_thickness_param(E_iso, n_ism, Gamma0, jet.T0) << std::endl;

    Coord coord = auto_grid(jet, t_obs, con::pi / 2, theta_v, z, 0.5, 100, 20);
    auto [f_shock, r_shock] = generate_shock_pair(coord, medium, jet, rad_fwd, rad_rvs);
    // auto f_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);
    // auto f_shock = generate_fwd_shock(coord, medium, jet, eps_e_rs, eps_B_rs);
    // auto r_shock = generate_fwd_shock(coord, medium, jet, eps_e, eps_B);

    auto elec = generate_syn_electrons(f_shock);
    auto elec_rs = generate_syn_electrons(r_shock);
    auto photons = generate_syn_photons(f_shock, elec);
    auto photons_rs = generate_syn_photons(r_shock, elec_rs);

    Observer obs;

    Real lumi_dist = 1;

    obs.observe(coord, f_shock, lumi_dist, z);

    auto flux = obs.specific_flux(t_obs, 1e17 * unit::Hz, photons);

    obs.observe(coord, r_shock, lumi_dist, z);

    auto flux_rs = obs.specific_flux(t_obs, 1e17 * unit::Hz, photons_rs);

    /*write_npz("frs/flux", flux);
    write_npz("frs/flux_rs", flux_rs);
    write_npz("frs/t_obs", t_obs);

    write_npz("frs/coord", coord);
    write_npz("frs/f_shock", f_shock);
    write_npz("frs/r_shock", r_shock);
    write_npz("frs/elec", elec);
    write_npz("frs/elec_rs", elec_rs);
    write_npz("frs/photons", photons);
    write_npz("frs/photons_rs", photons_rs);*/
}

int main() {
    for (double theta_c = 5; theta_c <= 30; theta_c += 0.01) {
        test_ic(theta_c);
    }
    return 0;
    /*for (Real i = 1; i < 30;) {
        test_ic(i);
        i += 0.01;
    }

    return 0;
    test_FRS();*/

    double xi[] = {0.001, 0.01, 0.1, 1, 2, 3, 5, 10, 100};
    // double xi[] = {100};
    double sigma[] = {0, 0.01, 0.05, 0.1, 1, 10, 100};

    // double sigma[] = {0};

    for (auto x : xi) {
        for (auto s : sigma) {
            test_reverse_shock(x, s);
        }
    }

    double xi2[] = {0.001, 0.01, 0.1, 1, 5, 10, 100, 1000, 10000};
    Array sigma2 = xt::logspace(std::log10(1e-5), std::log10(100), 100);

    //double sigma2[] = {0, 0, 0, 0, 0, 0};

    for (auto x : xi2) {
        std::ofstream out("rshock-data/crossing-time-" + std::to_string(x) + ".txt");
        for (auto s : sigma2) {
            auto [t_cross, duration, t_dec] = test_reverse_shock(x, s, false);
            out << x << ' ' << s << ' ' << t_cross << ' ' << duration << ' ' << t_dec << std::endl;
        }
    }

    return 0;
}
