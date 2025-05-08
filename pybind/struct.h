#pragma once

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "utilities.h"

struct Params {
    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double theta_v{0};
    double theta_w{con::pi / 2};
    double p{2.3};
    double p_r{p};
    double eps_e{0.1};
    double eps_B{0.01};
    double eps_e_r{eps_e};
    double eps_B_r{eps_B};
    double n_ism{1};
    double A_star{0.01};
    double xi{1};
    double xi_r{xi};
    double k_jet{2};
};

struct ConfigParams {
    double lumi_dist{1e26};
    double z{0};
    std::string medium{"ism"};
    std::string jet{"tophat"};
    std::function<double(double, double)> eps_k_profile;
    std::function<double(double, double)> Gamma0_profile;
    std::function<double(double, double, double)> rho_profile;
    std::function<double(double, double, double)> m_profile; 
    std::optional<double> T0{std::nullopt};
    size_t t_grid{24};
    size_t phi_grid{24};
    size_t theta_grid{24};
    double rtol{1e-5};
    bool spreading{false};
};