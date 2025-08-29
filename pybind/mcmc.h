//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#pragma once

#include <iostream>
#include <vector>

#include "afterglow.h"
#include "macros.h"
#include "mesh.h"
#include "pybind.h"
#include "utilities.h"
struct LightCurveData {
    double nu{0};
    Array t;
    Array Fv_obs;
    Array Fv_err;
    Array Fv_model;

    double estimate_chi2() const;
};

struct SpectrumData {
    double t{0};
    Array nu;
    Array Fv_obs;
    Array Fv_err;
    Array Fv_model;

    double estimate_chi2() const;
};

struct MultiBandData {
    std::vector<LightCurveData> light_curve;
    std::vector<SpectrumData> spectrum;

    double estimate_chi2() const;

    void add_light_curve(double nu, PyArray const& t, PyArray const& Fv_obs, PyArray const& Fv_err);
    void add_spectrum(double t, PyArray const& nu, PyArray const& Fv_obs, PyArray const& Fv_err);
};

struct Params {
    double E_iso{1e52};
    double Gamma0{300};
    double theta_c{0.1};
    double theta_cG{0.1};
    double theta_v{0};
    double theta_w{con::pi / 2};
    double p{2.3};
    double eps_e{0.1};
    double eps_B{0.01};
    double n_ism{1};
    double A_star{0.01};
    double xi_e{1};
    double k_jet{2};
    double kG_jet{2};
};

struct ConfigParams {
    double lumi_dist{1e26};
    double z{0};
    std::string medium{"ism"};
    std::string jet{"tophat"};
    Real phi_resol{0.25};
    Real theta_resol{1};
    Real t_resol{3};
    double rtol{1e-5};
    bool reverse_shock{false};
    bool forward_SSC{false};
    bool reverse_SSC{false};
};

struct MultiBandModel {
    MultiBandModel() = delete;
    MultiBandModel(MultiBandData const& data);

    void configure(ConfigParams const& param);
    double estimate_chi2(Params const& param);
    PyGrid light_curves(Params const& param, PyArray const& t, PyArray const& nu);
    PyGrid spectra(Params const& param, PyArray const& nu, PyArray const& t);

   private:
    void build_system(Params const& param, Array const& t_eval, Observer& obs, SynElectronGrid& electrons,
                      SynPhotonGrid& photons);
    MultiBandData obs_data;
    ConfigParams config;
    Array t_eval;
};