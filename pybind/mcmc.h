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
#include "utilities.h"
#include "struct.h"

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
    using List = std::vector<double>;

    std::vector<LightCurveData> light_curve;
    std::vector<SpectrumData> spectrum;

    double estimate_chi2() const;
    void add_light_curve(double nu, List const& t, List const& Fv_obs, List const& Fv_err);
    void add_spectrum(double t, List const& nu, List const& Fv_obs, List const& Fv_err);
};

struct MultiBandModel {
    using List = std::vector<double>;
    using Grid = std::vector<std::vector<double>>;
    MultiBandModel() = delete;
    MultiBandModel(MultiBandData const& data);

    void configure(ConfigParams const& param);
    double estimate_chi2(Params const& param);
    Grid light_curves(Params const& param, List const& t, List const& nu);
    Grid spectra(Params const& param, List const& nu, List const& t);

   private:
    void build_system(Params const& param, Array const& t_eval, Observer& obs, SynElectronGrid& electrons,
                      SynPhotonGrid& photons);
    MultiBandData obs_data;
    ConfigParams config;
    // SynElectronGrid electrons;
    // SynPhotonGrid photons;
    // Observer obs;
    Array t_eval;
};