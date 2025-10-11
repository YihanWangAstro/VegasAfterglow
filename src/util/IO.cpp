//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include "IO.h"

#include <fstream>

#include "macros.h"

void write_csv(std::string const& filename, Array const& array, Real unit) {
    std::ofstream file(filename + ".csv");
    for (size_t i = 0; i < array.size(); ++i) {
        file << array(i) / unit << "\n";
    }
}

void write_csv(std::string const& filename, MeshGrid const& grid, Real unit) {
    std::ofstream file(filename + ".csv");
    const auto shape = grid.shape();
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            file << grid(i, j) / unit;
            if (j + 1 != shape[1])
                file << ",";
        }
        file << "\n";
    }
}

void write_csv(std::string const& filename, MeshGrid3d const& grid3d, Real unit) {
    std::ofstream file(filename + ".csv");
    if (!file) {
        throw std::runtime_error("Failed to open file: " + filename + ".csv");
    }
    file << "i,j,k,value\n"; // CSV header
    const auto shape = grid3d.shape();
    for (size_t i = 0; i < shape[0]; ++i) {
        for (size_t j = 0; j < shape[1]; ++j) {
            for (size_t k = 0; k < shape[2]; ++k) {
                file << i << "," << j << "," << k << "," << (grid3d(i, j, k) / unit) << "\n";
            }
        }
    }
}

#ifndef NO_XTENSOR_IO
void write_npz(std::string const& filename, Coord const& coord) {
    xt::dump_npz(filename + ".npz", "t_src", xt::eval(coord.t / unit::sec), false, false);
    xt::dump_npz(filename + ".npz", "theta", coord.theta, false, true);
    xt::dump_npz(filename + ".npz", "phi", coord.phi, false, true);
}

void write_npz(std::string const& filename, SynPhotonGrid const& ph) {
    auto shape = ph.shape();
    xt::xarray<Real> nu_a = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> nu_m = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> nu_c = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> nu_M = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> I_nu_peak = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> regime = xt::zeros<Real>({shape[0], shape[1], shape[2]});

    for (size_t i = 0; i < shape[0]; ++i)
        for (size_t j = 0; j < shape[1]; ++j)
            for (size_t k = 0; k < shape[2]; ++k) {
                nu_a(i, j, k) = ph(i, j, k).nu_a / unit::Hz;
                nu_m(i, j, k) = ph(i, j, k).nu_m / unit::Hz;
                nu_c(i, j, k) = ph(i, j, k).nu_c / unit::Hz;
                nu_M(i, j, k) = ph(i, j, k).nu_M / unit::Hz;
                I_nu_peak(i, j, k) = ph(i, j, k).I_nu_max / (unit::erg / unit::cm2 / unit::sec / unit::Hz);
                regime(i, j, k) = ph(i, j, k).regime;
            }

    xt::dump_npz(filename + ".npz", "nu_a", nu_a, false, false);
    xt::dump_npz(filename + ".npz", "nu_m", nu_m, false, true);
    xt::dump_npz(filename + ".npz", "nu_c", nu_c, false, true);
    xt::dump_npz(filename + ".npz", "nu_Max", nu_M, false, true);
    xt::dump_npz(filename + ".npz", "I_nu_peak", I_nu_peak, false, true);
    xt::dump_npz(filename + ".npz", "regime", regime, false, true);
}

void write_npz(std::string const& filename, SynElectronGrid const& e) {
    auto shape = e.shape();
    xt::xarray<Real> gamma_a = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> gamma_m = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> gamma_c = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> gamma_M = xt::zeros<Real>({shape[0], shape[1], shape[2]});
    xt::xarray<Real> N_e = xt::zeros<Real>({shape[0], shape[1], shape[2]});

    for (size_t i = 0; i < shape[0]; ++i)
        for (size_t j = 0; j < shape[1]; ++j)
            for (size_t k = 0; k < shape[2]; ++k) {
                gamma_a(i, j, k) = e(i, j, k).gamma_a;
                gamma_m(i, j, k) = e(i, j, k).gamma_m;
                gamma_c(i, j, k) = e(i, j, k).gamma_c;
                gamma_M(i, j, k) = e(i, j, k).gamma_M;
                N_e(i, j, k) = e(i, j, k).N_e;
            }

    xt::dump_npz(filename + ".npz", "gamma_a", gamma_a, false, false);
    xt::dump_npz(filename + ".npz", "gamma_m", gamma_m, false, true);
    xt::dump_npz(filename + ".npz", "gamma_c", gamma_c, false, true);
    xt::dump_npz(filename + ".npz", "gamma_Max", gamma_M, false, true);
    xt::dump_npz(filename + ".npz", "N_e", N_e, false, true);
}

void write_npz(std::string const& filename, Shock const& shock) {
    xt::dump_npz(filename + ".npz", "Gamma", shock.Gamma, false, false);
    xt::dump_npz(filename + ".npz", "Gamma_th", shock.Gamma_th, false, true);
    xt::dump_npz(filename + ".npz", "B", xt::eval(shock.B / unit::Gauss), false, true);
    xt::dump_npz(filename + ".npz", "t_comv", xt::eval(shock.t_comv / unit::sec), false, true);
    xt::dump_npz(filename + ".npz", "r", xt::eval(shock.r / unit::cm), false, true);
    xt::dump_npz(filename + ".npz", "theta", shock.theta, false, true);
    xt::dump_npz(filename + ".npz", "N_p", xt::eval(shock.N_p), false, true);
}
#endif
