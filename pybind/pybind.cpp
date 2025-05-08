//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "mcmc.h"
#include "emission.h"

namespace py = pybind11;
PYBIND11_MODULE(VegasAfterglowC, m) {
    // Parameters for MCMC modeling
    py::class_<Params>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("E_iso", &Params::E_iso)
        .def_readwrite("Gamma0", &Params::Gamma0)
        .def_readwrite("theta_c", &Params::theta_c)
        .def_readwrite("theta_v", &Params::theta_v)
        .def_readwrite("theta_w", &Params::theta_w)
        .def_readwrite("p", &Params::p)
        .def_readwrite("p_r", &Params::p_r)
        .def_readwrite("eps_e", &Params::eps_e)
        .def_readwrite("eps_B", &Params::eps_B)
        .def_readwrite("eps_e_r", &Params::eps_e_r)
        .def_readwrite("eps_B_r", &Params::eps_B_r)
        .def_readwrite("n_ism", &Params::n_ism)
        .def_readwrite("A_star", &Params::A_star)
        .def_readwrite("xi", &Params::xi)
        .def_readwrite("xi_r", &Params::xi_r)
        .def_readwrite("k_jet", &Params::k_jet)
        .def("__repr__", [](const Params &p) {
            return "<Params E_iso=" + std::to_string(p.E_iso) + ", Gamma0=" + std::to_string(p.Gamma0) +
                ", theta_c=" + std::to_string(p.theta_c) + ", theta_v=" + std::to_string(p.theta_v) +
                ", theta_w=" + std::to_string(p.theta_w) + 
                ", p=" + std::to_string(p.p) + ", p_r=" + std::to_string(p.p_r) +
                ", eps_e=" + std::to_string(p.eps_e) + ", eps_e_r=" + std::to_string(p.eps_e_r) +
                ", eps_B=" + std::to_string(p.eps_B) + ", eps_B_r=" + std::to_string(p.eps_B_r) +
                ", n_ism=" + std::to_string(p.n_ism) + ", A_star=" + std::to_string(p.A_star) +
                ", xi=" + std::to_string(p.xi) + ", xi_r=" + std::to_string(p.xi_r) +
                ", k_jet=" + std::to_string(p.k_jet) + ">";
        });
    // Parameters for modeling that are not used in the MCMC
    py::class_<ConfigParams>(m, "Setups")
        .def(py::init<>())
        .def_readwrite("lumi_dist", &ConfigParams::lumi_dist)
        .def_readwrite("z", &ConfigParams::z)
        .def_readwrite("medium", &ConfigParams::medium)
        .def_readwrite("jet", &ConfigParams::jet)
        // Allow both f(theta, phi, r) and f(r) profiles
        .def_property("rho_profile",
            [](const ConfigParams &c) { return c.rho_profile; },
            [](ConfigParams &c, py::function f) {
                int nargs = f.attr("__code__").attr("co_argcount").cast<int>();
                if (nargs == 3) {
                    c.rho_profile = [f](double theta, double phi, double r) {
                        return f(theta, phi, r).cast<double>() * unit::g / unit::cm3;
                    };
                } else if (nargs == 1) {
                    c.rho_profile = [f](double /*theta*/, double /*phi*/, double r) {
                        return f(r).cast<double>() * unit::g / unit::cm3;
                    };
                } else {
                    throw std::invalid_argument("rho_profile must be of either f(r) or f(theta, phi, r)");
                }
            }
        )
        .def_property("m_profile",
            [](const ConfigParams &c) { return c.m_profile; },
            [](ConfigParams &c, py::function f) {
                int nargs = f.attr("__code__").attr("co_argcount").cast<int>();
                if (nargs == 3) {
                    c.m_profile = [f](double theta, double phi, double r) {
                        return f(theta, phi, r).cast<double>() * unit::g;
                    };
                } else if (nargs == 1) {
                    c.m_profile = [f](double /*theta*/, double /*phi*/, double r) {
                        return f(r).cast<double>() * unit::g;
                    };
                } else {
                    throw std::invalid_argument("m_profile must of either f(r) or f(theta, phi, r)");
                }
            }
        )
        .def_property("eps_k_profile",
            [](const ConfigParams &c) { return c.eps_k_profile; },
            [](ConfigParams &c, py::function f) {
                int nargs = f.attr("__code__").attr("co_argcount").cast<int>();
                if (nargs == 2) {
                    c.eps_k_profile = [f](double theta, double phi) {
                        return f(theta, phi).cast<double>() * unit::erg;
                    };
                } else if (nargs == 1) {
                    c.eps_k_profile = [f](double theta, double /*phi*/) {
                        return f(theta).cast<double>() * unit::erg;
                    };
                } else {
                    throw std::invalid_argument("eps_k_profile must take either 1 (theta) or 2 (theta, phi) arguments");
                }
            }
        )
        .def_property("Gamma0_profile",
            [](const ConfigParams &c) { return c.Gamma0_profile; },
            [](ConfigParams &c, py::function f) {
                int nargs = f.attr("__code__").attr("co_argcount").cast<int>();
                if (nargs == 2) {
                    c.Gamma0_profile = [f](double theta, double phi) {
                        return f(theta, phi).cast<double>();
                    };
                } else if (nargs == 1) {
                    c.Gamma0_profile = [f](double theta, double /*phi*/) {
                        return f(theta).cast<double>();
                    };
                } else {
                    throw std::invalid_argument("Gamma0_profile must take either 1 (theta) or 2 (theta, phi) arguments");
                }
            }
        )
        .def_readwrite("T0", &ConfigParams::T0)
        .def_readwrite("t_grid", &ConfigParams::t_grid)
        .def_readwrite("phi_grid", &ConfigParams::phi_grid)
        .def_readwrite("theta_grid", &ConfigParams::theta_grid)
        .def_readwrite("rtol", &ConfigParams::rtol)
        .def_readwrite("spreading", &ConfigParams::spreading)
        .def("__repr__", [](const ConfigParams &c) {
            return "<ConfigParams lumi_dist=" + std::to_string(c.lumi_dist) + ", z=" + std::to_string(c.z) +
                   ", medium='" + c.medium + "', jet='" + c.jet + "', t_grid_num=" + std::to_string(c.t_grid) +
                   ", phi_grid_num=" + std::to_string(c.phi_grid) + ", theta_grid_num=" + std::to_string(c.theta_grid) +
                   ", rtol=" + std::to_string(c.rtol) + ", spreading=" + std::to_string(c.spreading) + ">";
        });

    // MultiBandData bindings
    py::class_<MultiBandData>(m, "ObsData")
        .def(py::init<>())
        .def("add_light_curve", &MultiBandData::add_light_curve, py::arg("nu_cgs"), py::arg("t_cgs"),
             py::arg("Fnu_cgs"), py::arg("Fnu_err"))
        .def("add_spectrum", &MultiBandData::add_spectrum, py::arg("t_cgs"), py::arg("nu_cgs"), py::arg("Fnu_cgs"),
             py::arg("Fnu_err"));

    // MultiBandModel bindings
    py::class_<MultiBandModel>(m, "VegasMC")
        .def(py::init<MultiBandData const &>(), py::arg("obs_data"))
        .def("set", &MultiBandModel::configure, py::arg("param"))
        .def("estimate_chi2", &MultiBandModel::estimate_chi2, py::arg("param"),
             py::call_guard<py::gil_scoped_release>())
        .def("light_curves", &MultiBandModel::light_curves, py::arg("param"), py::arg("t_cgs"), py::arg("nu_cgs"),
             py::call_guard<py::gil_scoped_release>())
        .def("spectra", &MultiBandModel::spectra, py::arg("param"), py::arg("nu_cgs"), py::arg("t_cgs"),
             py::call_guard<py::gil_scoped_release>());

    py::class_<Emission>(m, "Emission")
        .def(py::init<const Params&, const ConfigParams&>())
        .def("lc", &Emission::lc, py::arg("nu"), py::arg("t"))
        .def("lc_r", &Emission::lc_r, py::arg("nu"), py::arg("t"))
        .def("spec", &Emission::spec, py::arg("nu"), py::arg("t"));

}