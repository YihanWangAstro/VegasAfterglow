//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#define FORCE_IMPORT_ARRAY // numpy C api loading must before any xtensor-python headers
#include "pybind.h"

#include "pymodel.h"

PYBIND11_MODULE(VegasAfterglowC, m) {
    xt::import_numpy();
    // Jet bindings
    py::object zero2d_fn = py::cpp_function(func::zero_2d);
    py::object zero3d_fn = py::cpp_function(func::zero_3d);

    //========================================================================================================
    //                                 Model bindings
    //========================================================================================================
    py::class_<PyMagnetar>(m, "Magnetar").def(py::init<Real, Real, Real>(), py::arg("L0"), py::arg("t0"), py::arg("q"));

    m.def("TophatJet", &PyTophatJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("GaussianJet", &PyGaussianJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("PowerLawJet", &PyPowerLawJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"), py::arg("k_e"),
          py::arg("k_g"), py::arg("spreading") = false, py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("PowerLawWing", &PyPowerLawWing, py::arg("theta_c"), py::arg("E_iso_w"), py::arg("Gamma0_w"), py::arg("k_e"),
          py::arg("k_g"), py::arg("spreading") = false, py::arg("duration") = 1);

    m.def("TwoComponentJet", &PyTwoComponentJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("theta_w"), py::arg("E_iso_w"), py::arg("Gamma0_w"), py::arg("spreading") = false,
          py::arg("duration") = 1, py::arg("magnetar") = py::none());

    m.def("StepPowerLawJet", &PyStepPowerLawJet, py::arg("theta_c"), py::arg("E_iso"), py::arg("Gamma0"),
          py::arg("E_iso_w"), py::arg("Gamma0_w"), py::arg("k_e"), py::arg("k_g"), py::arg("spreading") = false,
          py::arg("duration") = 1, py::arg("magnetar") = py::none());

    py::class_<Ejecta>(m, "Ejecta")
        .def(py::init<BinaryFunc, BinaryFunc, BinaryFunc, TernaryFunc, TernaryFunc, bool, Real>(), py::arg("E_iso"),
             py::arg("Gamma0"), py::arg("sigma0") = zero2d_fn, py::arg("E_dot") = zero3d_fn,
             py::arg("M_dot") = zero3d_fn, py::arg("spreading") = false, py::arg("duration") = 1);

    // Medium bindings
    m.def("ISM", &PyISM, py::arg("n_ism"));

    m.def("Wind", &PyWind, py::arg("A_star"), py::arg("n_ism") = 0, py::arg("n0") = con::inf, py::arg("k") = 2);

    py::class_<Medium>(m, "Medium").def(py::init<TernaryFunc>(), py::arg("rho"));

    // Observer bindings
    py::class_<PyObserver>(m, "Observer")
        .def(py::init<Real, Real, Real, Real>(), py::arg("lumi_dist"), py::arg("z"), py::arg("theta_obs"),
             py::arg("phi_obs") = 0);

    // Radiation bindings
    py::class_<PyRadiation>(m, "Radiation")
        .def(py::init<Real, Real, Real, Real, bool, bool, bool>(), py::arg("eps_e"), py::arg("eps_B"), py::arg("p"),
             py::arg("xi_e") = 1, py::arg("ssc_cooling") = false, py::arg("ssc") = false, py::arg("kn") = false);

    // Model bindings
    py::class_<PyModel>(m, "Model")
        .def(py::init<Ejecta, Medium, PyObserver, PyRadiation, std::optional<PyRadiation>, std::tuple<Real, Real, Real>,
                      Real, bool>(),
             py::arg("jet"), py::arg("medium"), py::arg("observer"), py::arg("fwd_rad"),
             py::arg("rvs_rad") = py::none(), py::arg("resolutions") = std::make_tuple(0.3, 1, 10),
             py::arg("rtol") = 1e-6, py::arg("axisymmetric") = true)

        .def("flux_density_grid", &PyModel::flux_density_grid, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux_density", &PyModel::flux_density, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux", &PyModel::flux, py::arg("t"), py::arg("nu_min"), py::arg("nu_max"), py::arg("num_nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux_density_exposures", &PyModel::flux_density_exposures, py::arg("t"), py::arg("nu"),
             py::arg("expo_time"), py::arg("num_points") = 10, py::call_guard<py::gil_scoped_release>())

        .def("details", &PyModel::details, py::arg("t_min"), py::arg("t_max"), py::call_guard<py::gil_scoped_release>())

        .def("medium", &PyModel::medium, py::arg("phi"), py::arg("theta"), py::arg("r"),
             py::call_guard<py::gil_scoped_release>())

        .def("jet_E_iso", &PyModel::jet_E_iso, py::arg("phi"), py::arg("theta"),
             py::call_guard<py::gil_scoped_release>())

        .def("jet_Gamma0", &PyModel::jet_Gamma0, py::arg("phi"), py::arg("theta"),
             py::call_guard<py::gil_scoped_release>())
#ifdef AFTERGLOW_PROFILE
        .def_static("profile_data", &PyModel::profile_data, "Get per-stage timing from the last computation (ms)")
        .def_static("profile_reset", &PyModel::profile_reset, "Reset profiling counters")
#endif
        ;

    py::class_<Flux>(m, "Flux").def(py::init<>()).def_readonly("sync", &Flux::sync).def_readonly("ssc", &Flux::ssc);

    py::class_<PyFlux>(m, "FluxDict")
        .def(py::init<>())
        .def_readonly("fwd", &PyFlux::fwd)
        .def_readonly("rvs", &PyFlux::rvs)
        .def_readonly("total", &PyFlux::total);

    py::class_<PyShock>(m, "ShockDetails")
        .def(py::init<>())
        .def_readonly("t_comv", &PyShock::t_comv)
        .def_readonly("t_obs", &PyShock::t_obs)
        .def_readonly("Gamma", &PyShock::Gamma)
        .def_readonly("Gamma_th", &PyShock::Gamma_th)
        .def_readonly("B_comv", &PyShock::B_comv)
        .def_readonly("r", &PyShock::r)
        .def_readonly("theta", &PyShock::theta)
        .def_readonly("N_p", &PyShock::N_p)
        .def_readonly("N_e", &PyShock::N_e)
        .def_readonly("gamma_m", &PyShock::gamma_m)
        .def_readonly("gamma_c", &PyShock::gamma_c)
        .def_readonly("gamma_M", &PyShock::gamma_M)
        .def_readonly("gamma_a", &PyShock::gamma_a)
        .def_readonly("gamma_m_hat", &PyShock::gamma_m_hat)
        .def_readonly("gamma_c_hat", &PyShock::gamma_c_hat)
        .def_readonly("nu_m", &PyShock::nu_m)
        .def_readonly("nu_c", &PyShock::nu_c)
        .def_readonly("nu_M", &PyShock::nu_M)
        .def_readonly("nu_a", &PyShock::nu_a)
        .def_readonly("nu_m_hat", &PyShock::nu_m_hat)
        .def_readonly("nu_c_hat", &PyShock::nu_c_hat)
        .def_readonly("Y_T", &PyShock::Y_T)
        .def_readonly("I_nu_max", &PyShock::I_nu_max)
        .def_readonly("Doppler", &PyShock::Doppler);

    py::class_<PyDetails>(m, "SimulationDetails")
        .def(py::init<>())
        .def_readonly("phi", &PyDetails::phi)
        .def_readonly("theta", &PyDetails::theta)
        .def_readonly("t_src", &PyDetails::t_src)
        .def_readonly("fwd", &PyDetails::fwd)
        .def_readonly("rvs", &PyDetails::rvs);

    // Pure physics models only - MCMC/inference removed
    // For parameter estimation, use redback: https://github.com/nikhil-sarin/redback
}
