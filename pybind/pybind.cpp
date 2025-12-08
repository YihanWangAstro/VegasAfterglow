//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#define FORCE_IMPORT_ARRAY // numpy C api loading must before any xtensor-python headers
#include "pybind.h"

#include "mcmc.h"
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
             py::call_guard<py::gil_scoped_release>());

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

    //========================================================================================================
    //                                 MCMC bindings
    //========================================================================================================
    // Parameters for MCMC modeling
    py::class_<Params>(m, "ModelParams")
        .def(py::init<>())
        .def_readwrite("theta_v", &Params::theta_v)

        .def_readwrite("n_ism", &Params::n_ism)
        .def_readwrite("n0", &Params::n0)
        .def_readwrite("A_star", &Params::A_star)
        .def_readwrite("k_m", &Params::k_m)

        .def_readwrite("E_iso", &Params::E_iso)
        .def_readwrite("Gamma0", &Params::Gamma0)
        .def_readwrite("theta_c", &Params::theta_c)
        .def_readwrite("k_e", &Params::k_e)
        .def_readwrite("k_g", &Params::k_g)
        .def_readwrite("tau", &Params::duration)

        .def_readwrite("E_iso_w", &Params::E_iso_w)
        .def_readwrite("Gamma0_w", &Params::Gamma0_w)
        .def_readwrite("theta_w", &Params::theta_w)

        .def_readwrite("L0", &Params::L0)
        .def_readwrite("t0", &Params::t0)
        .def_readwrite("q", &Params::q)

        .def_readwrite("p", &Params::p)
        .def_readwrite("eps_e", &Params::eps_e)
        .def_readwrite("eps_B", &Params::eps_B)
        .def_readwrite("xi_e", &Params::xi_e)

        .def_readwrite("p_r", &Params::p_r)
        .def_readwrite("eps_e_r", &Params::eps_e_r)
        .def_readwrite("eps_B_r", &Params::eps_B_r)
        .def_readwrite("xi_e_r", &Params::xi_e_r);

    // Parameters for modeling that are not used in the MCMC
    py::class_<ConfigParams>(m, "Setups")
        .def(py::init<>())
        .def_readwrite("lumi_dist", &ConfigParams::lumi_dist)
        .def_readwrite("z", &ConfigParams::z)
        .def_readwrite("medium", &ConfigParams::medium)
        .def_readwrite("jet", &ConfigParams::jet)
        .def_readwrite("t_resol", &ConfigParams::t_resol)
        .def_readwrite("phi_resol", &ConfigParams::phi_resol)
        .def_readwrite("theta_resol", &ConfigParams::theta_resol)
        .def_readwrite("rtol", &ConfigParams::rtol)
        .def_readwrite("rvs_shock", &ConfigParams::rvs_shock)
        .def_readwrite("fwd_ssc", &ConfigParams::fwd_ssc)
        .def_readwrite("rvs_ssc", &ConfigParams::rvs_ssc)
        .def_readwrite("ssc_cooling", &ConfigParams::ssc_cooling)
        .def_readwrite("kn", &ConfigParams::kn)
        .def_readwrite("magnetar", &ConfigParams::magnetar)
        .def(py::pickle(
            [](const ConfigParams& self) {
                return py::make_tuple(self.lumi_dist, self.z, self.medium, self.jet, self.t_resol, self.phi_resol,
                                      self.theta_resol, self.rtol, self.rvs_shock, self.fwd_ssc, self.rvs_ssc,
                                      self.ssc_cooling, self.kn, self.magnetar);
            },
            [](py::tuple state) {
                ConfigParams cfg;
                cfg.lumi_dist = state[0].cast<double>();
                cfg.z = state[1].cast<double>();
                cfg.medium = state[2].cast<std::string>();
                cfg.jet = state[3].cast<std::string>();
                cfg.t_resol = state[4].cast<Real>();
                cfg.phi_resol = state[5].cast<Real>();
                cfg.theta_resol = state[6].cast<Real>();
                cfg.rtol = state[7].cast<double>();
                cfg.rvs_shock = state[8].cast<bool>();
                cfg.fwd_ssc = state[9].cast<bool>();
                cfg.rvs_ssc = state[10].cast<bool>();
                cfg.ssc_cooling = state[11].cast<bool>();
                cfg.kn = state[12].cast<bool>();
                cfg.magnetar = state[13].cast<bool>();
                return cfg;
            }));

    // FluxData bindings (for integrated flux bands)
    py::class_<FluxData>(m, "FluxData")
        .def(py::init<>())
        .def_property_readonly("t", [](const FluxData& self) { return PyArray(self.t); })
        .def_property_readonly("nu", [](const FluxData& self) { return PyArray(self.nu); })
        .def_property_readonly("Fv_obs", [](const FluxData& self) { return PyArray(self.Fv_obs); })
        .def_property_readonly("Fv_err", [](const FluxData& self) { return PyArray(self.Fv_err); })
        .def_property_readonly("weights", [](const FluxData& self) { return PyArray(self.weights); })
        .def(py::pickle(
            [](const FluxData& self) {
                return py::make_tuple(PyArray(self.t), PyArray(self.nu), PyArray(self.Fv_obs), PyArray(self.Fv_err),
                                      PyArray(self.weights));
            },
            [](py::tuple state) {
                FluxData fd;
                PyArray py_t = state[0].cast<PyArray>();
                PyArray py_nu = state[1].cast<PyArray>();
                PyArray py_Fv_obs = state[2].cast<PyArray>();
                PyArray py_Fv_err = state[3].cast<PyArray>();
                PyArray py_weights = state[4].cast<PyArray>();

                fd.t = xt::eval(py_t);
                fd.nu = xt::eval(py_nu);
                fd.Fv_obs = xt::eval(py_Fv_obs);
                fd.Fv_err = xt::eval(py_Fv_err);
                fd.Fv_model = xt::zeros<Real>({fd.t.size()});
                fd.weights = xt::eval(py_weights);
                return fd;
            }));

    // MultiBandData bindings
    py::class_<MultiBandData>(m, "ObsData")
        .def(py::init<>())

        .def("add_flux_density", &MultiBandData::add_flux_density, py::arg("nu"), py::arg("t"), py::arg("f_nu"),
             py::arg("err"), py::arg("weights") = py::none(),
             "Add flux density data. All quantities in CGS units: nu [Hz], t [s], f_nu [erg/cm²/s/Hz]")

        .def("add_flux", &MultiBandData::add_flux, py::arg("nu_min"), py::arg("nu_max"), py::arg("num_points"),
             py::arg("t"), py::arg("flux"), py::arg("err"), py::arg("weights") = py::none(),
             "Add broadband flux data. All quantities in CGS units: nu [Hz], t [s], flux [erg/cm²/s]")

        .def("add_spectrum", &MultiBandData::add_spectrum, py::arg("t"), py::arg("nu"), py::arg("f_nu"), py::arg("err"),
             py::arg("weights") = py::none(),
             "Add spectrum data. All quantities in CGS units: t [s], nu [Hz], f_nu [erg/cm²/s/Hz]")

        .def_static("logscale_screen", &MultiBandData::logscale_screen, py::arg("t"), py::arg("data_density"))

        .def("data_points_num", &MultiBandData::data_points_num)

        // Expose internal data as numpy arrays (read-only)
        .def_property_readonly("times", [](const MultiBandData& self) { return PyArray(self.times); })
        .def_property_readonly("frequencies", [](const MultiBandData& self) { return PyArray(self.frequencies); })
        .def_property_readonly("fluxes", [](const MultiBandData& self) { return PyArray(self.fluxes); })
        .def_property_readonly("errors", [](const MultiBandData& self) { return PyArray(self.errors); })
        .def_property_readonly("weights", [](const MultiBandData& self) { return PyArray(self.weights); })
        .def_readonly("t_min", &MultiBandData::t_min)
        .def_readonly("t_max", &MultiBandData::t_max)
        .def_readonly("flux_data", &MultiBandData::flux_data)

        // Pickle support
        .def(py::pickle(
            [](MultiBandData& self) {
                // Ensure data arrays are filled before pickling
                self.fill_data_arrays();

                // Serialize consolidated arrays and flux_data
                py::list flux_data_list;
                for (const auto& fd : self.flux_data) {
                    flux_data_list.append(py::make_tuple(PyArray(fd.t), PyArray(fd.nu), PyArray(fd.Fv_obs),
                                                         PyArray(fd.Fv_err), PyArray(fd.weights)));
                }
                return py::make_tuple(PyArray(self.times), PyArray(self.frequencies), PyArray(self.fluxes),
                                      PyArray(self.errors), PyArray(self.weights), self.t_min, self.t_max,
                                      flux_data_list);
            },
            [](py::tuple state) {
                MultiBandData data;
                // Explicitly copy data from pytensor to xtensor using xt::eval to ensure deep copy
                PyArray py_times = state[0].cast<PyArray>();
                PyArray py_frequencies = state[1].cast<PyArray>();
                PyArray py_fluxes = state[2].cast<PyArray>();
                PyArray py_errors = state[3].cast<PyArray>();
                PyArray py_weights = state[4].cast<PyArray>();

                data.times = xt::eval(py_times);
                data.frequencies = xt::eval(py_frequencies);
                data.fluxes = xt::eval(py_fluxes);
                data.errors = xt::eval(py_errors);
                data.weights = xt::eval(py_weights);
                data.model_fluxes = xt::zeros<Real>({data.times.size()});
                data.t_min = state[5].cast<double>();
                data.t_max = state[6].cast<double>();

                py::list flux_data_list = state[7].cast<py::list>();
                for (auto item : flux_data_list) {
                    auto tup = item.cast<py::tuple>();
                    FluxData fd;
                    PyArray py_t = tup[0].cast<PyArray>();
                    PyArray py_nu = tup[1].cast<PyArray>();
                    PyArray py_Fv_obs = tup[2].cast<PyArray>();
                    PyArray py_Fv_err = tup[3].cast<PyArray>();
                    PyArray py_fd_weights = tup[4].cast<PyArray>();

                    fd.t = xt::eval(py_t);
                    fd.nu = xt::eval(py_nu);
                    fd.Fv_obs = xt::eval(py_Fv_obs);
                    fd.Fv_err = xt::eval(py_Fv_err);
                    fd.Fv_model = xt::zeros<Real>({fd.t.size()});
                    fd.weights = xt::eval(py_fd_weights);
                    data.flux_data.push_back(std::move(fd));
                }
                return data;
            }));

    // MultiBandModel bindings
    py::class_<MultiBandModel>(m, "VegasMC")
        .def(py::init<MultiBandData const&>(), py::arg("obs_data"))

        .def("set", &MultiBandModel::configure, py::arg("param"))

        .def("estimate_chi2", &MultiBandModel::estimate_chi2, py::arg("param"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux_density_grid", &MultiBandModel::flux_density_grid, py::arg("param"), py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux", &MultiBandModel::flux, py::arg("param"), py::arg("t"), py::arg("nu_min"), py::arg("nu_max"),
             py::arg("num_points"), py::call_guard<py::gil_scoped_release>());
}
