//              __     __                            _      __  _                     _
//              \ \   / /___   __ _   __ _  ___     / \    / _|| |_  ___  _ __  __ _ | |  ___ __      __
//               \ \ / // _ \ / _` | / _` |/ __|   / _ \  | |_ | __|/ _ \| '__|/ _` || | / _ \\ \ /\ / /
//                \ V /|  __/| (_| || (_| |\__ \  / ___ \ |  _|| |_|  __/| |  | (_| || || (_) |\ V  V /
//                 \_/  \___| \__, | \__,_||___/ /_/   \_\|_|   \__|\___||_|   \__, ||_| \___/  \_/\_/
//                            |___/                                            |___/

#define FORCE_IMPORT_ARRAY // numpy C api loading must before any xtensor-python headers
#include "pybind.h"

#include <array>
#include <cstdint>
#include <cstdio>

#include "pymodel.h"

// ---------------------------------------------------------------------------
// NativeFunc support: GIL-free callables from numba @cfunc + bound params
// ---------------------------------------------------------------------------
// Uses C++20 template lambdas to build tight (fn + params) lambdas.
// The index pack simultaneously deduces the function pointer type via
// decltype((void)I, Real{})..., initializes params, and unpacks the call.
constexpr size_t MAX_NATIVE_PARAMS = 10;

// ---------------------------------------------------------------------------
// Callback exception guard: user Python callbacks are invoked deep inside
// noexcept ODE right-hand sides, where a propagating py::error_already_set
// would std::terminate the interpreter. The wrappers below never throw --
// they capture the error thread-locally, return NaN (poisoning the solve),
// and the PyModel entry points rethrow before returning to Python.
// ---------------------------------------------------------------------------
namespace callback_guard {
    inline std::exception_ptr& slot() {
        thread_local std::exception_ptr err;
        return err;
    }

    inline Real capture_current() noexcept {
        try {
            py::gil_scoped_acquire gil;
            try {
                throw;
            } catch (py::error_already_set& e) {
                slot() = std::make_exception_ptr(
                    std::runtime_error(std::string("user-defined jet/medium callback raised:\n") + e.what()));
            } catch (std::exception& e) {
                slot() = std::make_exception_ptr(
                    std::runtime_error(std::string("user-defined jet/medium callback failed: ") + e.what()));
            } catch (...) {
                slot() = std::make_exception_ptr(std::runtime_error("user-defined jet/medium callback failed"));
            }
        } catch (...) { // never propagate from a noexcept context
        }
        return std::numeric_limits<Real>::quiet_NaN();
    }

    template <typename F>
    BinaryFunc guard_binary(F f) {
        return [f = std::move(f)](Real a, Real b) noexcept -> Real {
            try {
                return f(a, b);
            } catch (...) {
                return capture_current();
            }
        };
    }

    template <typename F>
    TernaryFunc guard_ternary(F f) {
        return [f = std::move(f)](Real a, Real b, Real c) noexcept -> Real {
            try {
                return f(a, b, c);
            } catch (...) {
                return capture_current();
            }
        };
    }
} // namespace callback_guard

void rethrow_callback_error() {
    if (callback_guard::slot()) {
        auto err = callback_guard::slot();
        callback_guard::slot() = nullptr;
        std::rethrow_exception(err);
    }
}

template <size_t N>
BinaryFunc make_native_binary(uintptr_t addr, std::vector<Real> const& p) {
    return [&]<size_t... I>(std::index_sequence<I...>) -> BinaryFunc {
        auto fn = reinterpret_cast<Real (*)(Real, Real, decltype((void)I, Real{})...)>(addr);
        std::array<Real, N> params = {p[I]...};
        return [fn, params](Real phi, Real theta) noexcept -> Real {
            (void)params;
            return fn(phi, theta, params[I]...);
        };
    }(std::make_index_sequence<N>{});
}

template <size_t... Ns>
BinaryFunc dispatch_binary(uintptr_t addr, std::vector<Real> const& p, std::index_sequence<Ns...>) {
    using maker_t = BinaryFunc (*)(uintptr_t, std::vector<Real> const&);
    static constexpr maker_t table[] = {make_native_binary<Ns>...};
    if (p.size() >= sizeof...(Ns)) {
        throw std::runtime_error("NativeFunc: too many bound parameters (max " + std::to_string(sizeof...(Ns) - 1) +
                                 ")");
    }
    return table[p.size()](addr, p);
}

template <size_t N>
TernaryFunc make_native_ternary(uintptr_t addr, std::vector<Real> const& p) {
    return [&]<size_t... I>(std::index_sequence<I...>) -> TernaryFunc {
        auto fn = reinterpret_cast<Real (*)(Real, Real, Real, decltype((void)I, Real{})...)>(addr);
        std::array<Real, N> params = {p[I]...};
        return [fn, params](Real phi, Real theta, Real r) noexcept -> Real {
            (void)params;
            return fn(phi, theta, r, params[I]...);
        };
    }(std::make_index_sequence<N>{});
}

template <size_t... Ns>
TernaryFunc dispatch_ternary(uintptr_t addr, std::vector<Real> const& p, std::index_sequence<Ns...>) {
    using maker_t = TernaryFunc (*)(uintptr_t, std::vector<Real> const&);
    static constexpr maker_t table[] = {make_native_ternary<Ns>...};
    if (p.size() >= sizeof...(Ns)) {
        throw std::runtime_error("NativeFunc: too many bound parameters (max " + std::to_string(sizeof...(Ns) - 1) +
                                 ")");
    }
    return table[p.size()](addr, p);
}

// --- Converters: detect NativeFunc and dispatch to GIL-free path ---
// A NativeFunc's cfunc signature is (runtime args..., bound params...); the
// raw function pointer is reinterpret_cast to that exact arity, so a mismatch
// (e.g. a two-argument cfunc used as a medium rho) would silently read bound
// parameters from the wrong registers. Validate the total argument count.
void check_native_arity(py::object const& obj, size_t runtime_args, size_t n_params, char const* what) {
    const auto n_args = obj.attr("n_args").cast<size_t>();
    AFTERGLOW_REQUIRE(n_args == runtime_args + n_params, std::string(what) + ": cfunc takes " + std::to_string(n_args) +
                                                             " arguments but " + std::to_string(runtime_args) +
                                                             " runtime argument(s) + " + std::to_string(n_params) +
                                                             " bound parameter(s) were provided");
}

BinaryFunc to_binary_func(py::object const& obj, py::object const& native_type) {
    if (!native_type.is_none() && py::isinstance(obj, native_type)) {
        const auto addr = obj.attr("address").cast<uintptr_t>();
        const auto params = obj.attr("params").cast<std::vector<Real>>();
        check_native_arity(obj, 2, params.size(), "NativeFunc used as (phi, theta) function");
        return dispatch_binary(addr, params, std::make_index_sequence<MAX_NATIVE_PARAMS + 1>{});
    }
    return callback_guard::guard_binary(obj.cast<BinaryFunc>());
}

TernaryFunc to_ternary_func(py::object const& obj, py::object const& native_type) {
    if (!native_type.is_none() && py::isinstance(obj, native_type)) {
        const auto addr = obj.attr("address").cast<uintptr_t>();
        const auto params = obj.attr("params").cast<std::vector<Real>>();
        check_native_arity(obj, 3, params.size(), "NativeFunc used as (phi, theta, r) or (phi, theta, t) function");
        return dispatch_ternary(addr, params, std::make_index_sequence<MAX_NATIVE_PARAMS + 1>{});
    }
    return callback_guard::guard_ternary(obj.cast<TernaryFunc>());
}

/// The VegasAfterglow.native NativeFunc type, resolved once per process (py::none() when the
/// module is unavailable). Deliberately leaked: a static py::object would be destroyed after
/// interpreter finalization.
static py::handle native_func_type() {
    static py::handle type = [] {
        try {
            py::object t = py::module_::import("VegasAfterglow.native").attr("NativeFunc");
            return t.release();
        } catch (...) {
            return py::handle();
        }
    }();
    return type;
}

PYBIND11_MODULE(VegasAfterglowC, m) {
    xt::import_numpy();

    // Build-flavor introspection: golden baselines must be generated with the
    // exact-libm build (AFTERGLOW_FAST_MATH off); regenerate.py checks this.
#ifdef AFTERGLOW_FAST_MATH
    m.attr("fast_math_enabled") = true;
#else
    m.attr("fast_math_enabled") = false;
#endif

    // Jet bindings

    //========================================================================================================
    //                                 Model bindings
    //========================================================================================================
    py::class_<PyMagnetar>(m, "Magnetar")
        .def(py::init<Real, Real, Real>(), py::arg("L0"), py::arg("t0"), py::arg("q"))
        .def_readonly("L0", &PyMagnetar::L0)
        .def_readonly("t0", &PyMagnetar::t0)
        .def_readonly("q", &PyMagnetar::q)
        .def("__repr__", &PyMagnetar::repr);

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
        .def(
            py::init([](py::object E_iso_obj, py::object Gamma0_obj, py::object sigma0_obj, py::object E_dot_obj,
                        py::object M_dot_obj, bool spreading, Real duration) {
                AFTERGLOW_REQUIRE(std::isfinite(duration) && duration > 0,
                                  std::string("duration must be positive and finite, got ") + std::to_string(duration));
                AFTERGLOW_REQUIRE(PyCallable_Check(E_iso_obj.ptr()), "E_iso must be callable: (phi, theta) -> erg");
                AFTERGLOW_REQUIRE(PyCallable_Check(Gamma0_obj.ptr()),
                                  "Gamma0 must be callable: (phi, theta) -> Lorentz factor");
                AFTERGLOW_REQUIRE(sigma0_obj.is_none() || PyCallable_Check(sigma0_obj.ptr()),
                                  "sigma0 must be callable: (phi, theta) -> magnetization");
                AFTERGLOW_REQUIRE(E_dot_obj.is_none() || PyCallable_Check(E_dot_obj.ptr()),
                                  "E_dot must be callable: (phi, theta, t) -> erg/s");
                AFTERGLOW_REQUIRE(M_dot_obj.is_none() || PyCallable_Check(M_dot_obj.ptr()),
                                  "M_dot must be callable: (phi, theta, t) -> g/s");

                const py::object native_type =
                    native_func_type() ? py::reinterpret_borrow<py::object>(native_func_type()) : py::none();

                auto eps_k = to_binary_func(E_iso_obj, native_type);
                auto Gamma0 = to_binary_func(Gamma0_obj, native_type);
                // Default sigma0/E_dot/M_dot to C++ zero functions (GIL-free)
                auto sigma0 =
                    sigma0_obj.is_none() ? BinaryFunc(func::zero_2d) : to_binary_func(sigma0_obj, native_type);
                auto E_dot = E_dot_obj.is_none() ? TernaryFunc(func::zero_3d) : to_ternary_func(E_dot_obj, native_type);
                auto M_dot = M_dot_obj.is_none() ? TernaryFunc(func::zero_3d) : to_ternary_func(M_dot_obj, native_type);

                // Probe the profiles on-axis to reject obviously broken functions at construction
                // time (same contract as the typed jet factories, which validate their scalars).
                const Real eps_probe = eps_k(0, 0);
                rethrow_callback_error();
                AFTERGLOW_REQUIRE(std::isfinite(eps_probe) && eps_probe >= 0,
                                  std::string("E_iso(phi=0, theta=0) must be finite and non-negative, got ") +
                                      std::to_string(eps_probe));
                const Real Gamma_probe = Gamma0(0, 0);
                rethrow_callback_error();
                AFTERGLOW_REQUIRE(std::isfinite(Gamma_probe) && Gamma_probe >= 1,
                                  std::string("Gamma0(phi=0, theta=0) must be finite and >= 1, got ") +
                                      std::to_string(Gamma_probe));
                const Real sigma_probe = sigma0(0, 0);
                AFTERGLOW_REQUIRE(std::isfinite(sigma_probe) && sigma_probe >= 0,
                                  std::string("sigma0(phi=0, theta=0) must be finite and non-negative, got ") +
                                      std::to_string(sigma_probe));

                return Ejecta(std::move(eps_k), std::move(Gamma0), std::move(sigma0), std::move(E_dot),
                              std::move(M_dot), spreading, duration);
            }),
            py::arg("E_iso"), py::arg("Gamma0"), py::arg("sigma0") = py::none(), py::arg("E_dot") = py::none(),
            py::arg("M_dot") = py::none(), py::arg("spreading") = false, py::arg("duration") = 1);

    // Jet bindings — register concrete types for std::variant support
    auto tophat_type [[maybe_unused]] = py::class_<TophatJet>(m, "_TophatJet");
    auto gaussian_type [[maybe_unused]] = py::class_<GaussianJet>(m, "_GaussianJet");
    auto powerlaw_type [[maybe_unused]] = py::class_<PowerLawJet>(m, "_PowerLawJet");

    // Medium bindings — register concrete types for std::variant support
    auto ism_type [[maybe_unused]] = py::class_<ISM>(m, "_ISM");
    auto wind_type [[maybe_unused]] = py::class_<Wind>(m, "_Wind");
    py::class_<Medium>(m, "Medium")
        .def(py::init([](py::object rho_obj, bool isotropic) {
                 AFTERGLOW_REQUIRE(PyCallable_Check(rho_obj.ptr()), "rho must be callable: (phi, theta, r) -> g/cm^3");
                 const py::object native_type =
                     native_func_type() ? py::reinterpret_borrow<py::object>(native_func_type()) : py::none();
                 auto rho = to_ternary_func(rho_obj, native_type);

                 // Probe the density on-axis at a typical afterglow radius to reject obviously
                 // broken functions at construction time (mirrors the typed ISM/Wind factories).
                 const Real rho_probe = rho(0, 0, 1e17);
                 rethrow_callback_error();
                 AFTERGLOW_REQUIRE(std::isfinite(rho_probe) && rho_probe >= 0,
                                   std::string("rho(phi=0, theta=0, r=1e17 cm) must be finite and non-negative, "
                                               "got ") +
                                       std::to_string(rho_probe));

                 // ``isotropic`` is a required declaration: an anisotropic medium forces
                 // per-cell blast-wave dynamics (no angular broadcast, an order-of-magnitude
                 // cost), while wrongly declaring an angular profile isotropic would give
                 // wrong physics. Sampling cannot *prove* isotropy for a black-box function,
                 // but it can disprove it — so a True declaration is cross-checked against
                 // an angular probe grid and rejected when the profile visibly varies.
                 if (isotropic) {
                     constexpr std::array thetas = {0.11, 0.31, 0.52, 0.73, 0.94, 1.15, 1.36, 1.57};
                     constexpr std::array phis = {0.7, 1.9, 3.4, 5.3};
                     constexpr std::array radii = {1e13, 1e14, 1e15, 1e16, 1e17, 1e18, 1e19, 1e20};
                     auto fmt = [](Real v) {
                         char buf[32];
                         std::snprintf(buf, sizeof(buf), "%.6g", v);
                         return std::string(buf);
                     };
                     for (const Real r : radii) {
                         const Real ref = rho(0, 0, r);
                         auto check = [&](Real phi, Real theta) {
                             const Real val = rho(phi, theta, r);
                             rethrow_callback_error();
                             AFTERGLOW_REQUIRE(
                                 val == ref,
                                 std::string("Medium was declared isotropic=True but rho varies with angle: "
                                             "rho(phi=") +
                                     fmt(phi) + ", theta=" + fmt(theta) + ", r=" + fmt(r) + " cm) = " + fmt(val) +
                                     " != rho(0, 0, r) = " + fmt(ref) +
                                     ". Pass isotropic=False for angle-dependent density profiles.");
                         };
                         for (const Real theta : thetas) {
                             check(0, theta);
                         }
                         for (const Real phi : phis) {
                             check(phi, 0);
                         }
                         // joint samples catch profiles separable in neither angle alone
                         check(1.3, 0.62);
                         check(4.1, 1.24);
                     }
                 }

                 return Medium(std::move(rho), isotropic);
             }),
             py::arg("rho"), py::arg("isotropic"),
             "Custom density profile rho(phi, theta, r[cm]) -> g/cm^3. ``isotropic`` must be\n"
             "declared explicitly: True keeps the fast broadcast dynamics (the declaration is\n"
             "cross-checked against an angular probe grid), False solves every angular cell.");

    // Factory functions return MediumVariant (ISM/Wind for optimized path, Medium for fallback)
    m.def(
        "ISM", [](Real n_ism) -> MediumVariant { return PyISM(n_ism); }, py::arg("n_ism"));

    m.def(
        "Wind",
        [](Real A_star, std::optional<Real> n_ism, std::optional<Real> n0, Real k_m) -> MediumVariant {
            return PyWind(A_star, n_ism, n0, k_m);
        },
        py::arg("A_star"), py::arg("n_ism") = py::none(), py::arg("n0") = py::none(), py::arg("k_m") = 2);

    // Observer bindings
    py::class_<PyObserver>(m, "Observer")
        .def(py::init<Real, Real, Real, Real>(), py::arg("lumi_dist"), py::arg("z"), py::arg("theta_obs"),
             py::arg("phi_obs") = 0)
        .def_property_readonly("lumi_dist", &PyObserver::lumi_dist_cgs)
        .def_readonly("z", &PyObserver::z)
        .def_readonly("theta_obs", &PyObserver::theta_obs)
        .def_readonly("phi_obs", &PyObserver::phi_obs)
        .def("__repr__", &PyObserver::repr);

    // Radiation bindings
    py::class_<PyRadiation>(m, "Radiation")
        .def(py::init<Real, Real, Real, Real, bool, bool>(), py::arg("eps_e"), py::arg("eps_B"), py::arg("p"),
             py::arg("xi_e") = 1, py::arg("ssc") = false, py::arg("kn") = false)
        .def_property_readonly("eps_e", [](PyRadiation const& r) { return r.rad.eps_e; })
        .def_property_readonly("eps_B", [](PyRadiation const& r) { return r.rad.eps_B; })
        .def_property_readonly("p", [](PyRadiation const& r) { return r.rad.p; })
        .def_property_readonly("xi_e", [](PyRadiation const& r) { return r.rad.xi_e; })
        .def_readonly("ssc", &PyRadiation::ssc)
        .def_readonly("kn", &PyRadiation::kn)
        .def("__repr__", &PyRadiation::repr);

    // Model bindings
    py::class_<PyModel>(m, "Model")
        // NOTE: pybind11's std::variant caster cannot LOAD JetVariant/MediumVariant (their
        // alternatives hold const members, so the variant is not copy-assignable and the
        // caster's load path does not compile) — hence the explicit isinstance dispatch.
        .def(py::init([](py::object jet_obj, py::object medium_obj, PyObserver observer, PyRadiation fwd_rad,
                         std::optional<PyRadiation> rvs_rad, std::optional<std::tuple<Real, Real, Real>> resolutions,
                         Real rtol, bool axisymmetric, bool radiative_fireball) -> PyModel {
                 auto jet = [&]() -> JetVariant {
                     if (py::isinstance<TophatJet>(jet_obj)) {
                         return jet_obj.cast<TophatJet>();
                     }
                     if (py::isinstance<GaussianJet>(jet_obj)) {
                         return jet_obj.cast<GaussianJet>();
                     }
                     if (py::isinstance<PowerLawJet>(jet_obj)) {
                         return jet_obj.cast<PowerLawJet>();
                     }
                     if (py::isinstance<Ejecta>(jet_obj)) {
                         return jet_obj.cast<Ejecta>();
                     }
                     throw py::type_error("jet must be TophatJet, GaussianJet, PowerLawJet, or Ejecta");
                 }();

                 auto medium = [&]() -> MediumVariant {
                     if (py::isinstance<ISM>(medium_obj)) {
                         return medium_obj.cast<ISM>();
                     }
                     if (py::isinstance<Wind>(medium_obj)) {
                         return medium_obj.cast<Wind>();
                     }
                     if (py::isinstance<Medium>(medium_obj)) {
                         return medium_obj.cast<Medium>();
                     }
                     throw py::type_error("medium must be ISM, Wind, or Medium");
                 }();

                 return PyModel(std::move(jet), std::move(medium), observer, fwd_rad, rvs_rad, resolutions, rtol,
                                axisymmetric, radiative_fireball);
             }),
             py::arg("jet"), py::arg("medium"), py::arg("observer"), py::arg("fwd_rad"),
             py::arg("rvs_rad") = py::none(), py::arg("resolutions") = py::none(),
             py::arg("rtol") = defaults::solver::dynamics_rtol, py::arg("axisymmetric") = true,
             py::arg("radiative_fireball") = true)

        .def("flux_density_grid", &PyModel::flux_density_grid, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux_density", &PyModel::flux_density, py::arg("t"), py::arg("nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux", &PyModel::flux, py::arg("t"), py::arg("nu_min"), py::arg("nu_max"), py::arg("num_nu"),
             py::call_guard<py::gil_scoped_release>())

        .def("flux_density_exposures", &PyModel::flux_density_exposures, py::arg("t"), py::arg("nu"),
             py::arg("expo_time"), py::arg("num_points") = 10, py::call_guard<py::gil_scoped_release>())

        .def("sky_image", &PyModel::sky_image, py::arg("t_obs"), py::arg("nu_obs"), py::arg("fov"),
             py::arg("npixel") = 256, py::call_guard<py::gil_scoped_release>())

        .def("details", &PyModel::details, py::arg("t_min"), py::arg("t_max"), py::call_guard<py::gil_scoped_release>())

        .def("medium", &PyModel::medium, py::arg("phi"), py::arg("theta"), py::arg("r"),
             py::call_guard<py::gil_scoped_release>())

        .def("jet_E_iso", &PyModel::jet_E_iso, py::arg("phi"), py::arg("theta"),
             py::call_guard<py::gil_scoped_release>())

        .def("jet_Gamma0", &PyModel::jet_Gamma0, py::arg("phi"), py::arg("theta"),
             py::call_guard<py::gil_scoped_release>())

        .def_property_readonly("observer", &PyModel::get_observer)
        .def_property_readonly("fwd_rad", &PyModel::get_fwd_rad)
        .def_property_readonly("rvs_rad", &PyModel::get_rvs_rad)
        .def_property_readonly("resolutions", &PyModel::get_resolutions)
        .def_property_readonly("radiative_fireball", &PyModel::get_radiative_fireball)
        .def_property_readonly("rtol", &PyModel::get_rtol)
        .def_property_readonly("axisymmetric", &PyModel::get_axisymmetric)
        .def("__repr__", &PyModel::repr)
#ifdef AFTERGLOW_PROFILE
        .def_static("profile_data", &PyModel::profile_data, "Get per-stage timing from the last computation (ms)")
        .def_static("profile_counters", &PyModel::profile_counters,
                    "Get profiling counters (rays, ODE steps, RHS evals)")
        .def_static("profile_reset", &PyModel::profile_reset, "Reset profiling counters")
#endif
        ;

    py::class_<PySkyImage>(m, "SkyImage")
        .def_readonly("image", &PySkyImage::image)
        .def_property_readonly("extent", [](PySkyImage const& s) { return py::array_t<Real>(4, s.extent.data()); })
        .def_readonly("pixel_solid_angle", &PySkyImage::pixel_solid_angle)
        .def("__repr__", &PySkyImage::repr);

    py::class_<Flux>(m, "Flux")
        .def(py::init<>())
        .def_readonly("sync", &Flux::sync)
        .def_readonly("ssc", &Flux::ssc)
        .def("__repr__", &Flux::repr);

    py::class_<PyFlux>(m, "FluxDict")
        .def(py::init<>())
        .def_readonly("fwd", &PyFlux::fwd)
        .def_readonly("rvs", &PyFlux::rvs)
        .def_readonly("total", &PyFlux::total)
        .def("__repr__", &PyFlux::repr);

    // Spectrum evaluators and grid accessors for details() callable interface
    py::class_<SpectrumEvaluator>(m, "_SpectrumEvaluator")
        .def("__call__", &SpectrumEvaluator::operator(), py::arg("nu_comv"));

    py::class_<YEvaluator>(m, "_YEvaluator").def("__call__", &YEvaluator::operator(), py::arg("gamma"));

    py::class_<SynSpectrumGrid>(m, "_SynSpectrumGrid")
        .def(
            "__getitem__",
            [](SynSpectrumGrid const& g, py::tuple idx) {
                auto& ph = (*g.grid_)(idx[0].cast<size_t>(), idx[1].cast<size_t>(), idx[2].cast<size_t>());
                const Real I_unit = unit::flux_den_cgs;
                return SpectrumEvaluator{
                    [&ph, I_unit](Real nu_comv) { return ph.compute_I_nu(nu_comv * unit::Hz) / I_unit; }};
            },
            py::return_value_policy::reference_internal);

    py::class_<ICSpectrumGrid>(m, "_ICSpectrumGrid")
        .def(
            "__getitem__",
            [](ICSpectrumGrid& g, py::tuple idx) {
                auto& ph = (*g.grid_)(idx[0].cast<size_t>(), idx[1].cast<size_t>(), idx[2].cast<size_t>());
                const Real I_unit = unit::flux_den_cgs;
                return SpectrumEvaluator{
                    [&ph, I_unit](Real nu_comv) { return ph.compute_I_nu(nu_comv * unit::Hz) / I_unit; }};
            },
            py::return_value_policy::reference_internal);

    py::class_<YSpectrumGrid>(m, "_YSpectrumGrid")
        .def(
            "__getitem__",
            [](YSpectrumGrid const& g, py::tuple idx) {
                auto& ph = (*g.grid_)(idx[0].cast<size_t>(), idx[1].cast<size_t>(), idx[2].cast<size_t>());
                return YEvaluator{[&ph](Real gamma) { return ph.Ys.gamma_spectrum(gamma); }};
            },
            py::return_value_policy::reference_internal);

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
        .def_readonly("Doppler", &PyShock::Doppler)
        .def_property_readonly(
            "sync_spectrum",
            [](PyShock& self) -> py::object {
                if (self.syn_photons_.size() == 0) {
                    return py::none();
                }
                return py::cast(SynSpectrumGrid{&self.syn_photons_});
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly(
            "ssc_spectrum",
            [](PyShock& self) -> py::object {
                if (self.ic_photons_.size() == 0) {
                    return py::none();
                }
                return py::cast(ICSpectrumGrid{&self.ic_photons_});
            },
            py::return_value_policy::reference_internal)
        .def_property_readonly(
            "Y_spectrum",
            [](PyShock& self) -> py::object {
                if (self.syn_photons_.size() == 0) {
                    return py::none();
                }
                return py::cast(YSpectrumGrid{&self.syn_photons_});
            },
            py::return_value_policy::reference_internal)
        .def("__repr__", &PyShock::repr);

    py::class_<PyDetails>(m, "SimulationDetails")
        .def(py::init<>())
        .def_readonly("phi", &PyDetails::phi)
        .def_readonly("theta", &PyDetails::theta)
        .def_readonly("t_src", &PyDetails::t_src)
        .def_readonly("fwd", &PyDetails::fwd)
        .def_readonly("rvs", &PyDetails::rvs)
        .def("__repr__", &PyDetails::repr);

    //========================================================================================================
    //                                 Utilities
    //========================================================================================================
    // Standalone utility: logscale_screen
    m.def("logscale_screen", &logscale_screen, py::arg("data"), py::arg("data_density"),
          "Select indices that subsample an array to ~data_density points per decade in log-space.");
}
