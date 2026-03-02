.. _advanced-fitting:

Advanced Fitting
================

VegasAfterglow's ``Fitter`` is built directly on the ``Model`` class, giving you full control over the fitting process. You can customize:

- **Priors**: Use informative or custom prior distributions
- **Likelihood**: Implement non-standard likelihood functions (e.g., upper limits, systematic uncertainties)
- **Jet profiles**: Define arbitrary angular energy and Lorentz factor distributions
- **Medium profiles**: Define custom circumburst density structures

The ``Fitter`` evaluates models using Python's ``ThreadPoolExecutor``, where each thread constructs a ``Model`` and calls ``flux_density()`` with the GIL released during C++ computation. This provides near-native parallelism with full Python flexibility.

.. _custom-priors:

Custom Priors
--------------

Pass a ``priors`` dictionary to ``fit()`` to apply custom prior distributions. The same interface works for all samplers (emcee, dynesty, etc.).

Keys are parameter labels (using the ``log10_`` prefix for ``LOG``-scale parameters), and values are ``bilby.core.prior.Prior`` objects. Parameters not in the dict automatically get uniform priors based on the ``ParamDef`` bounds.

.. code-block:: python

    import bilby

    custom_priors = {
        "p": bilby.core.prior.Gaussian(
            mu=2.3, sigma=0.1,
            minimum=2.1, maximum=2.8,
            name="p", latex_label=r"$p$",
        ),
        "log10_eps_e": bilby.core.prior.Gaussian(
            mu=-1.0, sigma=0.5,
            minimum=-3, maximum=np.log10(0.5),
            name="log10_eps_e",
            latex_label=r"$\log_{10}(\epsilon_e)$",
        ),
    }

    # Works with emcee
    result = fitter.fit(
        params,
        sampler="emcee",
        nsteps=10000,
        priors=custom_priors,
    )

    # Same priors work with dynesty
    result = fitter.fit(
        params,
        sampler="dynesty",
        nlive=1000,
        priors=custom_priors,
    )

.. important::
    When using ``LOG``-scale parameters, the prior keys must use the ``log10_`` prefix (e.g., ``log10_eps_e``, not ``eps_e``), since the sampler operates in log-space.

.. _custom-likelihood:

Custom Likelihood
------------------

The default likelihood is Gaussian: :math:`\ln\mathcal{L} = -\chi^2/2`. Override this with a ``log_likelihood_fn`` that maps :math:`\chi^2` to log-likelihood:

.. code-block:: python

    def student_t_likelihood(chi2):
        """Student-t likelihood with nu=5 degrees of freedom.

        More robust to outliers than Gaussian.
        """
        nu = 5
        n_data = 50  # number of data points
        return -0.5 * (nu + n_data) * np.log(1 + chi2 / nu)

    result = fitter.fit(
        params,
        sampler="emcee",
        log_likelihood_fn=student_t_likelihood,
    )

**Upper Limits**

For non-detections, a common approach is to weight the likelihood so data points
above the upper limit are penalized:

.. code-block:: python

    def likelihood_with_upper_limits(chi2):
        """Standard Gaussian likelihood; upper limits encoded in data weights."""
        return -0.5 * chi2

    # Upper limits: set flux to 0 and error to the upper limit value,
    # so chi2 = (0 - model)^2 / upper_limit^2, which penalizes models above the limit
    fitter.add_flux_density(
        nu=1e10, t=[1e5], f_nu=[0.0], err=[3e-29],
        weights=[1.0],
    )

Custom Jet Profiles
--------------------

The ``Fitter`` accepts a custom jet factory function via the ``jet`` keyword argument.
This lets you define arbitrary angular profiles using the same ``Ejecta`` class used for
direct model calculations (see :doc:`/examples/index`).

**How it works:**

1. Define functions for the energy and Lorentz factor angular profiles
2. Wrap them in a *factory function* that takes ``params`` (the current MCMC sample) and returns an ``Ejecta``
3. Pass the factory to ``Fitter(jet=...)``

The fitter calls your factory on every MCMC step with fresh parameter values.

.. code-block:: python

    from VegasAfterglow import Ejecta, Fitter, ParamDef, Scale

    def double_gaussian_jet(params):
        """Double-Gaussian jet: narrow core + wide wing."""

        def E_iso_func(phi, theta):
            core = params.E_iso * np.exp(-(theta / params.theta_c) ** 2)
            wing = params.E_iso_w * np.exp(-(theta / params.theta_w) ** 2)
            return core + wing

        def Gamma_func(phi, theta):
            return 1 + (params.Gamma0 - 1) * np.exp(-(theta / params.theta_c) ** 2)

        return Ejecta(E_iso=E_iso_func, Gamma0=Gamma_func, duration=params.tau)

    mc_params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    50,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.15, Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("E_iso_w", 1e48,  1e52,  Scale.log),
        ParamDef("theta_w",  0.1,   0.5,  Scale.linear),
        ParamDef("tau",        1,   1e3,  Scale.log),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet=double_gaussian_jet)
    fitter.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)

    result = fitter.fit(mc_params, sampler="emcee", nsteps=10000)

The ``Ejecta`` constructor supports these optional keyword arguments:

- ``sigma0``: Magnetization profile ``sigma0(phi, theta)`` (default: 0)
- ``E_dot``: Energy injection rate ``E_dot(phi, theta, t)`` in erg/s (default: 0)
- ``M_dot``: Mass injection rate ``M_dot(phi, theta, t)`` in g/s (default: 0)
- ``spreading``: Enable lateral spreading (default: False)
- ``duration``: Jet duration in seconds, for reverse shock (default: 1)

.. note::
    When using a custom ``jet``, parameter validation is automatically skipped since the fitter cannot know which parameters your factory requires. Ensure your factory function handles all necessary parameters.

**Energy and Mass Injection**

The ``Ejecta`` class also supports time-dependent energy and mass injection via ``E_dot`` and ``M_dot``.
These are functions of ``(phi, theta, t)`` returning injection rates in erg/s and g/s respectively:

.. code-block:: python

    def jet_with_injection(params):
        """Custom jet with energy injection (spin-down luminosity)."""

        def E_iso_func(phi, theta):
            return params.E_iso * np.exp(-(theta / params.theta_c) ** 2)

        def Gamma_func(phi, theta):
            return 1 + (params.Gamma0 - 1) * np.exp(-(theta / params.theta_c) ** 2)

        def E_dot_func(phi, theta, t):
            """Spin-down energy injection: L(t) = L0 / (1 + t/t_sd)^2."""
            return params.L0 / (1 + t / params.t_sd) ** 2

        def M_dot_func(phi, theta, t):
            """Mass injection: decaying mass loading."""
            return params.M_dot0 * np.exp(-t / params.t_sd)

        return Ejecta(
            E_iso=E_iso_func, Gamma0=Gamma_func,
            E_dot=E_dot_func, M_dot=M_dot_func,
            duration=params.tau,
        )

    mc_params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    50,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.15, Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("L0",      1e44,  1e48,  Scale.log),       # Injection luminosity [erg/s]
        ParamDef("t_sd",      10,  1e4,   Scale.log),       # Spin-down timescale [s]
        ParamDef("M_dot0",  1e20,  1e26,  Scale.log),       # Mass injection rate [g/s]
        ParamDef("tau",        1,   1e3,  Scale.log),
        # ... other standard params (n_ism, p, eps_e, eps_B, xi_e)
    ]

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet=jet_with_injection)

See :doc:`/examples/models` for more user-defined jet examples.

Custom Medium Profiles
-----------------------

Similarly, pass a custom medium factory via the ``medium`` keyword argument to define custom density profiles using the ``Medium`` class:

.. code-block:: python

    from VegasAfterglow import Medium, Fitter, ParamDef, Scale

    def exponential_medium(params):
        """Exponential density profile: rho = mp * n_ism * exp(-r / r_scale)."""
        mp = 1.67e-24  # proton mass [g]

        def density_func(phi, theta, r):
            return mp * params.n_ism * np.exp(-r / params.r_scale)

        return Medium(rho=density_func)

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, medium=exponential_medium)

The ``Medium`` class takes a callable ``rho(phi, theta, r)`` that returns the mass density [g/cm³] at position (phi, theta, r). See :doc:`/examples/models` for more user-defined medium examples.

**Custom MCMC Parameters**

When using custom jet or medium functions, you can define arbitrary MCMC parameters beyond the standard set. Simply include them in the ``ParamDef`` list -- the fitter automatically uses a Python-based parameter transformer when it detects non-standard parameter names:

.. code-block:: python

    mc_params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    50,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.15, Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("r_scale", 1e16,  1e20,  Scale.log),     # custom parameter
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, medium=exponential_medium)
    fitter.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)
    result = fitter.fit(mc_params, sampler="emcee", nsteps=10000)

The custom parameter ``r_scale`` is accessible inside the factory via ``params.r_scale``, and it is sampled by the MCMC just like any standard parameter. Standard parameters (``theta_v``, ``eps_e``, etc.) retain their defaults and are still used internally for observer geometry and radiation physics.

.. note::
    Standard ``ModelParams`` fields (like ``theta_v``, ``eps_e``, ``eps_B``, ``p``, ``xi_e``) must still be included in the ``ParamDef`` list -- they are needed by the ``Observer`` and ``Radiation`` objects that the fitter constructs internally.

**Combining Custom Jet and Medium**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet=double_gaussian_jet,
        medium=exponential_medium,
    )
    fitter.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)
    fitter.add_flux_density(nu=2.4e17, t=t_xray, f_nu=flux_xray, err=err_xray)

    result = fitter.fit(
        mc_params,
        sampler="emcee",
        nsteps=20000,
        nburn=5000,
        log_likelihood_fn=student_t_likelihood,
    )


Speeding Up Custom Profiles with ``@gil_free``
------------------------------------------------

The custom jet and medium examples above use plain Python callbacks. Each time C++
evaluates your profile function, it must cross the Python↔C++ boundary — acquiring the
GIL, calling into the interpreter, and returning. During blast wave evolution, this
happens hundreds of times per model evaluation, adding overhead even in single-threaded
mode and preventing true parallelism in multi-threaded MCMC.

The ``@gil_free`` decorator compiles your profile function to native machine code
using `numba <https://numba.pydata.org/>`_, so C++ calls it directly as a C function
pointer — no Python interpreter, no GIL, no boundary crossing:

.. code-block:: bash

    pip install numba

The two differences from plain Python callbacks are:

1. Add the ``@gil_free`` decorator
2. Pass MCMC parameters as explicit function arguments (after the spatial coordinates)
   instead of capturing them from the enclosing scope — you can add as many parameters
   as you need

Here is a complete example — a tophat jet with ``@gil_free``:

.. code-block:: python

    import math
    from VegasAfterglow import Fitter, Ejecta, ParamDef, Scale, gil_free

    # Decorate profile functions with @gil_free.
    # First 2 args (phi, theta) are spatial coordinates from C++.
    # Additional args are MCMC parameters, bound by keyword below.
    @gil_free
    def tophat_energy(phi, theta, E_iso, theta_c):
        return E_iso if theta <= theta_c else 0.0

    @gil_free
    def tophat_gamma(phi, theta, Gamma0, theta_c):
        return Gamma0 if theta <= theta_c else 1.0

    # Factory: bind MCMC parameters by keyword each step
    def jet_factory(mc_params):
        return Ejecta(
            E_iso=tophat_energy(E_iso=mc_params.E_iso, theta_c=mc_params.theta_c),
            Gamma0=tophat_gamma(Gamma0=mc_params.Gamma0, theta_c=mc_params.theta_c),
        )

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet=jet_factory, medium="wind")
    fitter.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)

    mc_params = [
        ParamDef("E_iso",    1e51,  1e54,  Scale.log,    1e53),
        ParamDef("Gamma0",      5,  1000,  Scale.log,      20),
        ParamDef("theta_c",  0.01,   0.5,  Scale.linear,  0.1),
        ParamDef("p",           2,     3,  Scale.linear,  2.3),
        ParamDef("eps_e",    1e-2,   0.3,  Scale.log,    0.05),
        ParamDef("eps_B",    1e-4,   0.3,  Scale.log,    0.03),
        ParamDef("xi_e",     1e-3,   0.1,  Scale.log,    0.01),
        ParamDef("A_star",   1e-3,    10,  Scale.log,    0.05),
    ]

    result = fitter.fit(mc_params, sampler="emcee", nsteps=10000)

The same decorator works for custom medium density functions (3 spatial coordinates
``phi, theta, r``) and for energy/mass injection functions (``phi, theta, t``):

.. code-block:: python

    @gil_free
    def custom_wind(phi, theta, r, A_star):
        return A_star * 5e11 * 1.67e-24 / (r * r)

    def medium_factory(mc_params):
        return Medium(rho=custom_wind(A_star=mc_params.A_star))

    @gil_free
    def spindown_injection(phi, theta, t, L0, t_sd):
        return L0 / (1.0 + t / t_sd) ** 2

.. tip::
    Functions decorated with ``@gil_free`` must use the ``math`` module
    (not ``numpy``) and only simple arithmetic — no Python objects, arrays, or closures.
    If you need more complex logic, use the plain Python callback approach instead.

.. note::
    Built-in jet types (``"tophat"``, ``"gaussian"``, ``"powerlaw"``, etc.) are already
    implemented in C++ and do not need this decorator.


Performance Notes
------------------

**Threading vs. Multiprocessing**

The ``Fitter`` uses ``ThreadPoolExecutor`` (threads, not processes) because:

1. **GIL is released** during the C++ ``flux_density()`` computation, so threads run truly in parallel
2. **No pickling overhead**: Model objects stay in the same process
3. **Shared memory**: All threads share observation data without copying

For typical afterglow models, the Python overhead (object creation, GIL acquisition) is <5% of total time for synchrotron-only models and <0.1% for SSC models.

**Parallelism Tips**

- **Emcee**: ``npool`` defaults to the number of CPU cores.
- **Dynesty**: ``npool`` controls thread count. The ``queue_size`` is automatically optimized.
- **Custom jet/medium with ``@gil_free``**: Full thread parallelism is preserved. All profile evaluations run as native C function calls without the GIL.
- **Custom jet/medium with plain Python callbacks**: Python callbacks require the GIL, which serializes the angular profile evaluation across threads. The blast wave evolution and radiation computation still run without the GIL, but the profile evaluation becomes a bottleneck. Use ``@gil_free`` to eliminate this overhead.

Troubleshooting
================

For comprehensive troubleshooting help including MCMC convergence issues, data selection problems, memory optimization, and performance tuning, see :doc:`/troubleshooting`.

.. seealso::
   :doc:`/validation` for code validation and comparison with other afterglow codes. :doc:`/parameter_reference` for the complete parameter reference.
