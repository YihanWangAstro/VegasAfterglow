Model Configurations
====================

.. seealso::
   :doc:`/parameter_reference` for all available parameters. :doc:`/examples/index` for direct model calculation examples.

TopHat + ISM (Default)
-----------------------

The default configuration uses a top-hat jet in a uniform ISM environment with forward shock synchrotron emission:

.. code-block:: python

    # Basic configuration
    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")

    # Basic parameter set
    params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),     # Isotropic energy in erg
        ParamDef("Gamma0",    10,   500,  Scale.log),     # Lorentz factor
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),  # Opening angle in radians
        ParamDef("theta_v",    0,     0,  Scale.fixed),   # Viewing angle (on-axis) in radians
        ParamDef("n_ism",   1e-3,   100,  Scale.log),     # Number density in cm^-3
        ParamDef("p",        2.1,   2.8,  Scale.linear),  # Electron spectral index
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),     # Electron energy fraction
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),     # Magnetic energy fraction
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),  # Fraction of accelerated electrons
    ]

Jet Structure Variations
-------------------------

**Power-law Structured Jet**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="powerlaw", medium="ism")

    params = [
        # Basic jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.3,  Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),  # Allow off-axis viewing

        # Power-law structure parameters
        ParamDef("k_e",      1.5,   3.0,  Scale.linear),  # Energy power-law index, default 2.0 if not specified
        ParamDef("k_g",      1.5,   3.0,  Scale.linear),  # Lorentz factor power-law, default 2.0 if not specified

        # Medium and microphysics (same as default)
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Gaussian Structured Jet**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="gaussian", medium="ism")

    params = [
        # Basic parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.2,  Scale.linear),  # Gaussian width parameter
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Two-Component Jet**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="two_component", medium="ism")

    params = [
        # Narrow component
        ParamDef("E_iso",   1e50,  1e53,  Scale.log),     # Core energy
        ParamDef("Gamma0",   100,   500,  Scale.log),     # Core Lorentz factor
        ParamDef("theta_c", 0.01,   0.1,  Scale.linear),  # Core angle

        # Wide component
        ParamDef("E_iso_w", 1e49,  1e52,  Scale.log),     # Wide energy in erg
        ParamDef("Gamma0_w",  10,   100,  Scale.log),     # Wide Lorentz factor
        ParamDef("theta_w",  0.1,   0.5,  Scale.linear),  # Wide angle in radians

        # Observation and medium (same as default)
        ParamDef("theta_v",    0,   0.3,  Scale.linear),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Step Power-law Jet**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="step_powerlaw", medium="ism")

    params = [
        # Core component (uniform)
        ParamDef("E_iso",   1e51,  1e54,  Scale.log),     # Core energy
        ParamDef("Gamma0",    50,   500,  Scale.log),     # Core Lorentz factor
        ParamDef("theta_c", 0.01,   0.1,  Scale.linear),  # Core boundary

        # Wing component (power-law)
        ParamDef("E_iso_w", 1e49,  1e52,  Scale.log),     # Wing energy scale
        ParamDef("Gamma0_w",  10,   100,  Scale.log),     # Wing Lorentz factor
        ParamDef("k_e",      1.5,   3.0,  Scale.linear),  # Energy power-law
        ParamDef("k_g",      1.5,   3.0,  Scale.linear),  # Lorentz factor power-law

        # Standard parameters (same as default)
        ParamDef("theta_v",    0,   0.3,  Scale.linear),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

Medium Type Variations
-----------------------

**Stellar Wind Medium**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),

        # Wind medium parameter (replaces n_ism)
        ParamDef("A_star",  1e-3,   1.0,  Scale.log),     # Wind parameter

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Stratified Medium: ISM-to-Wind**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),

        # Stratified medium parameters
        ParamDef("A_star",  1e-5,   0.1,  Scale.log),     # Wind strength (outer)
        ParamDef("n0",      1e-3,    10,  Scale.log),     # ISM density (inner) in cm^-3

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Stratified Medium: Wind-to-ISM**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),

        # Stratified medium (wind -> ISM)
        ParamDef("A_star",  1e-3,   1.0,  Scale.log),     # Inner wind strength
        ParamDef("n_ism",   1e-3,   100,  Scale.log),     # Outer ISM density

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Stratified Medium: ISM-Wind-ISM**

.. code-block:: python

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),

        # Three-zone stratified medium
        ParamDef("A_star",  1e-4,   0.1,  Scale.log),     # Wind parameter (middle)
        ParamDef("n_ism",   1e-3,   100,  Scale.log),     # Outer ISM density
        ParamDef("n0",      1e-2,    20,  Scale.log),     # Inner ISM density

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

.. important::
    **Stratified Medium Physics:**

    - **A_star = 0**: Pure ISM with density n_ism
    - **n0 = infinity**: Pure wind profile from center
    - **A_star > 0, n0 < infinity**: ISM-wind-ISM stratification
    - **A_star > 0, n0 = infinity**: Wind-ISM stratification

    **Density Profile:** Inner (r < r1): n = n0, Middle (r1 < r < r2): n proportional to A_star/r^2, Outer (r > r2): n = n_ism

Reverse Shock
--------------

**Basic Reverse Shock**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        rvs_shock=True,
    )

    params = [
        # Standard jet and medium parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),

        # Jet duration (important for reverse shock)
        ParamDef("tau",        1,   1e6,  Scale.log),     # Jet duration in seconds

        # Forward shock microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),

        # Reverse shock microphysics (can be different)
        ParamDef("p_r",      2.1,   2.8,  Scale.linear),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.log),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.log),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.linear),
    ]

**Reverse Shock with Structured Jet**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="gaussian", medium="ism",
        rvs_shock=True,
    )

    params = [
        # Gaussian jet parameters
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    50,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.2,  Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("tau",        1,   1e6,  Scale.log),

        # Forward + reverse shock microphysics
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
        ParamDef("p_r",      2.1,   2.8,  Scale.linear),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.log),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.log),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.linear),
    ]

Inverse Compton Radiation
--------------------------

For the physical details of inverse Compton and Klein-Nishina corrections, see :doc:`/physics`.

**Forward Shock Inverse Compton**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        fwd_ssc=True,
        kn=True,
    )

    params = [
        # Standard parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

**Reverse Shock Inverse Compton**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        rvs_shock=True,
        fwd_ssc=True,
        rvs_ssc=True,
        kn=True,
    )

    params = [
        # Standard parameters with reverse shock
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("tau",        1,   100,  Scale.log),

        # Forward + reverse microphysics
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
        ParamDef("p_r",      2.1,   2.8,  Scale.linear),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.log),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.log),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.linear),
    ]

Energy Injection
-----------------

**Magnetar Spin-down Injection**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        magnetar=True,
    )

    params = [
        # Standard jet and medium parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.5,  Scale.linear),
        ParamDef("theta_v",    0,     0,  Scale.fixed),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),

        # Magnetar injection parameters
        ParamDef("L0",      1e42,  1e48,  Scale.log),     # Initial luminosity [erg/s]
        ParamDef("t0",        10,  1000,  Scale.log),     # Spin-down timescale [s]
        ParamDef("q",        1.5,   3.0,  Scale.linear),  # Power-law index

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

.. note::
    **Magnetar Injection Profile:** L(t) = L0 x (1 + t/t0)^(-q) for theta < theta_c


**Magnetar with Structured Jet**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="powerlaw", medium="ism",
        magnetar=True,
    )

    params = [
        # Power-law jet with magnetar
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    10,   500,  Scale.log),
        ParamDef("theta_c", 0.01,   0.3,  Scale.linear),
        ParamDef("k_e",      1.5,   3.0,  Scale.linear),
        ParamDef("k_g",      1.5,   3.0,  Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),

        # Magnetar parameters
        ParamDef("L0",      1e42,  1e48,  Scale.log),
        ParamDef("t0",        10,  1000,  Scale.log),
        ParamDef("q",        1.5,   3.0,  Scale.linear),

        # Standard microphysics
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),
    ]

Complex Model Combinations
----------------------------

**Full Physics: Structured Jet + Stratified Medium + Reverse Shock + IC + Magnetar**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="gaussian", medium="wind",
        rvs_shock=True,
        fwd_ssc=True,
        rvs_ssc=True,
        kn=True,
        magnetar=True,
    )

    params = [
        # Gaussian jet
        ParamDef("E_iso",   1e50,  1e54,  Scale.log),
        ParamDef("Gamma0",    50,   500,  Scale.log),
        ParamDef("theta_c", 0.02,   0.2,  Scale.linear),
        ParamDef("theta_v",    0,   0.5,  Scale.linear),
        ParamDef("tau",        1,   100,  Scale.log),

        # Stratified medium
        ParamDef("A_star",  1e-4,   1.0,  Scale.log),
        ParamDef("n_ism",   1e-3,   100,  Scale.log),
        ParamDef("n0",      1e-2,    50,  Scale.log),

        # Magnetar injection
        ParamDef("L0",      1e42,  1e48,  Scale.log),
        ParamDef("t0",        10,  1000,  Scale.log),
        ParamDef("q",        1.5,   3.0,  Scale.linear),

        # Forward shock microphysics
        ParamDef("p",        2.1,   2.8,  Scale.linear),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.log),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.log),
        ParamDef("xi_e",     0.1,   1.0,  Scale.linear),

        # Reverse shock microphysics
        ParamDef("p_r",      2.1,   2.8,  Scale.linear),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.log),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.log),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.linear),
    ]

.. warning::
    **Complex Model Considerations:**
    - Use coarser resolution initially: ``resolution=(0.2, 0.7, 7)``
    - Increase live points: ``nlive=1000+``
    - Lower dlogz for better convergence: ``dlogz=0.05``
    - Consider parameter degeneracies in interpretation
