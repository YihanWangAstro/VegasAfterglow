MCMC Parameter Fitting
======================

.. contents:: Table of Contents
   :local:
   :depth: 3

Overview
--------

VegasAfterglow provides a comprehensive MCMC framework for parameter estimation of gamma-ray burst (GRB) afterglow models. The framework supports all built-in jet models, ambient medium configurations, radiation processes, and complex multi-wavelength datasets.

Key Features:

- **Multi-threaded MCMC engine** for fast parameter estimation
- **All model types supported**: TophatJet, GaussianJet, PowerLawJet, TwoComponentJet, StepPowerLawJet
- **Complete physics**: Forward shock, reverse shock, synchrotron, inverse Compton, magnetar injection
- **Flexible data handling**: Light curves and spectra with optional weighting

Basic Setup
-----------

Setting up Observational Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ObsData`` class handles multi-wavelength observational data:

.. code-block:: python

    import numpy as np
    import pandas as pd
    from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale

    # Create data container
    data = ObsData()

    # Method 1: Add data directly
    t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
    flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # erg/cm²/s/Hz
    flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]   # Error bars

    # Add light curve at R-band frequency
    data.add_flux_density(nu=4.84e14, t=t_data,
                         f_nu=flux_data, err=flux_err)  # All quantities in CGS units

    # Optional: Add weights for systematic uncertainties, normalization handled internally
    weights = np.ones(len(t_data))  # Equal weights
    data.add_flux_density(nu=2.4e17, t=t_data,
                         f_nu=flux_data, err=flux_err,
                         weights=weights)  # All quantities in CGS units

    # Method 2: Add frequency-integrated light curve (broadband flux)
    # For instruments with wide frequency coverage (e.g., BAT, LAT, Fermi)
    nu_min = 1e17  # Lower frequency bound [Hz]
    nu_max = 1e19  # Upper frequency bound [Hz]
    num_points = 5  # Number of frequency points for integration

    data.add_flux(nu_min=nu_min, nu_max=nu_max,
                         num_points=num_points, t=t_data,
                         flux=flux_data, err=flux_err,
                         weights=weights)  # All quantities in CGS units

    # Method 3: Load from CSV files
    bands = [2.4e17, 4.84e14, 1.4e14]  # X-ray, optical, near-IR
    lc_files = ["data/xray.csv", "data/optical.csv", "data/nir.csv"]

    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        data.add_flux_density(nu=nu, t=df["t"],
                             f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

    # Add spectra at specific times
    spec_times = [1000, 10000]  # seconds
    spec_files = ["data/spec_1000s.csv", "data/spec_10000s.csv"]

    for t, fname in zip(spec_times, spec_files):
        df = pd.read_csv(fname)
        data.add_spectrum(t=t, nu=df["nu"],
                          f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

Data Selection and Optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Smart Data Subsampling with logscale_screen**

For large datasets or densely sampled observations, using all available data points can lead to computational inefficiency and biased parameter estimation. The ``logscale_screen`` method provides intelligent data subsampling that maintains the essential information content while reducing computational overhead.

.. code-block:: python

    # Example: Large dense dataset
    t_dense = np.logspace(2, 7, 1000)  # 1000 time points
    flux_dense = np.random.lognormal(-60, 0.5, 1000)  # Dense flux measurements
    flux_err_dense = 0.1 * flux_dense

    # Subsample using logarithmic screening
    # This selects ~50-100 representative points across 5 decades in time
    indices = ObsData.logscale_screen(t_dense, num_order=5)

    # Add only the selected subset
    data.add_flux_density(nu=5e14,
                         t=t_dense[indices],
                         f_nu=flux_dense[indices],
                         err=flux_err_dense[indices])

**Why logscale_screen is Important:**

1. **Prevents Oversampling Bias**: Dense data clusters can dominate the χ² calculation, causing the MCMC to over-fit specific frequency bands or time periods.

2. **Computational Efficiency**: Reduces the number of model evaluations needed during MCMC sampling, significantly improving performance.

3. **Preserves Information**: Unlike uniform thinning, logarithmic sampling maintains representation across all temporal/spectral scales.

4. **Balanced Multi-band Fitting**: Ensures each frequency band contributes proportionally to the parameter constraints.

**Data Selection Guidelines:**

- **Target 10-30 points per frequency band** for balanced constraints
- **Avoid >100 points in any single band** unless scientifically justified
- **Maintain temporal coverage** across all evolutionary phases
- **Weight systematic uncertainties** appropriately using the weights parameter

.. warning::
    **Common Data Selection Pitfalls:**

    - **Optical-heavy datasets**: Dense optical coverage can bias parameters toward optical-dominant solutions
    - **Late-time clustering**: Too many late-time points can over-constrain decay slopes at the expense of early physics
    - **Single-epoch spectra**: Broadband spectra at one time can dominate multi-epoch light curves in χ² space

    **Solution**: Use ``logscale_screen`` for manual temporal reduction of over-sampled bands.

Global Configuration
^^^^^^^^^^^^^^^^^^^^

The ``Setups`` class defines fixed model properties:

.. code-block:: python

    cfg = Setups()

    # Source properties
    cfg.lumi_dist = 3.364e28  # Luminosity distance [cm]
    cfg.z = 1.58              # Redshift

    # Model selection (see sections below for all options)
    cfg.medium = "wind"       # Ambient medium type
    cfg.jet = "powerlaw"      # Jet structure type

    # Physics options
    cfg.rvs_shock = True      # Include reverse shock
    cfg.fwd_ssc = True        # Forward shock inverse Compton
    cfg.rvs_ssc = False       # Reverse shock inverse Compton
    cfg.ssc_cooling = True     # IC cooling effects
    cfg.kn = True             # Klein-Nishina corrections
    cfg.magnetar = True       # Magnetar energy injection

    # Numerical parameters
    cfg.rtol = 1e-5           # Numerical tolerance

Model Configurations
--------------------

Basic Setup (Default)
^^^^^^^^^^^^^^^^^^^^^

The default configuration uses a top-hat jet in a uniform ISM environment with forward shock synchrotron emission:

.. code-block:: python

    # Basic configuration
    cfg = Setups()
    cfg.medium = "ism"        # Uniform ISM density
    cfg.jet = "tophat"        # Top-hat jet structure

    # Basic parameter set
    params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),     # Isotropic energy in erg
        ParamDef("Gamma0",    10,   500,  Scale.LOG),     # Lorentz factor
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),  # Opening angle in radians
        ParamDef("theta_v",    0,     0,  Scale.FIXED),   # Viewing angle (on-axis) in radians
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),     # Number density in cm^-3
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),  # Electron spectral index
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),     # Electron energy fraction
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),     # Magnetic energy fraction
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),  # Fraction of accelerated electrons
    ]

Jet Structure Variations
^^^^^^^^^^^^^^^^^^^^^^^^

**Power-law Structured Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"        # Default ISM medium
    cfg.jet = "powerlaw"      # Power-law structured jet

    params = [
        # Basic jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.3,  Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),  # Allow off-axis viewing

        # Power-law structure parameters
        ParamDef("k_e",      1.5,   3.0,  Scale.LINEAR),  # Energy power-law index, default 2.0 if not specified
        ParamDef("k_g",      1.5,   3.0,  Scale.LINEAR),  # Lorentz factor power-law, default 2.0 if not specified

        # Medium and microphysics (same as default)
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Gaussian Structured Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "gaussian"      # Gaussian structured jet

    params = [
        # Basic parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.2,  Scale.LINEAR),  # Gaussian width parameter
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Two-Component Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "two_component"  # Two-component jet

    params = [
        # Narrow component
        ParamDef("E_iso",   1e50,  1e53,  Scale.LOG),     # Core energy
        ParamDef("Gamma0",   100,   500,  Scale.LOG),     # Core Lorentz factor
        ParamDef("theta_c", 0.01,   0.1,  Scale.LINEAR),  # Core angle

        # Wide component
        ParamDef("E_iso_w", 1e49,  1e52,  Scale.LOG),     # Wide energy in erg
        ParamDef("Gamma0_w",  10,   100,  Scale.LOG),     # Wide Lorentz factor
        ParamDef("theta_w",  0.1,   0.5,  Scale.LINEAR),  # Wide angle in radians

        # Observation and medium (same as default)
        ParamDef("theta_v",    0,   0.3,  Scale.LINEAR),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Step Power-law Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "step_powerlaw"  # Step power-law jet

    params = [
        # Core component (uniform)
        ParamDef("E_iso",   1e51,  1e54,  Scale.LOG),     # Core energy
        ParamDef("Gamma0",    50,   500,  Scale.LOG),     # Core Lorentz factor
        ParamDef("theta_c", 0.01,   0.1,  Scale.LINEAR),  # Core boundary

        # Wing component (power-law)
        ParamDef("E_iso_w", 1e49,  1e52,  Scale.LOG),     # Wing energy scale
        ParamDef("Gamma0_w",  10,   100,  Scale.LOG),     # Wing Lorentz factor
        ParamDef("k_e",      1.5,   3.0,  Scale.LINEAR),  # Energy power-law
        ParamDef("k_g",      1.5,   3.0,  Scale.LINEAR),  # Lorentz factor power-law

        # Standard parameters (same as default)
        ParamDef("theta_v",    0,   0.3,  Scale.LINEAR),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

Medium Type Variations
^^^^^^^^^^^^^^^^^^^^^^

**Stellar Wind Medium**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "wind"       # Stellar wind medium
    cfg.jet = "tophat"        # Default jet structure

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),

        # Wind medium parameter (replaces n_ism)
        ParamDef("A_star",  1e-3,   1.0,  Scale.LOG),     # Wind parameter

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Stratified Medium: ISM-to-Wind**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "wind"       # Use wind for stratified models
    cfg.jet = "tophat"        # Default jet structure

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),

        # Stratified medium parameters
        ParamDef("A_star",  1e-5,   0.1,  Scale.LOG),     # Wind strength (outer)
        ParamDef("n0",      1e-3,    10,  Scale.LOG),     # ISM density (inner) in cm^-3

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Stratified Medium: Wind-to-ISM**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "wind"
    cfg.jet = "tophat"

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),

        # Stratified medium (wind → ISM)
        ParamDef("A_star",  1e-3,   1.0,  Scale.LOG),     # Inner wind strength
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),     # Outer ISM density

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Stratified Medium: ISM-Wind-ISM**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "wind"
    cfg.jet = "tophat"

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),

        # Three-zone stratified medium
        ParamDef("A_star",  1e-4,   0.1,  Scale.LOG),     # Wind parameter (middle)
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),     # Outer ISM density
        ParamDef("n0",      1e-2,    20,  Scale.LOG),     # Inner ISM density

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

.. important::
    **Stratified Medium Physics:**

    - **A_star = 0**: Pure ISM with density n_ism
    - **n0 = ∞**: Pure wind profile from center
    - **A_star > 0, n0 < ∞**: ISM-wind-ISM stratification
    - **A_star > 0, n0 = ∞**: Wind-ISM stratification

    **Density Profile:** Inner (r < r₁): n = n0, Middle (r₁ < r < r₂): n ∝ A_star/r², Outer (r > r₂): n = n_ism

Reverse Shock
^^^^^^^^^^^^^

**Basic Reverse Shock**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"        # Default medium
    cfg.jet = "tophat"        # Default jet
    cfg.rvs_shock = True      # Enable reverse shock

    params = [
        # Standard jet and medium parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),

        # Jet duration (important for reverse shock)
        ParamDef("tau",        1,   1e6,  Scale.LOG),     # Jet duration in seconds

        # Forward shock microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),

        # Reverse shock microphysics (can be different)
        ParamDef("p_r",      2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.LINEAR),
    ]

**Reverse Shock with Structured Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "gaussian"      # Structured jet example
    cfg.rvs_shock = True

    params = [
        # Gaussian jet parameters
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    50,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.2,  Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("tau",        1,   1e6,  Scale.LOG),

        # Forward + reverse shock microphysics
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
        ParamDef("p_r",      2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.LINEAR),
    ]

Inverse Compton Radiation
^^^^^^^^^^^^^^^^^^^^^^^^^

**Forward Shock Inverse Compton**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"        # Default medium
    cfg.jet = "tophat"        # Default jet
    cfg.fwd_ssc = True        # Forward shock SSC
    cfg.ssc_cooling = True     # IC cooling effects
    cfg.kn = True             # Klein-Nishina corrections

    params = [
        # Standard parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

**Reverse Shock Inverse Compton**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "tophat"
    cfg.rvs_shock = True      # Enable reverse shock
    cfg.fwd_ssc = True        # Forward shock SSC
    cfg.rvs_ssc = True        # Reverse shock SSC
    cfg.ssc_cooling = True
    cfg.kn = True

    params = [
        # Standard parameters with reverse shock
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("tau",        1,   100,  Scale.LOG),

        # Forward + reverse microphysics
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
        ParamDef("p_r",      2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.LINEAR),
    ]

Energy Injection
^^^^^^^^^^^^^^^^

**Magnetar Spin-down Injection**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"        # Default medium
    cfg.jet = "tophat"        # Default jet
    cfg.magnetar = True       # Enable magnetar injection

    params = [
        # Standard jet and medium parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),

        # Magnetar injection parameters
        ParamDef("L0",      1e42,  1e48,  Scale.LOG),     # Initial luminosity [erg/s]
        ParamDef("t0",        10,  1000,  Scale.LOG),     # Spin-down timescale [s]
        ParamDef("q",        1.5,   3.0,  Scale.LINEAR),  # Power-law index

        # Standard microphysics (same as default)
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

.. note::
    **Magnetar Injection Profile:** L(t) = L0 × (1 + t/t0)^(-q) for θ < θc


**Magnetar with Structured Jet**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "ism"
    cfg.jet = "powerlaw"      # Structured jet
    cfg.magnetar = True

    params = [
        # Power-law jet with magnetar
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.3,  Scale.LINEAR),
        ParamDef("k_e",      1.5,   3.0,  Scale.LINEAR),
        ParamDef("k_g",      1.5,   3.0,  Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),

        # Magnetar parameters
        ParamDef("L0",      1e42,  1e48,  Scale.LOG),
        ParamDef("t0",        10,  1000,  Scale.LOG),
        ParamDef("q",        1.5,   3.0,  Scale.LINEAR),

        # Standard microphysics
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
    ]

Complex Model Combinations
^^^^^^^^^^^^^^^^^^^^^^^^^^

**Full Physics: Structured Jet + Stratified Medium + Reverse Shock + IC + Magnetar**

.. code-block:: python

    cfg = Setups()
    cfg.medium = "wind"       # Stratified medium
    cfg.jet = "gaussian"      # Structured jet
    cfg.rvs_shock = True      # Reverse shock
    cfg.fwd_ssc = True        # Forward SSC
    cfg.rvs_ssc = True        # Reverse SSC
    cfg.ssc_cooling = True     # IC cooling
    cfg.kn = True             # Klein-Nishina
    cfg.magnetar = True       # Energy injection

    params = [
        # Gaussian jet
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    50,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.2,  Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("tau",        1,   100,  Scale.LOG),

        # Stratified medium
        ParamDef("A_star",  1e-4,   1.0,  Scale.LOG),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("n0",      1e-2,    50,  Scale.LOG),

        # Magnetar injection
        ParamDef("L0",      1e42,  1e48,  Scale.LOG),
        ParamDef("t0",        10,  1000,  Scale.LOG),
        ParamDef("q",        1.5,   3.0,  Scale.LINEAR),

        # Forward shock microphysics
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),

        # Reverse shock microphysics
        ParamDef("p_r",      2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e_r", 1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B_r", 1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e_r",   0.1,   1.0,  Scale.LINEAR),
    ]

.. warning::
    **Complex Model Considerations:**
    - Use coarser resolution initially: ``resolution=(0.2, 0.7, 7)``
    - Increase MCMC steps: ``total_steps=30000+``
    - More burn-in: ``burn_frac=0.4``
    - Consider parameter degeneracies in interpretation

Running MCMC
------------

Basic MCMC Execution
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Check data statistics before MCMC
    print(f"Total data points: {data.data_points_num()}")

    # Create fitter object
    fitter = Fitter(data, cfg, num_workers=8)  # Use 8 CPU cores

    # Run MCMC
    result = fitter.fit(
        param_defs=params,
        resolution=(0.3, 1, 10),     # (phi, theta, time) resolution
        total_steps=20000,           # Total MCMC steps
        burn_frac=0.3,               # Burn-in fraction
        thin=1                      # Thinning factor
    )


Analyzing Results
-----------------

Parameter Constraints
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Print best-fit parameters
    top_k_data = []
    for i in range(result.top_k_params.shape[0]):
        row = {'Rank': i+1, 'chi^2': f"{-2*result.top_k_log_probs[i]:.2f}"}
        for name, val in zip(result.labels, result.top_k_params[i]):
            row[name] = f"{val:.4f}"
        top_k_data.append(row)

    top_k_df = pd.DataFrame(top_k_data)
    print("Top-k parameters:")
    print(top_k_df.to_string(index=False))

Model Predictions
^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Generate model predictions with best-fit parameters
    t_model = np.logspace(2, 8, 200)
    nu_model = np.array([1e9, 5e14, 2e17])  # Radio, optical, X-ray

    # Light curves at specific frequencies
    lc_model = fitter.flux_density_grid(result.best_params, t_model, nu_model)

    # Spectra at specific times
    nu_spec = np.logspace(8, 20, 100)
    times_spec = [1000, 10000]
    spec_model = fitter.flux_density_grid(result.best_params, times_spec, nu_spec)

    # Frequency-integrated flux (broadband light curves)
    # Useful for comparing with instruments like Swift/BAT, Fermi/LAT
    nu_min_broad = 1e17  # Lower frequency bound [Hz]
    nu_max_broad = 1e19  # Upper frequency bound [Hz]
    num_freq_points = 5  # Number of frequency points for integration

    flux_integrated = fitter.flux(result.best_params, t_model,
                                  nu_min_broad, nu_max_broad, num_freq_points)

Visualization
^^^^^^^^^^^^^

.. code-block:: python

    import matplotlib.pyplot as plt
    import corner

    # Corner plot for parameter correlations
    fig = corner.corner(
        flat_chain,
        labels=result.labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12}
    )
    plt.savefig("corner_plot.png", dpi=300, bbox_inches='tight')

    # Light curve comparison
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    colors = ['blue', 'orange', 'red']

    for i, (nu, color) in enumerate(zip(nu_model, colors)):
        ax = axes[i]

        # Plot data (if available)
        # ax.errorbar(t_data, flux_data, flux_err, fmt='o', color=color)

        # Plot model
        ax.loglog(t_model, lc_model[i], '-', color=color, linewidth=2)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Flux Density [erg/cm²/s/Hz]')
        ax.set_title(f'ν = {nu:.1e} Hz')

    plt.tight_layout()
    plt.savefig("lightcurve_fit.png", dpi=300, bbox_inches='tight')


Troubleshooting
---------------

For comprehensive troubleshooting help including MCMC convergence issues, data selection problems, memory optimization, and performance tuning, see :doc:`troubleshooting`.
