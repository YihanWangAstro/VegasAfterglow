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
    data.add_light_curve(nu_cgs=4.84e14, t_cgs=t_data, 
                         Fnu_cgs=flux_data, Fnu_err=flux_err)

    # Optional: Add weights for systematic uncertainties, normalization handled internally
    weights = np.ones(len(t_data))  # Equal weights
    data.add_light_curve(nu_cgs=2.4e17, t_cgs=t_data, 
                         Fnu_cgs=flux_data, Fnu_err=flux_err, 
                         weights=weights)

    # Method 2: Load from CSV files
    bands = [2.4e17, 4.84e14, 1.4e14]  # X-ray, optical, near-IR
    lc_files = ["data/xray.csv", "data/optical.csv", "data/nir.csv"]
    
    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], 
                             Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

    # Add spectra at specific times
    spec_times = [1000, 10000]  # seconds
    spec_files = ["data/spec_1000s.csv", "data/spec_10000s.csv"]
    
    for t, fname in zip(spec_times, spec_files):
        df = pd.read_csv(fname)
        data.add_spectrum(t_cgs=t, nu_cgs=df["nu"], 
                          Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

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
    cfg.fwd_SSC = True        # Forward shock inverse Compton
    cfg.rvs_SSC = False       # Reverse shock inverse Compton  
    cfg.IC_cooling = True     # IC cooling effects
    cfg.KN = True             # Klein-Nishina corrections
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
    cfg.jet = "twocomponent"  # Two-component jet
    
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
    cfg.jet = "steppowerlaw"  # Step power-law jet
    
    params = [
        # Core component (uniform)
        ParamDef("E_iso_c", 1e51,  1e54,  Scale.LOG),     # Core energy
        ParamDef("Gamma0_c",  50,   500,  Scale.LOG),     # Core Lorentz factor
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
    cfg.fwd_SSC = True        # Forward shock SSC
    cfg.IC_cooling = True     # IC cooling effects
    cfg.KN = True             # Klein-Nishina corrections
    
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
    cfg.fwd_SSC = True        # Forward shock SSC
    cfg.rvs_SSC = True        # Reverse shock SSC
    cfg.IC_cooling = True
    cfg.KN = True
    
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
    cfg.fwd_SSC = True        # Forward SSC
    cfg.rvs_SSC = True        # Reverse SSC
    cfg.IC_cooling = True     # IC cooling
    cfg.KN = True             # Klein-Nishina
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
    print("Best-fit parameters:")
    for name, val in zip(result.labels, result.best_params):
        print(f"  {name}: {val:.4e}")

    # Compute credible intervals
    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
    medians = np.median(flat_chain, axis=0)
    lower = np.percentile(flat_chain, 16, axis=0)
    upper = np.percentile(flat_chain, 84, axis=0)

    print("\nParameter constraints (68% credible intervals):")
    for i, name in enumerate(result.labels):
        print(f"  {name}: {medians[i]:.3e} "
              f"(+{upper[i]-medians[i]:.3e}, -{medians[i]-lower[i]:.3e})")

Model Predictions
^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Generate model predictions with best-fit parameters
    t_model = np.logspace(2, 8, 200)
    nu_model = np.array([1e9, 5e14, 2e17])  # Radio, optical, X-ray

    # Light curves
    lc_model = fitter.specific_flux(result.best_params, t_model, nu_model)

    # Spectra at specific times
    nu_spec = np.logspace(8, 20, 100)
    times_spec = [1000, 10000]
    spec_model = fitter.specific_flux(result.best_params, times_spec, nu_spec)

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

Common Issues
^^^^^^^^^^^^^

**Poor Convergence**
  - Increase burn-in fraction
  - Use more walkers
  - Check parameter ranges
  - Start with coarser resolution

**Memory Errors**
  - Reduce resolution
  - Decrease number of workers
  - Use thinning (thin > 1)

**Slow Performance**
  - Reduce resolution initially
  - Use fewer walkers for exploration
  - Check data size and complexity

**Failed Fits**
  - Verify data format and units
  - Check parameter bounds
  - Ensure model configuration matches data

For more troubleshooting help, see :doc:`troubleshooting`.
