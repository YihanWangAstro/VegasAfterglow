MCMC Parameter Fitting
======================

.. contents:: Table of Contents
   :local:
   :depth: 3

Overview
--------

.. note::
   MCMC fitting requires additional dependencies. Install them with:
   ``pip install VegasAfterglow[mcmc]``

VegasAfterglow provides a comprehensive MCMC framework for parameter estimation of gamma-ray burst (GRB) afterglow models. The framework supports all built-in jet models, ambient medium configurations, radiation processes, and complex multi-wavelength datasets.

Key Features:

- **Bilby/Dynesty nested sampling** for efficient parameter estimation
- **All model types supported**: TophatJet, GaussianJet, PowerLawJet, TwoComponentJet, StepPowerLawJet
- **Complete physics**: Forward shock, reverse shock, synchrotron, inverse Compton, magnetar injection
- **Flexible data handling**: Light curves and spectra with optional weighting

.. seealso::
   :doc:`parameter_reference` for parameter details and typical ranges. :doc:`physics` for the underlying physical models.

MCMC Best Practices
-------------------

Start Simple, Add Complexity Gradually
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The most common mistake in GRB afterglow fitting is starting with an overly complex model. **Always begin with the simplest physically motivated model** and only add complexity when the data clearly demands it.

**Recommended Progression:**

1. **Start with TopHat + ISM + Forward Shock Only**

   .. code-block:: python

       fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")
       # All other physics options default to False

   This gives you ~7-8 free parameters. Run MCMC and examine the residuals.

2. **Check if residuals suggest additional physics**

   Examine the fit quality and residual structure. Systematic deviations at specific times or frequencies may indicate missing physics.

3. **Add ONE component at a time**

   Never jump from a simple model to enabling everything. Each addition should be justified by improved fit statistics (e.g., Bayesian evidence comparison).

Why Complex Models Are Problematic
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**1. Parameter Degeneracies**

More parameters create more degeneracies. For example:

- ``E_iso`` and ``n_ism`` are degenerate in flux normalization
- ``eps_e`` and ``eps_B`` trade off in spectral shape
- ``theta_c`` and ``theta_v`` correlate for off-axis observers
- Reverse shock parameters can mimic forward shock with different microphysics

With a complex model, the MCMC may find multiple solutions that fit equally well but have very different physical interpretations.

**2. Computational Cost**

Each additional physics module increases computation time. A model with all physics enabled can be significantly slower than a basic forward-shock-only model.

**3. Overfitting Risk**

With enough free parameters, you can fit noise. A model that fits your data perfectly but has 15+ parameters may not be physically meaningful. Use Bayesian evidence (from dynesty) to compare models:

.. code-block:: python

    # Compare two models using log evidence
    result_simple = fitter_simple.fit(params_simple, sampler="dynesty", ...)
    result_complex = fitter_complex.fit(params_complex, sampler="dynesty", ...)

    # Bayes factor
    log_BF = result_complex.bilby_result.log_evidence - result_simple.bilby_result.log_evidence

    # Interpretation:
    # log_BF < 1: No preference (stick with simple model)
    # 1 < log_BF < 3: Weak preference for complex
    # 3 < log_BF < 5: Moderate preference
    # log_BF > 5: Strong preference for complex model

**4. Non-Physical Solutions**

Complex models can converge to non-physical parameter combinations. Always check:

- Is ``eps_e + eps_B < 1``? (energy conservation)
- Is ``p > 2``? (required for finite electron energy)
- Are microphysics parameters consistent between forward/reverse shocks?
- Does the inferred jet energy make sense given the host galaxy and redshift?

When to Use Complex Models
^^^^^^^^^^^^^^^^^^^^^^^^^^

Complex models are justified when:

- **Clear observational signatures** that simple models cannot explain after careful analysis
- **Comparison with similar GRBs** where the additional physics was robustly established

.. important::
    **Golden Rule**: If you cannot clearly explain WHY each physics component is needed based on your data, you probably don't need it.

Basic Setup
-----------

Fitter Constructor
^^^^^^^^^^^^^^^^^^

The ``Fitter`` constructor accepts keyword arguments for source properties, model selection, physics options, and numerical parameters:

.. code-block:: python

    fitter = Fitter(
        # Source properties
        z=1.58,                   # Redshift
        lumi_dist=3.364e28,       # Luminosity distance [cm]

        # Model selection (see sections below for all options)
        jet="powerlaw",           # Jet structure type
        medium="wind",            # Ambient medium type

        # Physics options
        rvs_shock=True,           # Include reverse shock
        fwd_ssc=True,             # Forward shock inverse Compton
        rvs_ssc=False,            # Reverse shock inverse Compton
        ssc_cooling=True,         # IC cooling effects
        kn=True,                  # Klein-Nishina corrections
        magnetar=True,            # Magnetar energy injection

        # Numerical parameters
        rtol=1e-5,                # Numerical tolerance
        resolution=(0.15, 0.5, 10),  # Grid resolution (phi, theta, t)
    )

Setting up Data and the Fitter
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Fitter`` class handles both model configuration and observational data. Add light curves, spectra, and broadband flux directly to the fitter:

.. code-block:: python

    import numpy as np
    import pandas as pd
    from VegasAfterglow import Fitter, ParamDef, Scale

    # Create the fitter with model configuration
    fitter = Fitter(
        z=1.58,
        lumi_dist=3.364e28,
        jet="tophat",
        medium="ism",
    )

    # Method 1: Add data directly
    t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
    flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # erg/cm²/s/Hz
    flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]   # Error bars

    # Add light curve at R-band frequency
    fitter.add_flux_density(nu=4.84e14, t=t_data,
                            f_nu=flux_data, err=flux_err)  # All quantities in CGS units

    # Optional: Add weights for systematic uncertainties, normalization handled internally
    weights = np.ones_like(t_data)  # Equal weights
    fitter.add_flux_density(nu=2.4e17, t=t_data,
                            f_nu=flux_data, err=flux_err,
                            weights=weights)  # All quantities in CGS units

    # Method 2: Add frequency-integrated light curve (broadband flux)
    # For instruments with wide frequency coverage (e.g., BAT, LAT, Fermi)
    nu_min = 1e17  # Lower frequency bound [Hz]
    nu_max = 1e19  # Upper frequency bound [Hz]
    num_points = 5  # Number of frequency points for integration

    fitter.add_flux(nu_min=nu_min, nu_max=nu_max,
                    num_points=num_points, t=t_data,
                    flux=flux_data, err=flux_err,
                    weights=weights)  # All quantities in CGS units

    # Method 3: Load from CSV files
    bands = [2.4e17, 4.84e14, 1.4e14]  # X-ray, optical, near-IR
    lc_files = ["data/xray.csv", "data/optical.csv", "data/nir.csv"]

    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        fitter.add_flux_density(nu=nu, t=df["t"],
                                f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

    # Add spectra at specific times
    spec_times = [1000, 10000]  # seconds
    spec_files = ["data/spec_1000s.csv", "data/spec_10000s.csv"]

    for t, fname in zip(spec_times, spec_files):
        df = pd.read_csv(fname)
        fitter.add_spectrum(t=t, nu=df["nu"],
                            f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

Data Selection and Optimization
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Smart Data Subsampling with logscale_screen**

For large datasets or densely sampled observations, using all available data points can lead to computational inefficiency and biased parameter estimation. The ``logscale_screen`` utility provides intelligent data subsampling that maintains the essential information content while reducing computational overhead.

.. code-block:: python

    from VegasAfterglow import logscale_screen

    # Example: Large dense dataset
    t_dense = np.logspace(2, 7, 1000)  # 1000 time points
    flux_dense = np.random.lognormal(-60, 0.5, 1000)  # Dense flux measurements
    flux_err_dense = 0.1 * flux_dense

    # Subsample using logarithmic screening
    # This selects ~5*5=25 representative points across 5 decades in time
    indices = logscale_screen(t_dense, data_density=5)

    # Add only the selected subset
    fitter.add_flux_density(nu=5e14,
                            t=t_dense[indices],
                            f_nu=flux_dense[indices],
                            err=flux_err_dense[indices])

**Why logscale_screen is Important:**

1. **Prevents Oversampling Bias**: Dense data clusters can dominate the chi-squared calculation, causing the MCMC to over-fit specific frequency bands or time periods.

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
    - **Single-epoch spectra**: Broadband spectra at one time can dominate multi-epoch light curves in chi-squared space

    **Solution**: Use ``logscale_screen()`` for manual temporal reduction of over-sampled bands.

**Handling Measurement Uncertainties**

Error bars strongly influence MCMC fitting. Improperly characterized uncertainties can lead to biased parameter estimates or misleading confidence intervals.

**The Problem with Small Error Bars**

Data points with very small uncertainties dominate the chi-squared calculation:

.. math::

    \chi^2 = \sum_i \frac{(F_{\rm obs,i} - F_{\rm model,i})^2}{\sigma_i^2}

A single point with sigma = 0.01 mJy contributes 100x more to chi-squared than a point with sigma = 0.1 mJy. This causes the MCMC to:

- Heavily weight high-precision data at the expense of other observations
- Potentially fit noise or systematic effects in the "precise" data
- Produce artificially tight parameter constraints that don't reflect true uncertainties

**Common Scenarios:**

1. **Heterogeneous data quality**: X-ray data from different telescopes may have vastly different calibration uncertainties
2. **Underestimated errors**: Published error bars sometimes only include statistical uncertainty, missing systematics
3. **Calibration offsets**: Different instruments may have 5-20% flux calibration offsets not captured in quoted errors

**Solutions:**

1. **Add systematic error floors**

   For data where quoted errors seem unrealistically small:

   .. code-block:: python

       # Add 10% systematic floor to error bars
       systematic_floor = 0.1  # 10% of flux
       err_with_floor = np.sqrt(err**2 + (systematic_floor * flux)**2)

       fitter.add_flux_density(nu=nu, t=t, f_nu=flux, err=err_with_floor)

2. **Use weights to balance bands**

   Down-weight bands with suspiciously small errors:

   .. code-block:: python

       # Reduce contribution of high-precision optical data
       optical_weight = 0.5  # Effectively doubles the error contribution
       fitter.add_flux_density(nu=5e14, t=t_optical, f_nu=flux_optical,
                              err=err_optical, weights=optical_weight * np.ones_like(t_optical))

3. **Inspect residuals by band**

   After fitting, check if residuals are consistent across bands. If one band shows systematic deviations while others don't, the errors for that band may be mischaracterized.

.. tip::
    **Rule of thumb**: If your best-fit chi-squared/dof >> 1, either the model is inadequate OR the error bars are underestimated. If chi-squared/dof << 1, the errors may be overestimated.

Model Configurations
--------------------

.. seealso::
   :doc:`parameter_reference` for all available parameters. :doc:`examples` for direct model calculation examples.

TopHat + ISM (Default)
^^^^^^^^^^^^^^^^^^^^^^

The default configuration uses a top-hat jet in a uniform ISM environment with forward shock synchrotron emission:

.. code-block:: python

    # Basic configuration
    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="powerlaw", medium="ism")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="gaussian", medium="ism")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="two_component", medium="ism")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="step_powerlaw", medium="ism")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

    params = [
        # Standard jet parameters (same as default)
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    10,   500,  Scale.LOG),
        ParamDef("theta_c", 0.01,   0.5,  Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),

        # Stratified medium (wind -> ISM)
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

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="wind")

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
    - **n0 = infinity**: Pure wind profile from center
    - **A_star > 0, n0 < infinity**: ISM-wind-ISM stratification
    - **A_star > 0, n0 = infinity**: Wind-ISM stratification

    **Density Profile:** Inner (r < r1): n = n0, Middle (r1 < r < r2): n proportional to A_star/r^2, Outer (r > r2): n = n_ism

Reverse Shock
^^^^^^^^^^^^^

**Basic Reverse Shock**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        rvs_shock=True,
    )

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

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="gaussian", medium="ism",
        rvs_shock=True,
    )

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

For the physical details of inverse Compton and Klein-Nishina corrections, see :doc:`physics`.

**Forward Shock Inverse Compton**

.. code-block:: python

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        fwd_ssc=True,
        ssc_cooling=True,
        kn=True,
    )

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

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        rvs_shock=True,
        fwd_ssc=True,
        rvs_ssc=True,
        ssc_cooling=True,
        kn=True,
    )

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

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
        magnetar=True,
    )

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

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="gaussian", medium="wind",
        rvs_shock=True,
        fwd_ssc=True,
        rvs_ssc=True,
        ssc_cooling=True,
        kn=True,
        magnetar=True,
    )

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
    - Increase live points: ``nlive=1000+``
    - Lower dlogz for better convergence: ``dlogz=0.05``
    - Consider parameter degeneracies in interpretation

Running MCMC
------------

Choosing a Sampler
^^^^^^^^^^^^^^^^^^

The choice of sampler significantly affects both the quality of results and computational efficiency. Here's guidance on when to use each:

**Emcee (Ensemble MCMC) - Recommended for Most Cases**

.. code-block:: python

    result = fitter.fit(params, sampler="emcee", nsteps=10000, nburn=2000, ...)

Use emcee when:

- You need **fast exploration** of parameter space
- The posterior is expected to be **unimodal** (single peak)
- You're doing **quick preliminary fits** or production runs

VegasAfterglow uses an optimized emcee configuration:

- **Custom proposal moves**: DEMove (70%) + DESnookerMove (30%) for better mixing than default stretch move
- **Thread-based parallelism**: Likelihood evaluations are parallelized via ``ThreadPoolExecutor`` with the GIL released during C++ computation
- **Automatic nwalkers**: Optimized based on parameter count and CPU cores

.. note::
    **Parallelization**: For emcee, set ``npool`` to the number of CPU cores. Each thread creates a ``Model`` and calls ``flux_density()`` with the GIL released during C++ computation, achieving near-native parallelism.

**Dynesty (Nested Sampling) - For Model Comparison**

.. code-block:: python

    result = fitter.fit(params, sampler="dynesty", nlive=1000, dlogz=0.5, ...)

Use dynesty when:

- You need **Bayesian evidence** for model comparison
- The posterior may be **multimodal** (multiple peaks)
- You want to compare different physics configurations

VegasAfterglow uses optimized dynesty settings:

- **rslice sampling**: More efficient than random walk for high dimensions
- **Thread pool**: Uses Python threading (not multiprocessing) for better performance with C++ backend
- **Optimal queue size**: Automatically calculated based on nlive and npool

.. note::
    For dynesty, ``npool`` controls the number of parallel threads. If not specified, it defaults to the number of CPU cores on your machine.

**Recommended Workflow**

1. **Quick exploration**: Use emcee with moderate steps (5000-10000) to find approximate parameter region
2. **Model comparison**: Use dynesty when comparing different physics (forward shock only vs. with reverse shock, etc.)
3. **Production run**: Use emcee with longer chains (20000-50000 steps) for final parameter constraints

Fitter.fit() Interface
^^^^^^^^^^^^^^^^^^^^^^

The ``Fitter.fit()`` method provides a unified interface for parameter estimation using different sampling algorithms via bilby.

**Method Signature**

.. code-block:: python

    def fit(
        self,
        param_defs: Sequence[ParamDef],      # Parameter definitions
        sampler: str = "emcee",              # Sampler algorithm
        resolution: Tuple = None,            # Override constructor resolution
        npool: int = None,                   # Number of parallel threads (default: all cores)
        top_k: int = 10,                     # Number of best fits to return
        outdir: str = "bilby_output",        # Output directory
        label: str = "afterglow",            # Run label
        clean: bool = True,                  # Clean up intermediate files
        resume: bool = False,                # Resume previous run
        log_prior_fn: Callable = None,       # Custom log-prior (emcee only)
        log_likelihood_fn: Callable = None,  # Custom log-likelihood
        priors: dict = None,                 # Custom bilby priors (dynesty only)
        **sampler_kwargs                     # Sampler-specific parameters
    ) -> FitResult

**Common Parameters**

- ``param_defs``: List of ``ParamDef`` objects defining parameters to fit
    - Required parameter for all samplers
    - See "Parameter Definition" section for details

- ``resolution``: Optional tuple of (phi_res, theta_res, t_res)
    - Overrides the resolution set in the ``Fitter`` constructor for this run
    - If ``None``, uses the constructor value (default: ``(0.15, 0.5, 10)``)
    - Example: ``(0.1, 1, 15)`` for higher accuracy

- ``sampler``: Sampling algorithm to use
    - ``"emcee"``: Affine-invariant MCMC ensemble sampler (recommended for speed)
    - ``"dynesty"``: Dynamic nested sampling (for Bayesian evidence)
    - ``"nestle"``: Alternative nested sampling implementation
    - ``"cpnest"``: Nested sampling with parallel tempering
    - ``"pymultinest"``: MultiNest nested sampling algorithm
    - ``"ultranest"``: Nested sampling with slice sampling
    - Other samplers supported by bilby - see `bilby documentation <https://lscsoft.docs.ligo.org/bilby/samplers.html>`_ for details

- ``npool``: Number of parallel threads
    - Set to number of CPU cores for best performance
    - Default: number of CPU cores
    - Example: ``npool=8`` on an 8-core machine

- ``top_k``: Number of best-fit parameter sets to extract
    - Useful for exploring multiple local maxima
    - Default: 10
    - Results sorted by log probability (highest first)

- ``outdir``: Directory for output files
    - Default: ``"bilby_output"``
    - Contains checkpoint files, result files, and plots

- ``label``: Run identifier
    - Default: ``"afterglow"``
    - Used in output filenames: ``{label}_result.json``

- ``clean``: Whether to remove intermediate files
    - Default: ``True``
    - Set ``False`` to keep checkpoint files for inspection

- ``resume``: Resume from previous run
    - Default: ``False``
    - Requires matching ``outdir`` and ``label`` from previous run

**Emcee-Specific Parameters** (``sampler="emcee"``)

- ``nsteps``: Number of MCMC steps per walker
    - Default: 5000
    - More steps = better convergence but longer runtime
    - Typical range: 5000-50000

- ``nburn``: Number of burn-in steps to discard
    - Default: 1000
    - Should be long enough to reach equilibrium
    - Rule of thumb: 10-20% of ``nsteps``

- ``thin``: Thinning factor (save every nth sample)
    - Default: 1 (save all samples)
    - Use to reduce autocorrelation in chain
    - Example: ``thin=10`` saves every 10th step

- ``nwalkers``: Number of ensemble walkers
    - Default: ``2 * n_params`` (automatically set)
    - More walkers = better exploration but more computation
    - Minimum: ``2 * n_params``

**Dynesty-Specific Parameters** (``sampler="dynesty"``)

- ``nlive``: Number of live points
    - Default: 500 (recommended: 1000 for production runs)
    - More points = better evidence estimate but slower
    - Typical range: 500-2000

- ``dlogz``: Evidence tolerance (stopping criterion)
    - Default: 0.1 (recommended: 0.5 for faster convergence)
    - Smaller values = more thorough sampling
    - Typical range: 0.01-0.5

- ``sample``: Sampling method
    - Default: ``"rwalk"`` (random walk)
    - Other options: ``"slice"``, ``"rslice"``, ``"unif"``
    - ``"rwalk"`` works well for most problems

- ``walks``: Number of random walk steps
    - Default: 100
    - Only relevant for ``sample="rwalk"``

- ``maxmcmc``: Maximum MCMC steps for slice sampling
    - Default: 5000
    - Only relevant for slice samplers

**Other Samplers**

VegasAfterglow supports all samplers available in bilby, including:

- ``"nestle"``: Nested sampling (similar to dynesty)
- ``"cpnest"``: Nested sampling with parallel tempering
- ``"pymultinest"``: MultiNest algorithm (requires separate installation)
- ``"ultranest"``: Advanced nested sampling with slice sampling
- ``"ptemcee"``: Parallel-tempered MCMC
- ``"kombine"``: Clustered KDE proposal MCMC
- ``"pypolychord"``: PolyChord nested sampling

For sampler-specific parameters, refer to the `bilby sampler documentation <https://bilby-dev.github.io/bilby/>`_.
Pass sampler-specific parameters via ``**sampler_kwargs`` in the ``fit()`` method.

**Example with alternative sampler:**

.. code-block:: python

    # Using nestle sampler
    result = fitter.fit(
        params,
        resolution=(0.15, 0.5, 10),
        sampler="nestle",
        nlive=1000,           # nestle-specific parameter
        method='multi',       # nestle-specific: 'classic', 'single', or 'multi'
        npool=8,
        top_k=10,
    )

**Return Value: FitResult**

The method returns a ``FitResult`` object with the following attributes:

- ``samples``: Posterior samples array
    - Shape: ``[n_samples, 1, n_params]``
    - In transformed space (log10 for LOG-scale parameters)

- ``labels``: Parameter names as used in sampling
    - LOG-scale parameters have ``log10_`` prefix
    - Example: ``["log10_E_iso", "log10_Gamma0", "theta_c", "p", ...]``

- ``latex_labels``: LaTeX-formatted labels for plotting
    - Example: ``["$\\log_{10}(E_{\\rm iso})$", "$\\log_{10}(\\Gamma_0)$", "$\\theta_c$", ...]``
    - Use these for corner plots and visualizations

- ``top_k_params``: Best-fit parameter values
    - Shape: ``[top_k, n_params]``
    - Sorted by log probability (highest first)
    - In transformed space (log10 for LOG-scale parameters)

- ``top_k_log_probs``: Log probabilities for top-k fits
    - Shape: ``[top_k]``
    - Can convert to chi-squared via: ``chi2 = -2 * log_prob``

- ``log_probs``: Log probabilities for all samples
    - Shape: ``[n_samples, 1]``

- ``bilby_result``: Full bilby Result object
    - Contains additional diagnostics and metadata
    - Use for advanced analysis and plotting
    - Access via: ``result.bilby_result.plot_corner()``

Basic MCMC Execution
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Create fitter object and add data (see above for data loading examples)
    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")
    # ... add data via fitter.add_flux_density(), fitter.add_spectrum(), etc.

    # Option 1 (Recommended): Nested sampling with dynesty (computes Bayesian evidence, robust for multimodal posteriors)
    result = fitter.fit(
        params,
        resolution=(0.15, 0.5, 10),     # Grid resolution (phi, theta, t)
        sampler="dynesty",             # Nested sampling algorithm
        nlive=1000,                    # Number of live points
        walks=100,                     # Number of random walks per live point
        dlogz=0.5,                     # Stopping criterion (evidence tolerance)
        npool=8,                       # Number of parallel threads
        top_k=10,                      # Number of best-fit parameters to return
    )

    # Option 2: MCMC with emcee (faster, good for unimodal posteriors, hard to converge for multimodal posteriors)
    result = fitter.fit(
        params,
        resolution=(0.15, 0.5, 10),     # Grid resolution (phi, theta, t)
        sampler="emcee",               # MCMC sampler
        nsteps=50000,                  # Number of steps per walker
        nburn=10000,                   # Burn-in steps to discard
        thin=1,                        # Save every nth sample
        npool=8,                       # Number of parallel threads
        top_k=10,                      # Number of best-fit parameters to return
    )

**Important Notes:**
    - Parameters with ``Scale.LOG`` are sampled as ``log10_<name>`` (e.g., ``log10_E_iso``)
    - The sampler works in log10 space for LOG-scale parameters, then transforms back to physical values
    - Use ``npool`` to parallelize likelihood evaluations across multiple CPU cores
    - ``result.latex_labels`` provides properly formatted labels for corner plots
    - For emcee, total samples = ``nwalkers * (nsteps - nburn) / thin``
    - For dynesty, sample count depends on convergence (not fixed)


Analyzing Results
-----------------

Parameter Constraints
^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Print best-fit parameters
    print("Top-k parameters:")
    header = f"{'Rank':>4s}  {'chi^2':>10s}  " + "  ".join(f"{name:>10s}" for name in result.labels)
    print(header)
    print("-" * len(header))
    for i in range(result.top_k_params.shape[0]):
        chi2 = -2 * result.top_k_log_probs[i]
        vals = "  ".join(f"{val:10.4f}" for val in result.top_k_params[i])
        print(f"{i+1:4d}  {chi2:10.2f}  {vals}")

Model Predictions
^^^^^^^^^^^^^^^^^

.. code-block:: python

    # Generate model predictions with best-fit parameters
    t_model = np.logspace(2, 8, 200)
    nu_model = np.array([1e9, 5e14, 2e17])  # Radio, optical, X-ray

    # Light curves at specific frequencies
    lc_model = fitter.flux_density_grid(result.top_k_params[0], t_model, nu_model)

    # Spectra at specific times
    nu_spec = np.logspace(8, 20, 100)
    times_spec = [1000, 10000]
    spec_model = fitter.flux_density_grid(result.top_k_params[0], times_spec, nu_spec)

    # Frequency-integrated flux (broadband light curves)
    # Useful for comparing with instruments like Swift/BAT, Fermi/LAT
    nu_min_broad = 1e17  # Lower frequency bound [Hz]
    nu_max_broad = 1e19  # Upper frequency bound [Hz]
    num_freq_points = 5  # Number of frequency points for integration

    flux_integrated = fitter.flux(result.top_k_params[0], t_model,
                                  nu_min_broad, nu_max_broad, num_freq_points)

Visualization
^^^^^^^^^^^^^

.. code-block:: python

    import matplotlib.pyplot as plt
    import corner

    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])

    # Corner plot for parameter correlations
    fig = corner.corner(
        flat_chain,
        labels=result.latex_labels,
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
        ax.set_ylabel('Flux Density [erg/cm^2/s/Hz]')
        ax.set_title(f'nu = {nu:.1e} Hz')

    plt.tight_layout()
    plt.savefig("lightcurve_fit.png", dpi=300, bbox_inches='tight')


Advanced Fitting
----------------

VegasAfterglow's ``Fitter`` is built directly on the ``Model`` class, giving you full control over the fitting process. You can customize:

- **Priors**: Use informative or custom prior distributions
- **Likelihood**: Implement non-standard likelihood functions (e.g., upper limits, systematic uncertainties)
- **Jet profiles**: Define arbitrary angular energy and Lorentz factor distributions
- **Medium profiles**: Define custom circumburst density structures

The ``Fitter`` evaluates models using Python's ``ThreadPoolExecutor``, where each thread constructs a ``Model`` and calls ``flux_density()`` with the GIL released during C++ computation. This provides near-native parallelism with full Python flexibility.

Custom Priors
^^^^^^^^^^^^^

**Emcee: Log-Prior Function**

For emcee, pass a ``log_prior_fn`` that takes an array of walker positions and returns log-prior values:

.. code-block:: python

    import numpy as np
    from VegasAfterglow import Fitter, ParamDef, Scale

    fitter = Fitter(
        z=1.58, lumi_dist=3.364e28,
        jet="tophat", medium="ism",
    )

    params = [
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

    fitter.add_flux_density(nu=4.84e14, t=t_data, f_nu=flux_data, err=flux_err)

    def log_prior(samples):
        """Gaussian prior on p centered at 2.3 with sigma=0.1.

        Args:
            samples: Array of shape (nwalkers, ndim)

        Returns:
            Array of shape (nwalkers,) with log-prior values
        """
        log_p = np.zeros(samples.shape[0])

        # p is the 6th free parameter (index 5 in sampler space)
        p_values = samples[:, 5]
        log_p += -0.5 * ((p_values - 2.3) / 0.1) ** 2

        return log_p

    result = fitter.fit(
        params,
        sampler="emcee",
        nsteps=10000,
        nburn=2000,
        log_prior_fn=log_prior,
    )

.. note::
    The ``log_prior_fn`` receives parameters in **sampler space**: ``LOG``-scale parameters are passed as ``log10(value)``. The prior is added to the log-likelihood, so returning ``-np.inf`` rejects the sample.

**Bilby/Dynesty: Prior Distributions**

For bilby-based samplers (dynesty, nestle, etc.), pass a ``priors`` dictionary mapping parameter labels to ``bilby.core.prior.Prior`` objects:

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

    result = fitter.fit(
        params,
        sampler="dynesty",
        nlive=1000,
        priors=custom_priors,
    )

Parameters not included in the ``priors`` dictionary automatically get uniform priors based on the bounds in ``ParamDef``.

.. important::
    When using ``LOG``-scale parameters, the prior keys must use the ``log10_`` prefix (e.g., ``log10_eps_e``, not ``eps_e``), since the sampler operates in log-space.

Custom Likelihood
^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^

The ``Fitter`` accepts a custom jet factory function via the ``jet`` keyword argument.
This lets you define arbitrary angular profiles using the same ``Ejecta`` class used for
direct model calculations (see :doc:`examples`).

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
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    50,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.15, Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("E_iso_w", 1e48,  1e52,  Scale.LOG),
        ParamDef("theta_w",  0.1,   0.5,  Scale.LINEAR),
        ParamDef("tau",        1,   1e3,  Scale.LOG),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
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
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    50,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.15, Scale.LINEAR),
        ParamDef("theta_v",    0,   0.5,  Scale.LINEAR),
        ParamDef("L0",      1e44,  1e48,  Scale.LOG),       # Injection luminosity [erg/s]
        ParamDef("t_sd",      10,  1e4,   Scale.LOG),       # Spin-down timescale [s]
        ParamDef("M_dot0",  1e20,  1e26,  Scale.LOG),       # Mass injection rate [g/s]
        ParamDef("tau",        1,   1e3,  Scale.LOG),
        # ... other standard params (n_ism, p, eps_e, eps_B, xi_e)
    ]

    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet=jet_with_injection)

See :doc:`examples` for more user-defined jet examples.

Custom Medium Profiles
^^^^^^^^^^^^^^^^^^^^^^

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

The ``Medium`` class takes a callable ``rho(phi, theta, r)`` that returns the mass density [g/cm³] at position (phi, theta, r). See :doc:`examples` for more user-defined medium examples.

**Custom MCMC Parameters**

When using custom jet or medium functions, you can define arbitrary MCMC parameters beyond the standard set. Simply include them in the ``ParamDef`` list -- the fitter automatically uses a Python-based parameter transformer when it detects non-standard parameter names:

.. code-block:: python

    mc_params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),
        ParamDef("Gamma0",    50,   500,  Scale.LOG),
        ParamDef("theta_c", 0.02,   0.15, Scale.LINEAR),
        ParamDef("theta_v",    0,     0,  Scale.FIXED),
        ParamDef("n_ism",   1e-3,   100,  Scale.LOG),
        ParamDef("r_scale", 1e16,  1e20,  Scale.LOG),     # custom parameter
        ParamDef("p",        2.1,   2.8,  Scale.LINEAR),
        ParamDef("eps_e",   1e-3,   0.5,  Scale.LOG),
        ParamDef("eps_B",   1e-5,   0.1,  Scale.LOG),
        ParamDef("xi_e",     0.1,   1.0,  Scale.LINEAR),
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
        log_prior_fn=log_prior,
        log_likelihood_fn=student_t_likelihood,
    )


Speeding Up Custom Profiles with ``@gil_free``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The custom jet and medium examples above use plain Python callbacks. These work correctly
for single-threaded fitting, but become slow in multi-threaded MCMC because Python's
Global Interpreter Lock prevents true parallel execution of Python callbacks.

The ``@gil_free`` decorator compiles your profile function to native code
using `numba <https://numba.pydata.org/>`_, enabling full thread parallelism:

.. code-block:: bash

    pip install numba

The two differences from plain Python callbacks are:

1. Add the ``@gil_free`` decorator
2. Pass MCMC parameters as explicit function arguments (after the spatial coordinates)
   instead of capturing them from the enclosing scope

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
        ParamDef("E_iso",    1e51,  1e54,  Scale.LOG,    1e53),
        ParamDef("Gamma0",      5,  1000,  Scale.LOG,      20),
        ParamDef("theta_c",  0.01,   0.5,  Scale.LINEAR,  0.1),
        ParamDef("p",           2,     3,  Scale.LINEAR,  2.3),
        ParamDef("eps_e",    1e-2,   0.3,  Scale.LOG,    0.05),
        ParamDef("eps_B",    1e-4,   0.3,  Scale.LOG,    0.03),
        ParamDef("xi_e",     1e-3,   0.1,  Scale.LOG,    0.01),
        ParamDef("A_star",   1e-3,    10,  Scale.LOG,    0.05),
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
^^^^^^^^^^^^^^^^^

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
---------------

For comprehensive troubleshooting help including MCMC convergence issues, data selection problems, memory optimization, and performance tuning, see :doc:`troubleshooting`.

.. seealso::
   :doc:`validation` for code validation and comparison with other afterglow codes. :doc:`parameter_reference` for the complete parameter reference.
