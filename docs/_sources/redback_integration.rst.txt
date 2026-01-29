Parameter Estimation with Redback
==================================

.. contents:: Table of Contents
   :local:
   :depth: 3

Overview
--------

VegasAfterglow is a **pure physics library** focused on fast, accurate afterglow calculations. For parameter estimation, MCMC fitting, and comprehensive transient analysis, VegasAfterglow models are integrated into `redback <https://github.com/nikhil-sarin/redback>`_.

Redback provides:

* **Unified interface** for all transient types (afterglows, kilonovae, supernovae, TDEs, etc.)
* **Multiple sampling algorithms** (dynesty, emcee, ultranest, etc.)
* **Data management** and open-access catalogs
* **Model comparison** and selection
* **Visualization tools** and plotting utilities
* **Prior specification** with bilby integration

Quick Start
-----------

Loading Data and Basic Fitting
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Redback provides multiple ways to load and work with GRB afterglow data. Here's a complete example:

.. code-block:: python

    import redback

    # Method 1: Load from Swift BAT+XRT catalogs (recommended for Swift GRBs)
    GRB = '070809'

    # Get flux/flux density data from Swift
    redback.get_data.get_bat_xrt_afterglow_data_from_swift(
        grb=GRB,
        data_mode="flux"  # or "flux_density", "magnitude", "luminosity"
    )
    # This creates a GRBDir with all available Swift data

    # Create transient object from Swift data
    afterglow = redback.afterglow.SGRB.from_swift_grb(
        name=GRB,
        data_mode='flux'
    )

    # Method 2: Load from open access catalogs
    afterglow = redback.afterglow.Afterglow.from_open_access_catalogue(
        GRB,
        data_mode='flux_density'
    )

    # Method 3: Load from your own data files
    import pandas as pd
    data = pd.read_csv('path/to/my_grb_data.csv')
    afterglow = redback.transient.Afterglow(
        name=GRB,
        data_mode='flux_density',
        time=data['time'].values,
        flux_density=data['flux'].values,
        flux_density_err=data['flux_err'].values,
        frequency=data['frequency'].values
    )

    # Fit with VegasAfterglow tophat model
    result = redback.fit_model(
        transient=afterglow,
        model='vegas_tophat',  # VegasAfterglow tophat jet model
        sampler='dynesty',
        nlive=1000
    )

    # Analyze results
    result.plot_corner()
    result.plot_lightcurve()
    result.plot_residuals()

**Data Modes and Formats:**

Redback supports multiple data modes for afterglow analysis:

* ``flux_density`` - Flux density [erg/cm²/s/Hz] at specific frequencies
* ``flux`` - Integrated flux [erg/cm²/s] over energy bands
* ``magnitude`` - Optical/NIR magnitudes in various filters
* ``luminosity`` - Source-frame luminosity (requires known redshift)

**Advanced Data Loading:**

For comprehensive data management capabilities, including:

* Multi-wavelength data from different instruments
* Custom data formats and preprocessing
* Data filtering and quality cuts
* Combining multiple datasets
* Time-domain sampling strategies

See the complete `redback data documentation <https://redback.readthedocs.io/en/latest/data.html>`_.

Available VegasAfterglow Models
--------------------------------

Redback includes multiple VegasAfterglow models for different jet structures and media:

Jet Models with ISM Medium
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``vegas_tophat`` - Uniform tophat jet with sharp boundary
* ``vegas_gaussian`` - Gaussian energy profile
* ``vegas_powerlaw`` - Power-law energy and Lorentz factor structure
* ``vegas_powerlaw_wing`` - Power-law wings without discontinuity
* ``vegas_step_powerlaw`` - Pencil beam core with power-law wings
* ``vegas_two_component`` - Narrow core + wide component

Jet Models with Wind Medium
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``vegas_tophat_wind`` - Tophat jet in stellar wind
* ``vegas_gaussian_wind`` - Gaussian jet in stellar wind
* ``vegas_powerlaw_wind`` - Power-law jet in stellar wind

Model Parameters
----------------

Each VegasAfterglow model requires specific parameters. Here are typical parameters for a tophat jet in ISM:

Common Parameters
^^^^^^^^^^^^^^^^^

* ``redshift`` - Source redshift
* ``thv`` - Viewing angle (radians)
* ``loge0`` - log10 isotropic equivalent energy (ergs)
* ``thc`` - Jet core opening angle (radians)
* ``lognism`` - log10 ISM number density (cm^-3)
* ``p`` - Electron power-law index
* ``logepse`` - log10 electron energy fraction
* ``logepsb`` - log10 magnetic field energy fraction
* ``g0`` - Initial Lorentz factor

Additional Parameters (via kwargs)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* ``spread`` - Enable jet spreading (default: False)
* ``reverse_shock`` - Enable reverse shock (default: False)
* ``reverse_logepse`` - Reverse shock electron fraction
* ``reverse_logepsb`` - Reverse shock magnetic fraction
* ``reverse_p`` - Reverse shock electron index

Example: Custom Priors
-----------------------

.. code-block:: python

    import bilby
    import redback

    # Load data
    afterglow = redback.afterglow.Afterglow.from_open_access_catalogue(
        'GRB070809', data_mode='flux_density'
    )

    # Define custom priors
    priors = bilby.prior.PriorDict()
    priors['redshift'] = bilby.prior.Uniform(0.01, 3, 'redshift')
    priors['thv'] = bilby.prior.Sine(maximum=np.pi/2, name='thv')
    priors['loge0'] = bilby.prior.Uniform(48, 54, 'loge0')
    priors['thc'] = bilby.prior.Uniform(0.01, 0.5, 'thc')
    priors['lognism'] = bilby.prior.Uniform(-5, 2, 'lognism')
    priors['p'] = bilby.prior.Uniform(2.0, 3.5, 'p')
    priors['logepse'] = bilby.prior.Uniform(-4, 0, 'logepse')
    priors['logepsb'] = bilby.prior.Uniform(-4, 0, 'logepsb')
    priors['g0'] = bilby.prior.Uniform(50, 1000, 'g0')

    # Run fit with custom priors
    result = redback.fit_model(
        transient=afterglow,
        model='vegas_tophat',
        sampler='dynesty',
        priors=priors,
        nlive=1000
    )

Example: Multi-band Fitting
----------------------------

.. code-block:: python

    import redback

    # Load multi-band data
    afterglow = redback.afterglow.Afterglow.from_open_access_catalogue(
        'GRB070809', data_mode='flux_density'
    )

    # Fit with VegasAfterglow model
    result = redback.fit_model(
        transient=afterglow,
        model='vegas_gaussian',  # Gaussian jet structure
        sampler='dynesty',
        nlive=2000,  # More live points for complex posteriors
        dlogz=0.1    # Tighter evidence tolerance
    )

    # Generate posterior predictions at new frequencies
    result.plot_multiband(frequencies=[1e9, 5e14, 1e18])

Example: Model Comparison
--------------------------

.. code-block:: python

    import redback

    afterglow = redback.afterglow.Afterglow.from_open_access_catalogue(
        'GRB070809', data_mode='flux_density'
    )

    # Fit multiple models
    models = ['vegas_tophat', 'vegas_gaussian', 'vegas_powerlaw']
    results = {}

    for model in models:
        results[model] = redback.fit_model(
            transient=afterglow,
            model=model,
            sampler='dynesty',
            nlive=1000
        )

    # Compare evidences
    for model, result in results.items():
        print(f"{model}: log(Z) = {result.log_evidence:.2f}")

    # Plot comparison
    redback.plot_model_comparison(results)

Advanced Features
-----------------

Reverse Shock Modeling
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import redback

    # Enable reverse shock with different microphysics
    kwargs = {
        'reverse_shock': True,
        'reverse_logepse': -1.5,  # Different from forward shock
        'reverse_logepsb': -2.5,
        'reverse_p': 2.2
    }

    result = redback.fit_model(
        transient=afterglow,
        model='vegas_tophat',
        model_kwargs=kwargs,
        sampler='dynesty',
        nlive=1000
    )

Jet Spreading
^^^^^^^^^^^^^

.. code-block:: python

    # Enable jet spreading dynamics
    kwargs = {'spread': True}

    result = redback.fit_model(
        transient=afterglow,
        model='vegas_tophat',
        model_kwargs=kwargs,
        sampler='dynesty',
        nlive=1000
    )

Best Practices
--------------

Start Simple
^^^^^^^^^^^^

Always begin with the simplest physically motivated model:

1. Start with ``vegas_tophat`` + ISM + forward shock only
2. Examine residuals and fit quality
3. Add complexity (reverse shock, different jet structure) only if justified

Choose Appropriate Priors
^^^^^^^^^^^^^^^^^^^^^^^^^^

* Use physical constraints to set prior ranges
* For viewing angle, use ``Sine`` prior (uniform in solid angle)
* For energies/densities, use ``Uniform`` in log space
* Check prior predictive distributions before fitting

Sampler Selection
^^^^^^^^^^^^^^^^^

* **dynesty**: Best for evidence calculation and multimodal posteriors
* **emcee**: Faster for unimodal posteriors, good for quick checks
* **ultranest**: Most robust for complex, multimodal posteriors

Computational Considerations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* VegasAfterglow is extremely fast (~1ms per light curve)
* Can fit 1000s of MCMC samples in minutes on laptop
* Use ``nlive=1000`` for quick runs, ``nlive=2000+`` for publication
* Parallel chains with ``npool`` parameter for faster sampling

Redback's Full Capabilities
----------------------------

Redback is a comprehensive transient analysis framework that provides everything you need for parameter estimation with VegasAfterglow models. Through the redback interface, you can:

**Data Management:**

* Load data from multiple sources: Swift, Fermi, BATSE, open catalogs, private data
* Support for all data types: flux density, flux, magnitude, luminosity
* Multi-wavelength data from radio to X-ray
* Automatic unit conversions and data preprocessing
* Filter and quality-cut data
* Combine datasets from multiple instruments
* Handle upper limits and non-detections

**Model Fitting:**

* Access all VegasAfterglow models (tophat, gaussian, power-law jets, etc.)
* Multiple sampling backends: dynesty, emcee, ultranest, nessai
* Custom prior specification with bilby
* Parallel sampling for faster computation
* Model comparison via Bayesian evidence
* Fit multiple transient types beyond afterglows

**Analysis & Visualization:**

* Corner plots for posterior distributions
* Light curve plots with model predictions
* Residual analysis
* Multi-band visualization
* Posterior predictive checks
* Evidence comparison tables
* Publication-ready figures

**Advanced Features:**

* Joint analysis of multiple transients
* Population inference
* Model selection and comparison
* Custom likelihood functions
* Systematic uncertainty modeling
* Time-dependent viewing angle effects
* Energy injection modeling

**Example: Complete Workflow**

.. code-block:: python

    import redback

    # 1. Get data (multiple methods available)
    GRB = 'GRB070809'

    # From Swift
    redback.get_data.get_bat_xrt_afterglow_data_from_swift(grb=GRB, data_mode="flux")
    afterglow = redback.afterglow.SGRB.from_swift_grb(name=GRB, data_mode='flux')

    # Or from your own data
    # import pandas as pd
    # data = pd.read_csv('my_grb_data.csv')
    # afterglow = redback.transient.Afterglow(
    #     name=GRB, data_mode='flux_density',
    #     time=data['time'].values, flux_density=data['flux'].values,
    #     flux_density_err=data['flux_err'].values, frequency=data['frequency'].values
    # )

    # 2. Fit model
    result = redback.fit_model(
        transient=afterglow,
        model='vegas_tophat',
        sampler='dynesty',
        nlive=2000
    )

    # 3. Analyze
    result.plot_corner()
    result.plot_lightcurve(save=True)
    result.print_summary()

    # 4. Make predictions
    result.plot_multiband(frequencies=[1e9, 1e14, 1e18])

For complete documentation on all these capabilities and more, see the comprehensive `redback documentation <https://redback.readthedocs.io/en/latest/?badge=latest>`_.

**Quick Links:**

* `Redback Installation <https://redback.readthedocs.io/en/latest/installation.html>`_
* `Data Loading Guide <https://redback.readthedocs.io/en/latest/data.html>`_
* `Fitting Tutorial <https://redback.readthedocs.io/en/latest/examples.html>`_
* `Model Library <https://redback.readthedocs.io/en/latest/models.html>`_
* `API Reference <https://redback.readthedocs.io/en/latest/api.html>`_

For VegasAfterglow physics details, see:

* :doc:`physics` - Detailed physics description
* :doc:`parameter_reference` - Complete parameter reference
* :doc:`examples` - VegasAfterglow usage examples
