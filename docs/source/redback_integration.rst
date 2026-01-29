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

Basic Parameter Estimation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import redback

    # Load GRB data from open access catalogs
    afterglow = redback.afterglow.Afterglow.from_open_access_catalogue(
        'GRB170817A', data_mode='flux_density'
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
        'GRB170817A', data_mode='flux_density'
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
        'GRB170817A', data_mode='flux_density'
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
        'GRB170817A', data_mode='flux_density'
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

Further Documentation
---------------------

For complete redback documentation, including:

* Data loading and management
* Custom model creation
* Result visualization
* Advanced sampling techniques

See the `redback documentation <https://redback.readthedocs.io/>`_.

For VegasAfterglow physics details, see:

* :doc:`physics` - Detailed physics description
* :doc:`parameter_reference` - Complete parameter reference
* :doc:`examples` - VegasAfterglow usage examples
