Running MCMC
============

Choosing a Sampler
-------------------

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
-----------------------

The ``Fitter.fit()`` method provides a unified interface for parameter estimation using different sampling algorithms via bilby.

**Method Signature**

.. code-block:: python

    # noqa: syntax
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
        log_likelihood_fn: Callable = None,  # Custom log-likelihood
        priors: dict = None,                 # Custom bilby priors (all samplers)
        **sampler_kwargs                     # Sampler-specific parameters
    ) -> FitResult

**Common Parameters**

- ``param_defs``: List of ``ParamDef`` objects defining parameters to fit
    - Required parameter for all samplers
    - See "Parameter Definition" section for details

- ``resolution``: Optional tuple of (phi_res, theta_res, t_res)
    - Overrides the resolution set in the ``Fitter`` constructor for this run
    - If ``None``, uses the constructor value (default: ``(0.1, 0.25, 10)``)
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
        resolution=(0.1, 0.25, 10),
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
---------------------

.. code-block:: python

    # Create fitter object and add data (see above for data loading examples)
    fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")
    # ... add data via fitter.add_flux_density(), fitter.add_spectrum(), etc.

    # Option 1 (Recommended): Nested sampling with dynesty (computes Bayesian evidence, robust for multimodal posteriors)
    result = fitter.fit(
        params,
        resolution=(0.1, 0.25, 10),      # Grid resolution (phi, theta, t)
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
        resolution=(0.1, 0.25, 10),      # Grid resolution (phi, theta, t)
        sampler="emcee",               # MCMC sampler
        nsteps=50000,                  # Number of steps per walker
        nburn=10000,                   # Burn-in steps to discard
        thin=1,                        # Save every nth sample
        npool=8,                       # Number of parallel threads
        top_k=10,                      # Number of best-fit parameters to return
    )

**Important Notes:**
    - Parameters with ``Scale.log`` are sampled as ``log10_<name>`` (e.g., ``log10_E_iso``)
    - The sampler works in log10 space for LOG-scale parameters, then transforms back to physical values
    - Use ``npool`` to parallelize likelihood evaluations across multiple CPU cores
    - ``result.latex_labels`` provides properly formatted labels for corner plots
    - For emcee, total samples = ``nwalkers * (nsteps - nburn) / thin``
    - For dynesty, sample count depends on convergence (not fixed)


Analyzing Results
==================

Parameter Constraints
----------------------

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
------------------

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
    from VegasAfterglow.units import band

    flux_integrated = fitter.flux(result.top_k_params[0], t_model,
                                  band=band("BAT"))

Visualization
--------------

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
