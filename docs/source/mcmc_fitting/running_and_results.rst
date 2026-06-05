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

``FitResult`` also exposes two convenience methods:

- ``summary(top_k=None)``: returns a formatted top-K best-fit table (rank, chi-squared, parameter values), renders cleanly in both ``print()`` and Jupyter last-line auto-display.
- ``save(path)``: write the fit to disk in bilby's native HDF5 (``.h5`` / ``.hdf5``) or JSON (``.json``) format. The resulting file is interoperable with bilby's tooling. The original ``Fitter`` configuration (constructor args, observation data, parameter definitions) is persisted alongside the samples so a later session can reload the full state with :py:meth:`Fitter.load` -- no manual reconfiguration required.

To reload a saved fit, use :py:meth:`Fitter.load` (see below). Raw bilby Result files without VegasAfterglow metadata can be opened directly with ``bilby.read_in_result(path)``.

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

``FitResult.summary(top_k=None)`` returns a ranked table of the best-fit parameter sets, sorted by log-likelihood. Each row shows the rank, the chi-squared (``-2 * log_prob``; exact for the default Gaussian likelihood, a generic deviance for custom likelihoods), and the sampler-space value of every parameter (``Scale.log`` parameters appear as ``log10_<name>``):

.. code-block:: python

    # All top-K rows stored on the result (top_k was set at fit() time)
    print(result.summary())

    # Show only the top 3
    print(result.summary(top_k=3))

    # In a Jupyter cell, last-line auto-display also works:
    # result.summary()

Example output::

    Rank       chi^2     log10_E_iso       log10_Gamma0          theta_c               p
    -----------------------------------------------------------------------------------
       1       42.18         53.4127             2.4912           0.0612          2.2143
       2       43.05         53.3984             2.5031           0.0598          2.2255
       3       44.71         53.4276             2.4798           0.0625          2.2079

The returned object renders identically under ``print(result.summary())`` and a bare ``result.summary()`` in a Jupyter cell, so you can paste either pattern into a notebook.

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

.. _diagnostic-plot:

Diagnostic Plot
----------------

``Fitter.draw_fit()`` produces a publication-style two-panel diagnostic figure of the fitted model overlaid on the observation data, in one call:

.. code-block:: python

    fig, (ax_top, ax_bot) = fitter.draw_fit()

The figure layout:

- **Top panel** — observed data with errorbars plus the model curves. When posterior samples are available, each band's central line is the **posterior median** trajectory with a shaded ``ci`` credible band (default 68%); if no posterior is available the curve falls back to the MAP (``best_params``) trajectory. Single-frequency light curves (added via ``add_flux_density``) appear on the left y-axis as ``F_nu``; band-integrated fluxes (added via ``add_flux``) appear on the right y-axis as ``F``. When both kinds are present the two log y-axes are matched to span the same number of decades so a power-law slope alpha has identical visual slope on both sides. Multiple bands are auto-shifted vertically (rank-based, relaxed to avoid overlap) and the legend reports the shift factors.
- **Bottom panel** — observer-frame evolution of the synchrotron break frequencies ``nu_a``, ``nu_m``, ``nu_c``. Observed light-curve frequencies are overlaid as horizontal dashed lines; observed bands appear as shaded ``axhspan`` regions. Circles mark where each break-frequency curve crosses an observed line; squares mark crossings with the band edges -- useful for reading off when the model transitions through each observed band.

The break frequencies are evaluated at the jet ``theta`` column closest to the line of sight (``theta_v``), so the Doppler boost reflects what the off-axis observer sees (this collapses to the jet axis for ``theta_v == 0``).

Common kwargs:

- ``best_params=None`` — sampler-space parameter array. Defaults to ``self.result.top_k_params[0]``. Used for the break-frequency panel and as the model trajectory when posterior samples are unavailable or ``n_samples == 0``.
- ``ci=0.68`` — credible interval for the shaded band (e.g. ``0.95`` for ~2σ). Ignored when ``n_samples == 0``.
- ``n_samples=100`` — posterior draws used to construct the credible band. Set ``n_samples=0`` to skip the band entirely and plot the MAP trajectory only.
- ``t_range=None`` — ``(t_min, t_max)`` for the model curves in seconds. Default extends one decade below ``tmin`` and two decades above ``tmax`` of the observed data so the curves reach the visible edges.
- ``n_t=200`` — number of points on the model time grid.
- ``auto_shift_gap=1.0`` — decades added between consecutive bands (frequency-sorted). Each band is shifted by ``(rank - (n-1)/2) * auto_shift_gap``, then a second pass pushes bands apart locally if their shifted log-flux ranges still overlap. Raise this for wider visual separation.
- ``shifts=None`` — dict overriding individual auto-computed shifts. Keys are ``("lc", nu_hz)`` for light curves or ``("band", (nu_min, nu_max))`` for band-integrated entries; values are ``log10`` multipliers.
- ``show_nu_panel=True`` — set ``False`` to drop the bottom break-frequency panel and get a single-panel data/model plot.
- ``resolution=None`` — ``(phi, theta, t)`` override applied just to this rendering call.
- ``fig=None, axes=None`` — pass an existing figure / ``(ax_top, ax_bot)`` tuple to draw into instead of letting the method allocate its own.

Legend entries pick up filter / instrument names from the optional ``label=`` argument on ``add_flux_density`` and ``add_flux``:

.. code-block:: python

    fitter.add_flux_density(nu=4.6e14, t=t_r, f_nu=f_r, err=e_r, label="r")
    fitter.add_flux_density(nu=4.6e14, t=t_VT_R, f_nu=f_VT_R, err=e_VT_R, label="VT_R")
    fitter.add_flux(band=(0.5e3 * u_keV, 4.0e3 * u_keV), t=t_x, flux=f_x, err=e_x, label="WXT")

Known filter codes (SDSS ``u/g/r/i/z``, Johnson ``UBVRI``, 2MASS ``JHK``, X-ray instruments like ``WXT``, ``XRT``, ``BAT``) get curated colors from a built-in registry; same-color filters at different frequencies (e.g. ``r`` and ``VT_R``) are disambiguated via lightness so they remain visually distinct.

.. _posterior-credible-bands:

Posterior Credible Bands
-------------------------

For custom plotting (or any analysis that needs a flux uncertainty), ``Fitter.flux_density_credible()`` and ``Fitter.flux_credible()`` are the posterior-predictive siblings of ``flux_density_grid`` and ``flux``. They draw ``n_samples`` posterior samples, evaluate the model at each (in parallel), and return the central ``ci`` percentile envelope plus the posterior median trajectory.

.. code-block:: python

    # 68% credible band on multi-band light curves
    t_out = np.logspace(2, 9, 200)
    nu_out = np.array([1e9, 5e14, 2.4e17])           # radio, optical, X-ray [Hz]
    cb = fitter.flux_density_credible(t_out, nu_out, ci=0.68, n_samples=200)

    # Drop straight into matplotlib
    for i, nu in enumerate(nu_out):
        ax.plot(t_out, cb.median[i, :])
        ax.fill_between(t_out, cb.lower[i, :], cb.upper[i, :], alpha=0.2)

    # 95% credible band on band-integrated flux (e.g. XRT band-pass)
    cb_xrt = fitter.flux_credible(t_out, band=(2.4e17, 2.4e18), ci=0.95, n_samples=200)
    ax.fill_between(t_out, cb_xrt.lower, cb_xrt.upper, alpha=0.2)

Returned ``SimpleNamespace`` with three arrays:

- ``cb.lower`` / ``cb.upper`` — the ``(50 - 50·ci, 50 + 50·ci)`` percentiles at each ``(nu, t)`` cell. Shape is ``(len(nu), len(t))`` for ``flux_density_credible`` and ``(len(t),)`` for ``flux_credible``.
- ``cb.median`` — same shape, the posterior median trajectory. **This is what you usually want for the central line**: by construction it lies inside ``[lower, upper]`` at every cell. The MAP / best-fit trajectory does *not* — for correlated posteriors it can sit outside the 1-D marginal band at some points, which is statistically expected (the 16/84 envelope is a per-cell marginal, the MAP is a single joint point estimate).

Common kwargs:

- ``ci=0.68`` — credible interval (``0.95`` for ~2σ, ``0.99`` for ~3σ).
- ``n_samples=200`` — posterior draws. Larger gives smoother bands at the cost of more model evaluations.
- ``n_workers=None`` — thread parallelism. Defaults to ``os.cpu_count()``; pass an int to throttle.
- ``rng=None`` — ``numpy.random.Generator`` for reproducible draw selection.
- ``resolution=None`` — ``(phi, theta, t)`` override applied to all draws.

Extinction caveat: ``flux_density_credible`` *does* apply host-galaxy extinction per draw (matching ``flux_density_grid``), so the band reflects ``A_V`` uncertainty when ``A_V`` is a sampled parameter. ``flux_credible`` does *not* apply extinction — same as ``flux`` and the band-integrated chi² evaluation path.

``Fitter.draw_fit()`` uses these internally to render its credible band; call them directly when you want full control over the figure.

Saving and Loading Fits
------------------------

MCMC fits can take minutes to hours. To avoid re-running the sampler when reopening a notebook, persist the result to disk using bilby's native HDF5 / JSON format:

.. code-block:: python

    # After fitting -- symmetric with Fitter.load below.
    fitter.save("ep_grb_fit.h5")   # .h5 / .hdf5 (recommended) or .json

    # Later session: a single line restores the fitter (constructor args,
    # observation data, parameter transformer) and the result. No manual
    # reconfiguration, no MCMC re-run.
    from VegasAfterglow import Fitter
    fitter = Fitter.load("ep_grb_fit.h5")
    result = fitter.result

    print(result.summary())                            # works
    result.bilby_result.plot_corner()                  # bilby tools work
    lc = fitter.flux_density_grid(result.top_k_params[0], t_out, nu_out)

If the original fit used a custom callable for ``jet``, ``medium``, or ``extinction``, pass the same callable as a keyword to ``Fitter.load`` -- callables don't round-trip through HDF5/JSON serialization:

.. code-block:: python

    fitter = Fitter.load("ep_grb_fit.h5", jet=my_custom_jet_factory)

If you have a raw bilby Result file from another tool (no VegasAfterglow metadata), open it directly with ``bilby.read_in_result(path)`` instead -- ``Fitter.load`` will raise a clear ``ValueError`` for those files.

Files written by ``fitter.save()`` are full bilby Result files. They can also be opened directly by bilby's tools — useful when sharing fits with collaborators who use bilby:

.. code-block:: python

    import bilby
    br = bilby.read_in_result(filename="ep_grb_fit.h5")
    br.plot_corner()

Pre-existing bilby Result files from other tools (no VegasAfterglow metadata) can be opened directly with ``bilby.read_in_result(path)`` for samples / corner inspection — ``Fitter.load`` requires the VegasAfterglow snapshot so it can't load those files for predictions; re-run ``.fit(...)`` if you need them.

Visualization
--------------

For the data/model overlay use :ref:`Diagnostic Plot <diagnostic-plot>` above (``fitter.draw_fit()``); for custom plots with posterior-predictive uncertainty, build them yourself with :ref:`flux_density_credible <posterior-credible-bands>`. For posterior correlations, the standard tool is the external ``corner`` package:

.. code-block:: python

    import matplotlib.pyplot as plt
    import corner

    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])

    fig = corner.corner(
        flat_chain,
        labels=result.latex_labels,
        quantiles=[0.16, 0.5, 0.84],
        show_titles=True,
        title_kwargs={"fontsize": 12},
    )
    plt.savefig("corner_plot.png", dpi=300, bbox_inches="tight")

bilby's own corner helper is also available via ``result.bilby_result.plot_corner()`` if you prefer the bilby styling.
