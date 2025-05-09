Quickstart
==========

This guide will help you get started with VegasAfterglow quickly. We'll cover basic installation, setting up a simple model, and running your first afterglow parameter estimation.

Installation
------------

The easiest way to install VegasAfterglow is via pip:

.. code-block:: bash

    pip install VegasAfterglow

For more detailed installation instructions, see the :doc:`installation` page.

Basic Usage
-----------

VegasAfterglow is designed to efficiently model gamma-ray burst (GRB) afterglows and perform Markov Chain Monte Carlo (MCMC) parameter estimation. Let's start by importing the necessary modules:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale

Preparing Observational Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, we need to create an instance to store observational data through the ``ObsData`` class. You can add light curves (specific flux vs. time) and spectra (specific flux vs. frequency) in multiple ways:

.. code-block:: python

    # Create an instance to store observational data
    data = ObsData()

    # Method 1: Add data directly from lists or numpy arrays
    
    # For light curves
    t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
    flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # Specific flux in erg/cm²/s/Hz
    flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]  # Specific flux error
    data.add_light_curve(nu_cgs=4.84e14, t_cgs=t_data, Fnu_cgs=flux_data, Fnu_err=flux_err)
    
    # For spectra
    nu_data = [1e14, 3e14, 5e14, 1e15]  # Frequencies in Hz
    spectrum_data = [2e-27, 1.8e-27, 1.6e-27, 1.2e-27]  # Specific flux values
    spectrum_err = [2e-29, 1.8e-29, 1.6e-29, 1.2e-29]  # Specific flux errors
    data.add_spectrum(t_cgs=3000, nu_cgs=nu_data, Fnu_cgs=spectrum_data, Fnu_err=spectrum_err)

    # Method 2: Load from CSV files
    import pandas as pd
    
    # Define your bands and files
    bands = [2.4e17, 4.84e14]  # Example: X-ray, optical R-band
    lc_files = ["data/x-ray.csv", "data/r-band.csv"]
    
    # Load light curves from files
    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

Setting up the Model Configuration
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``Setups`` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process:

.. code-block:: python

    cfg = Setups()
    
    # Source properties
    cfg.lumi_dist = 3.364e28  # Luminosity distance [cm]  
    cfg.z = 1.58              # Redshift
    
    # Physical model configuration
    cfg.medium = "wind"       # Ambient medium: "wind", "ISM" or "user"
    cfg.jet = "powerlaw"      # Jet structure: "powerlaw", "gaussian", "tophat" or "user"
    
    # Optional: Advanced grid settings (default is 24, 24, 24)
    # cfg.phi_num = 24        # Number of grid points in phi direction
    # cfg.theta_num = 24      # Number of grid points in theta direction
    # cfg.t_num = 24          # Number of time grid points

Defining MCMC Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

The ``ParamDef`` class is used to define the parameters for MCMC exploration. Each parameter requires a name, initial value, prior range, and sampling scale:

.. code-block:: python

    mc_params = [
        ParamDef("E_iso",    1e52,  1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
        ParamDef("Gamma0",     30,     5,  1000,  Scale.LOG),       # Lorentz factor at the core
        ParamDef("theta_c",   0.2,   0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
        ParamDef("theta_v",   0.,  None,  None,   Scale.FIXED),     # Viewing angle [rad]
        ParamDef("p",         2.5,     2,     3,  Scale.LINEAR),    # Shocked electron power law index
        ParamDef("eps_e",     0.1,  1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
        ParamDef("eps_B",    1e-2,  1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
        ParamDef("A_star",   0.01,  1e-3,     1,  Scale.LOG),       # Wind parameter
        ParamDef("xi",        0.5,  1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
    ]

The parameters you include depend on your model configuration:
- For "wind" medium: use ``A_star`` parameter 
- For "ISM" medium: use ``n_ism`` parameter instead
- Different jet structures may require different parameters

Scale types:
- ``Scale.LOG``: Sample in logarithmic space (log10) - ideal for parameters spanning multiple orders of magnitude
- ``Scale.LINEAR``: Sample in linear space - appropriate for parameters with narrower ranges
- ``Scale.FIXED``: Keep parameter fixed at the initial value - use for parameters you don't want to vary

Running the MCMC
^^^^^^^^^^^^^^^^

Initialize the ``Fitter`` class with your data and configuration, then run the MCMC process:

.. code-block:: python

    # Create the fitter object
    fitter = Fitter(data, cfg)
    
    # Run the MCMC fitting
    result = fitter.fit(
        param_defs=mc_params,          # Parameter definitions
        resolution=(24, 24, 24),       # Grid resolution (phi, theta, time)
        total_steps=10000,             # Total number of MCMC steps
        burn_frac=0.3,                 # Fraction of steps to discard as burn-in
        thin=1                         # Thinning factor
    )

Analyzing the Results
^^^^^^^^^^^^^^^^^^^^^^

Examine the posterior distribution to understand parameter constraints:

.. code-block:: python

    # Print best-fit parameters (maximum likelihood)
    print("Best-fit parameters:")
    for name, val in zip(result.labels, result.best_params):
        print(f"  {name}: {val:.4g}")
    
    # Compute median and credible intervals
    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
    medians = np.median(flat_chain, axis=0)
    lower = np.percentile(flat_chain, 16, axis=0)
    upper = np.percentile(flat_chain, 84, axis=0)
    
    print("\nParameter constraints (median and 68% credible intervals):")
    for i, name in enumerate(result.labels):
        print(f"  {name}: {medians[i]:.4g} (+{upper[i]-medians[i]:.4g}, -{medians[i]-lower[i]:.4g})")

Generating Model Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use samples from the posterior to generate model predictions with uncertainties:

.. code-block:: python

    # Define time and frequency ranges for predictions
    t_out = np.logspace(2, 9, 150)
    bands = [2.4e17, 4.84e14]  # X-ray and optical R-band
    
    # Generate light curves with the best-fit model
    lc_best = fitter.light_curves(result.best_params, t_out, bands)
    
    # Generate spectra at specific times
    nu_out = np.logspace(6, 20, 150)
    times = [3000]  # 3000 seconds
    spec_best = fitter.spectra(result.best_params, nu_out, times)
    
    # Plot the light curves
    plt.figure(figsize=(10, 6))
    for i, nu in enumerate(bands):
        plt.loglog(t_out, lc_best[:, i], label=f'ν = {nu:.2e} Hz')
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Visualizing Parameter Correlations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Corner plots are essential for visualizing parameter correlations and posterior distributions:

.. code-block:: python

    import corner
    
    def plot_corner(flat_chain, labels, filename="corner_plot.png"):
        fig = corner.corner(
            flat_chain,
            labels=labels,
            quantiles=[0.16, 0.5, 0.84],  # For median and ±1σ
            show_titles=True,
            title_kwargs={"fontsize": 14},
            label_kwargs={"fontsize": 14},
            truths=np.median(flat_chain, axis=0),  # Show median values
            truth_color='red',
            bins=30,
            fill_contours=True,
            levels=[0.16, 0.5, 0.68],  # 1σ and 2σ contours
            color='k'
        )
        fig.savefig(filename, dpi=300, bbox_inches='tight')
    
    # Create the corner plot
    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
    plot_corner(flat_chain, result.labels)

Checking MCMC Convergence
^^^^^^^^^^^^^^^^^^^^^^^^^

Trace plots help verify MCMC convergence:

.. code-block:: python

    def plot_trace(chain, labels, filename="trace_plot.png"):
        nsteps, nwalkers, ndim = chain.shape
        fig, axes = plt.subplots(ndim, figsize=(10, 2.5 * ndim), sharex=True)
        
        for i in range(ndim):
            for j in range(nwalkers):
                axes[i].plot(chain[:, j, i], alpha=0.5, lw=0.5)
            axes[i].set_ylabel(labels[i])
            
        axes[-1].set_xlabel("Step")
        plt.tight_layout()
        plt.savefig(filename, dpi=300)
    
    # Create the trace plot
    plot_trace(result.samples, result.labels)

Next Steps
----------

Now that you've run your first MCMC parameter estimation, here are some suggested next steps:

1. **Prior Ranges**: Adjust the prior ranges based on theoretical constraints for your GRB event
2. **Convergence Testing**: Experiment with different numbers of steps and check convergence metrics
3. **Model Comparison**: Try different physical models (e.g., wind vs. ISM medium, or different jet structures)
4. **Physical Interpretation**: Connect your parameter constraints with physical processes in GRB afterglows
5. **Examine Examples**: See the :doc:`examples` page for more detailed examples

For more advanced users, you can also check the C++ interface for creating custom problem generators (see :doc:`cpp_api`). 