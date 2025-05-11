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

VegasAfterglow is designed to efficiently model gamma-ray burst (GRB) afterglows and perform Markov Chain Monte Carlo (MCMC) parameter estimation. 

Direct Model Calculation
------------------------
Before diving into MCMC parameter estimation, you can directly use VegasAfterglow to generate light curves and spectra from a specific model. Let's start by importing the necessary modules:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model


Then, let's set up the physical components of our afterglow model, including the environment, jet, observer, and radiation parameters:

.. code-block:: python

    # 1. Define the circumburst environment (constant density ISM)
    medium = ISM(n_ism=1)

    # 2. Configure the jet structure (top-hat with opening angle, energy, and Lorentz factor)
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)

    # 3. Set observer parameters (distance, redshift, viewing angle)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)

    # 4. Define radiation microphysics parameters
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)

    # 5. Combine all components into a complete afterglow model
    model = Model(jet=jet, medium=medium, observer=obs, forward_rad=rad)

Light Curve Calculation
^^^^^^^^^^^^^^^^^^^^^^^

Now, let's compute and plot multi-wavelength light curves to see how the afterglow evolves over time:

.. code-block:: python

    # 1. Create logarithmic time array from 10² to 10⁸ seconds (100s to ~3yrs)
    times = np.logspace(2, 8, 200)  

    # 2. Define observing frequencies (radio, optical, X-ray bands in Hz)
    bands = np.array([1e9, 1e14, 1e17])  

    # 3. Calculate the afterglow emission at each time and frequency
    results = model.specific_flux(times, bands)

    # 4. Visualize the multi-wavelength light curves
    plt.figure(figsize=(4.8, 3.6), dpi=200)

    # 5. Plot each frequency band 
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

    # 6. Add annotations for important transitions
    plt.annotate('jet break', xy=(3e4, 1e-26), xytext=(3e3, 5e-28), 
                arrowprops=dict(arrowstyle='->'))
    plt.annotate(r'$\nu_m=\nu_a$', xy=(3e5, 3e-25), xytext=(5.5e4, 5e-24), 
                arrowprops=dict(arrowstyle='->'))
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('Light Curves')
    plt.tight_layout()
    plt.savefig('assets/quick-lc.png', dpi=300)

.. figure:: /_static/images/quick-lc.png
   :width: 600
   :align: center
   
   Running the light curve script will produce this figure showing the afterglow evolution across different frequencies.

Spectral Analysis
^^^^^^^^^^^^^^^^^

We can also examine how the broadband spectrum evolves at different times after the burst:

.. code-block:: python

    # 1. Define broad frequency range (10⁵ to 10²² Hz) 
    frequencies = np.logspace(5, 22, 200)  

    # 2. Select specific time epochs for spectral snapshots 
    epochs = np.array([1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8])

    # 3. Calculate spectra at each epoch
    results = model.spectra(frequencies, epochs)

    # 4. Plot broadband spectra at each epoch
    plt.figure(figsize=(4.8, 3.6), dpi=200)
    colors = plt.cm.viridis(np.linspace(0, 1, len(epochs)))

    for i, t in enumerate(epochs):
        exp = int(np.floor(np.log10(t)))
        base = t / 10**exp
        plt.loglog(frequencies, results['syn'][i,:], color=colors[i], 
                  label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

    # 5. Add vertical lines marking the bands from the light curve plot
    for i, band in enumerate(bands):
        plt.axvline(band, ls='--', color=f'C{i}')

    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend(ncol=2)
    plt.title('Synchrotron Spectra')
    plt.tight_layout()
    plt.savefig('assets/quick-spec.png', dpi=300)

.. figure:: /_static/images/quick-spec.png
   :width: 600
   :align: center
   
   The spectral analysis code will generate this visualization showing spectra at different times, with vertical lines indicating the frequencies calculated in the light curve example.

Parameter Estimation with MCMC
------------------------------

For more advanced analysis, VegasAfterglow provides powerful MCMC capabilities to fit model parameters to observational data. 

First, let's import the necessary modules:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    import pandas as pd
    import corner
    from VegasAfterglow import ObsData, Setups, Fitter, ParamDef, Scale

Preparing Data and Configuring the Model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

VegasAfterglow provides flexible options for loading observational data through the ``ObsData`` class. You can add light curves (specific flux vs. time) and spectra (specific flux vs. frequency) in multiple ways:

.. code-block:: python

    # Create an instance to store observational data
    data = ObsData()

    # Method 1: Add data directly from lists or numpy arrays
    
    # For light curves
    t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
    flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # Specific flux in erg/cm²/s/Hz
    flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]  # Specific flux error in erg/cm²/s/Hz
    data.add_light_curve(nu_cgs=4.84e14, t_cgs=t_data, Fnu_cgs=flux_data, Fnu_err=flux_err)

    # For spectra
    nu_data = [...]  # Frequencies in Hz
    spectrum_data = [...] # Specific flux values in erg/cm²/s/Hz
    spectrum_err = [...]   # Specific flux errors in erg/cm²/s/Hz
    data.add_spectrum(t_cgs=3000, nu_cgs=nu_data, Fnu_cgs=spectrum_data, Fnu_err=spectrum_err)

.. code-block:: python

    # Method 2: Load from CSV files
    data = ObsData()
    # Define your bands and files
    bands = [2.4e17, 4.84e14, 1.4e14]  # Example: X-ray, optical R-band
    lc_files = ["data/ep.csv", "data/r.csv", "data/vt-r.csv"]

    # Load light curves from files
    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        data.add_light_curve(nu_cgs=nu, t_cgs=df["t"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

    times = [3000] # Example: time in seconds
    spec_files = ["data/ep-spec.csv"]

    # Load spectra from files
    for t, fname in zip(times, spec_files):
        df = pd.read_csv(fname)
        data.add_spectrum(t_cgs=t, nu_cgs=df["nu"], Fnu_cgs=df["Fv_obs"], Fnu_err=df["Fv_err"])

.. note::
   The ``ObsData`` interface is designed to be flexible. You can mix and match different data sources, and add multiple light curves at different frequencies as well as multiple spectra at different times.

The ``Setups`` class defines the global properties and environment for your model. These settings remain fixed during the MCMC process:

.. code-block:: python

    cfg = Setups()

    # Source properties
    cfg.lumi_dist = 3.364e28    # Luminosity distance [cm]  
    cfg.z = 1.58               # Redshift

    # Physical model configuration
    cfg.medium = "wind"        # Ambient medium: "wind", "ISM" (Interstellar Medium) or "user" (user-defined)
    cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat" or "user" (user-defined)


These settings affect how the model is calculated but are not varied during the MCMC process.

Defining Parameters and Running MCMC
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``ParamDef`` class is used to define the parameters for MCMC exploration. Each parameter requires a name, initial value, prior range, and sampling scale:

.. code-block:: python

    mc_params = [
        ParamDef("E_iso",    1e52,  1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
        ParamDef("Gamma0",     30,     5,  1000,  Scale.LOG),       # Lorentz factor at the core
        ParamDef("theta_c",   0.2,   0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
        ParamDef("theta_v",    0.,  None,  None,  Scale.FIXED),     # Viewing angle [rad]
        ParamDef("p",         2.5,     2,     3,  Scale.LINEAR),    # Shocked electron power law index
        ParamDef("eps_e",     0.1,  1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
        ParamDef("eps_B",    1e-2,  1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
        ParamDef("A_star",   0.01,  1e-3,     1,  Scale.LOG),       # Wind parameter
        ParamDef("xi",        0.5,  1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
    ]

**Scale Types:**
    - ``Scale.LOG``: Sample in logarithmic space (log10) - ideal for parameters spanning multiple orders of magnitude
    - ``Scale.LINEAR``: Sample in linear space - appropriate for parameters with narrower ranges
    - ``Scale.FIXED``: Keep parameter fixed at the initial value - use for parameters you don't want to vary

**Parameter Choices:**
The parameters you include depend on your model configuration:
    - For "wind" medium: use ``A_star`` parameter 
    - For "ISM" medium: use ``n_ism`` parameter instead
    - Different jet structures may require different parameters

Initialize the ``Fitter`` class with your data and configuration, then run the MCMC process:

.. code-block:: python

    # Create the fitter object
    fitter = Fitter(data, cfg)

    # Run the MCMC fitting
    result = fitter.fit(
        param_defs=mc_params,          # Parameter definitions
        resolution=(0.1, 1, 5),        # Grid resolution (see more details in `Examples`)
        total_steps=10000,             # Total number of MCMC steps
        burn_frac=0.3,                 # Fraction of steps to discard as burn-in
        thin=1                         # Thinning factor
    )

The ``result`` object contains:
    - ``samples``: The MCMC chain samples (posterior distribution)
    - ``labels``: Parameter names
    - ``best_params``: Maximum likelihood parameter values

Analyzing Results and Generating Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Check the best-fit parameters and their uncertainties:

.. code-block:: python

    # Print best-fit parameters (maximum likelihood)
    print("Best-fit parameters:")
    for name, val in zip(result.labels, result.best_params):
        print(f"  {name}: {val:.4f}")

    # Compute median and credible intervals
    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
    medians = np.median(flat_chain, axis=0)
    lower = np.percentile(flat_chain, 16, axis=0)
    upper = np.percentile(flat_chain, 84, axis=0)

    print("\nParameter constraints (median and 68% credible intervals):")
    for i, name in enumerate(result.labels):
        print(f"  {name}: {medians[i]:.4f} (+{upper[i]-medians[i]:.4f}, -{medians[i]-lower[i]:.4f})")

Use the best-fit parameters to generate model predictions:

.. code-block:: python

    # Define time and frequency ranges for predictions
    t_out = np.logspace(2, 9, 150)
    bands = [2.4e17, 4.84e14, 1.4e14] 

    # Generate light curves with the best-fit model
    lc_best = fitter.light_curves(result.best_params, t_out, bands)

    nu_out = np.logspace(6, 20, 150)
    times = [3000]
    # Generate model spectra at the specified times using the best-fit parameters
    spec_best = fitter.spectra(result.best_params, nu_out, times)

Now you can plot the best-fit model:

.. code-block:: python

    def draw_bestfit(t, lc_fit, nu, spec_fit):
        # Create figure with two subplots
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(4.5, 7.5))
        
        # Plot light curves
        shifts = [1, 1, 200]
        colors = ['blue', 'orange', 'green']
        
        for i in range(len(lc_files)):
            df = pd.read_csv(lc_files[i])
            ax1.errorbar(df["t"], df["Fv_obs"] * shifts[i], df["Fv_err"] * shifts[i], 
                        fmt='o', color=colors[i], label=lc_files[i])
            ax1.plot(t, np.array(lc_fit[i]) * shifts[i], color=colors[i], lw=1)

        # Plot spectra
        for i in range(len(spec_files)):
            df = pd.read_csv(spec_files[i])
            ax2.errorbar(df["nu"], df["Fv_obs"] * shifts[i], df["Fv_err"] * shifts[i], 
                        fmt='o', color=colors[i], label=spec_files[i])
            ax2.plot(nu, np.array(spec_fit[0]) * shifts[i], color=colors[i], lw=1)

        # Configure axes
        for ax, xlabel, ylabel in [(ax1, 't [s]', r'$F_\nu$ [erg/cm$^2$/s/Hz]'),
                                  (ax2, r'$\nu$ [Hz]', r'$F_\nu$ [erg/cm$^2$/s/Hz]')]:
            ax.set_xscale('log'); ax.set_yscale('log')
            ax.set_xlabel(xlabel); ax.set_ylabel(ylabel)
            ax.legend()

        plt.tight_layout()

    draw_bestfit(t_out, lc_best, nu_out, spec_best)

Corner plots are essential for visualizing parameter correlations and posterior distributions:

.. code-block:: python

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
            smooth=1,
            fill_contours=True,
            levels=[0.16, 0.5, 0.68],  # 1σ and 2σ contours
            color='k'
        )
        fig.savefig(filename, dpi=300, bbox_inches='tight')

    # Create the corner plot
    flat_chain = result.samples.reshape(-1, result.samples.shape[-1])
    plot_corner(flat_chain, result.labels)

Next Steps
----------

See the :doc:`examples` page for more detailed examples

