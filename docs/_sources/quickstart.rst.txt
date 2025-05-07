Quickstart
==========

This guide will help you get started with VegasAfterglow quickly. We'll cover basic installation, setting up a simple model, and running your first afterglow simulation.

Installation
-----------

The easiest way to install VegasAfterglow is via pip:

.. code-block:: bash

    pip install VegasAfterglow

For more detailed installation instructions, see the :doc:`installation` page.

Basic Usage
----------

First, let's import the necessary modules:

.. code-block:: python

    import numpy as np
    import matplotlib.pyplot as plt
    from VegasAfterglow import Model, TophatJet, ISM, Synchrotron

Setting up a simple model
^^^^^^^^^^^^^^^^^^^^^^^

To create a basic GRB afterglow model, we need to define three core components:

1. The **Jet**: Describes the relativistic outflow from the central engine
2. The **Ambient Medium**: Defines the environment surrounding the GRB
3. The **Radiation Process**: Specifies the emission mechanism

.. code-block:: python

    # Create a simple on-axis top-hat jet
    # Parameters: theta_c (core angle), E_iso (isotropic energy), Gamma0 (Lorentz factor)
    jet = TophatJet(
        theta_c=0.1,        # Core angle in radians
        E_iso=1e53,         # Isotropic-equivalent energy in ergs
        Gamma0=300          # Initial bulk Lorentz factor
    )
    
    # Create a uniform ISM medium
    # Parameter: n_ism (number density in cm^-3)
    medium = ISM(n_ism=1.0)
    
    # Create a synchrotron radiation process
    # Parameters: epsilon_e, epsilon_B (energy fractions), p (electron index)
    radiation = Synchrotron(
        epsilon_e=0.1,      # Fraction of energy in electrons
        epsilon_B=0.01,     # Fraction of energy in magnetic field
        p=2.2               # Electron energy index
    )
    
    # Create the model
    model = Model(jet=jet, medium=medium, radiation=radiation)

Calculating light curves
^^^^^^^^^^^^^^^^^^^^^

Once we have set up our model, we can calculate light curves at different frequencies:

.. code-block:: python

    # Define observation times (in seconds)
    times = np.logspace(2, 7, 100)  # 100 seconds to 10^7 seconds
    
    # Define observation frequencies (in Hz)
    frequencies = np.array([
        1e9,   # 1 GHz (radio)
        5e14,  # Optical R-band
        1e17   # X-ray (~0.4 keV)
    ])
    
    # Calculate light curves
    results = model.calculate_light_curves(times, frequencies)

Plotting the results
^^^^^^^^^^^^^^^^^

Let's visualize the light curves:

.. code-block:: python

    plt.figure(figsize=(10, 6))
    
    labels = ['Radio (1 GHz)', 'Optical (R-band)', 'X-ray (0.4 keV)']
    
    for i, nu in enumerate(frequencies):
        plt.loglog(times, results[:, i], label=labels[i])
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('GRB Afterglow Light Curves')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Using a structured jet model
^^^^^^^^^^^^^^^^^^^^^^^^

Let's try a more complex model with a Gaussian structured jet:

.. code-block:: python

    # Create a Gaussian structured jet
    # Parameters: theta_c (core angle), E_iso (isotropic energy), Gamma0 (Lorentz factor)
    gaussian_jet = GaussianJet(
        theta_c=0.05,        # Core angle in radians
        E_iso=1e53,          # Isotropic-equivalent energy in ergs
        Gamma0=300           # Initial bulk Lorentz factor
    )
    
    # Update the model with the new jet
    model.set_jet(gaussian_jet)
    
    # Recalculate with the structured jet
    results_gaussian = model.calculate_light_curves(times, frequencies)
    
    # Compare with the original top-hat jet
    plt.figure(figsize=(10, 6))
    
    for i, nu in enumerate(frequencies):
        plt.loglog(times, results[:, i], '--', color=f'C{i}', label=f'TophatJet {labels[i]}')
        plt.loglog(times, results_gaussian[:, i], '-', color=f'C{i}', label=f'GaussianJet {labels[i]}')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('Comparison of Jet Models')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Using a wind medium
^^^^^^^^^^^^^^^

Let's change the ambient medium to a stellar wind profile:

.. code-block:: python

    # Create a stellar wind medium
    # Parameter: A_star (wind density parameter)
    wind = Wind(A_star=0.1)
    
    # Update the model
    model.set_medium(wind)
    
    # Calculate new results with wind medium
    results_wind = model.calculate_light_curves(times, frequencies)
    
    # Plot and compare
    plt.figure(figsize=(10, 6))
    
    for i, nu in enumerate(frequencies):
        plt.loglog(times, results_gaussian[:, i], '--', color=f'C{i}', label=f'ISM {labels[i]}')
        plt.loglog(times, results_wind[:, i], '-', color=f'C{i}', label=f'Wind {labels[i]}')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('Comparison of Ambient Media')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Computing a spectrum
^^^^^^^^^^^^^^^^^

We can also calculate the spectrum at a specific time:

.. code-block:: python

    # Define a wide range of frequencies
    frequencies_spectrum = np.logspace(9, 19, 100)  # 1 GHz to 10^19 Hz
    
    # Calculate spectrum at 1 day post-burst
    t_day = 86400  # seconds in a day
    spectrum = model.calculate_spectrum(t_day, frequencies_spectrum)
    
    # Plot the spectrum
    plt.figure(figsize=(10, 6))
    plt.loglog(frequencies_spectrum, spectrum)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.title(f'GRB Afterglow Spectrum at t = {t_day/86400:.1f} days')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()

Next Steps
---------

Now that you've created your first models, here are some ways to extend them:

- Try different parameter values to see how they affect the afterglow evolution
- Explore off-axis viewing with `GaussianJet(theta_c=0.05, E_iso=1e53, Gamma0=300, viewer_angle=0.2)`
- Use the `PowerLawJet` class to model jets with power-law profiles
- Add inverse Compton scattering with the `SynchrotronSelfCompton` class
- Compare models with observational data using the `Fitter` class

For more detailed examples, see the :doc:`examples` page. 