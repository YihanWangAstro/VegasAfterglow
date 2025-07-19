Examples
========

.. contents:: Table of Contents
   :local:
   :depth: 2

Basic Usage
-----------

Setting up a simple afterglow model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    import matplotlib.pyplot as plt
    import numpy as np
    from VegasAfterglow import ISM, TophatJet, Observer, Radiation, Model
    
    # Define the circumburst environment (constant density ISM)
    medium = ISM(n_ism=1)

    # Configure the jet structure (top-hat with opening angle, energy, and Lorentz factor)
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300)

    # Set observer parameters (distance, redshift, viewing angle)
    obs = Observer(lumi_dist=1e26, z=0.1, theta_obs=0)

    # Define radiation microphysics parameters
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, xi_e=1)

    # Combine all components into a complete afterglow model
    model = Model(jet=jet, medium=medium, observer=obs, forward_rad=rad, resolutions=(0.3,3,5))

    # Define time range for light curve calculation
    times = np.logspace(2, 8, 200)  

    # Define observing frequencies (radio, optical, X-ray bands in Hz)
    bands = np.array([1e9, 1e14, 1e17])  

    # Calculate the afterglow emission at each time and frequency
    results = model.specific_flux(times, bands)

    # Visualize the multi-wavelength light curves
    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band 
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()

    # Define broad frequency range (10⁵ to 10²² Hz) 
    frequencies = np.logspace(5, 22, 200)  

    # Select specific time epochs for spectral snapshots 
    epochs = np.array([1e2, 1e3, 1e4, 1e5 ,1e6, 1e7, 1e8])

    # Calculate spectra at each epoch
    results = model.specific_flux(epochs, frequencies)

    # Plot broadband spectra at each epoch
    plt.figure(figsize=(4.8, 3.6),dpi=200)
    colors = plt.cm.viridis(np.linspace(0,1,len(epochs)))

    for i, t in enumerate(epochs):
        exp = int(np.floor(np.log10(t)))
        base = t / 10**exp
        plt.loglog(frequencies, results['syn'][:,i], color=colors[i], label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

    # Add vertical lines marking the bands from the light curve plot
    for i, band in enumerate(bands):
        exp = int(np.floor(np.log10(band)))
        base = band / 10**exp
        plt.axvline(band,ls='--',color='C'+str(i))

    plt.xlabel('frequency (Hz)')
    plt.ylabel('flux density (erg/cm²/s/Hz)')
    plt.legend(ncol=2)
    plt.title('Synchrotron Spectra')


Ambient Media Models
--------------------

Wind Medium
^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Wind

    # Create a stellar wind medium
    wind = Wind(A_star=0.1)  # A* parameter

    #..other settings
    model = Model(medium=wind, ...)

User-Defined Medium
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Medium

    m_p = 1.6726219e-24  # Proton mass in grams

    # Define a custom density profile function
    def density(phi, theta, r):# r in cm, phi and theta in radians
        return m_p # n_ism =  1 cm^-3
        #return whatever density profile (g/cm^3) you want as a function of phi, theta, and r
        
    def swept_mass(phi, theta, r):# r in cm, phi and theta in radians
        return m_p * r * r * r / 3
        #return the swept up mass PER SOLID ANGLE (g) over r (\int_0^r rho*r*r dr).
        #you may keep the consistency of the mass profile with the density profile
        #the purpose of providing the extra mass profile is to reduce the extra computations.
    
    # Create a user-defined medium
    medium = Medium(rho=density, mass=swept_mass)

    #..other settings
    model = Model(medium=medium, ...)


Jet Models
----------

Gaussian Jet
^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import GaussianJet

    # Create a structured jet with Gaussian energy profile
    jet = GaussianJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300            # Initial Lorentz factor
    )

    #..other settings
    model = Model(jet=jet, ...)

Power-Law Jet
^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import PowerLawJet

    # Create a power-law structured jet
    jet = PowerLawJet(
        theta_c=0.05,         # Core angular size (radians)
        E_iso=1e53,           # Isotropic-equivalent energy (ergs)
        Gamma0=300,           # Initial Lorentz factor
        k=2.0                 # Power-law index
    )

    #..other settings
    model = Model(jet=jet, ...)

Two-Component Jet
^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import TwoComponentJet

    # Create a two-component jet
    jet = TwoComponentJet(
        theta_n=0.05,        # Narrow component angular size (radians)
        E_iso_n=1e53,        # Isotropic-equivalent energy of the narrow component (ergs)
        Gamma0_n=300,        # Initial Lorentz factor of the narrow component
        theta_w=0.1,         # Wide component angular size (radians)
        E_iso_w=1e52,        # Isotropic-equivalent energy of the wide component (ergs)
        Gamma0_w=100         # Initial Lorentz factor of the wide component
    )

    #..other settings
    model = Model(jet=jet, ...)

Jet with Spreading
^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import TophatJet

    jet = TophatJet(
        theta_c=0.05,        
        E_iso=1e53,          
        Gamma0=300,         
        spreading=True       # Enable spreading
    )

    #..other settings
    model = Model(jet=jet, ...)

.. note::
    The jet spreading (Lateral Expansion) is experimental and only works for the top-hat jet, Gaussian jet, and power-law jet with a jet core.
    The spreading prescription may not work for arbitrary user-defined jet structures.

Magnetar Spin-down
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Magnetar

    # Create a tophat jet with magnetar spin-down energy injection; Luminosity 1e46 erg/s, t_0 = 100 seconds, and q = 2
    jet = TophatJet(theta_c=0.05, E_iso=1e53, Gamma0=300, magnetar=Magnetar(L_0=1e46, t_0=100, q=2))

.. note::
    The magnetar spin-down injection is implemented in the default form L_0*(1+t/t_0)^(-q) for theta < theta_c. You can pass the `magnetar` argument to the power-law and Gaussian jet as well.


User-Defined Jet
^^^^^^^^^^^^^^^^

You may also define your own jet structure by providing the energy and lorentz factor profile.
Those two profiles are required to complete a jet structure. You may also provide the magnetization profile, enregy injection profile, and mass injection profile.
Those profiles are optional and will be set to zero function if not provided.

.. code-block:: python

    from VegasAfterglow import Ejecta

    # Define a custom energy profile function, required to complete the jet structure
    def E_iso_profile(phi, theta):
        return 1e53  # E_iso = 1e53 erg isotropic fireball
        #return whatever energy profile you want as a function of phi and theta in unit of erg [not erg per solid angle]

    # Define a custom lorentz factor profile function, required to complete the jet structure
    def Gamma0_profile(phi, theta):
        return 300 # Gamma0 = 300
        #return whatever lorentz factor profile you want as a function of phi and theta
    
    # Define a custom magnetization profile function, optional
    def sigma0_profile(phi, theta):
        return 0.1 # sigma = 0.1
        #return whatever magnetization profile you want as a function of phi and theta

    # Define a custom energy injection profile function, optional
    def E_dot_profile(phi, theta, t):
        return 1e46 * (1 + t / 100)**(-2) # L = 1e46 erg/s, t0 = 100 seconds
        #return whatever energy injection  profile you want as a function of phi, theta, and time in unit of erg/s [not erg/s per solid angle]

    # Define a custom mass injection profile function, optional
    def M_dot_profile(phi, theta, t):
        #return whatever mass injection profile you want as a function of phi, theta, and time in unit of g/s [not g/s per solid angle]

    # Create a user-defined jet
    jet = Ejecta(E_iso=E_iso_profile, Gamma0=Gamma0_profile, sigma0=sigma0_profile, E_dot=E_dot_profile, M_dot=M_dot_profile)

    #..other settings

    #if your jet is not axisymmetric, set axisymmetric to False
    model = Model(jet=jet, ..., axisymmetric=False, resolutions=(0.3, 5, 10))

    # the user-defined jet structure could be spiky, if the default resolution may not resolve the jet structure, if that is the case, you can try a finer resolution (phi_ppd, theta_ppd, t_ppd)
    # where phi_ppd is the number of points per degree in the phi direction, theta_ppd is the number of points per degree in the theta direction, and t_ppd is the number of points per decade in the time direction    .
    
.. note::
    Setting user-defined structured jet in the Python level is OK for light curve and spectrum calculation. However, it is not recommended for MCMC parameter fitting if you do care about the performance.
    The reason is that setting user-defined profiles in the Python level leads to a large overhead due to the Python-C++ inter-process communication.
    Users are recommended to set up the user-defined jet structure in the C++ level for MCMC parameter fitting for better performance.
          

Radiation Processes
-------------------

Reverse Shock Emission
^^^^^^^^^^^^^^^^^^^^^^    

.. code-block:: python

    from VegasAfterglow import Radiation

    #set the jet duration to be 100 seconds, the default is 1 second. The jet duration affects the reverse shock thickness (thin shell or thick shell).
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300, duration = 100)

    # Create a radiation model with self-Compton radiation
    fwd_rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    rvs_rad = Radiation(eps_e=1e-2, eps_B=1e-4, p=2.4)

    #..other settings
    model = Model(forward_rad=fwd_rad, reverse_rad=rvs_rad, resolutions=(0.5, 10, 10),...)

    times = np.logspace(2, 8, 200)  

    bands = np.array([1e9, 1e14, 1e17])  

    results = model.specific_flux(times, bands)
    
    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band 
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')
        plt.loglog(times, results['syn_rvs'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')#reverse shock synchrotron

.. note::
    You may increase the resolution of the grid to improve the accuracy of the reverse shock synchrotron radiation.


Inverse Compton Cooling
^^^^^^^^^^^^^^^^^^^^^^^    

.. code-block:: python

    from VegasAfterglow import Radiation

    # Create a radiation model with inverse Compton cooling (without Klein-Nishina correction) on synchrotron radiation
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, IC_cooling=True, KN=False)

    #..other settings
    model = Model(forward_rad=rad, ...)

Self-Synchrotron Compton Radiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Radiation

    # Create a radiation model with self-Compton radiation
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, SSC=True, KN=True, IC_cooling=True)

    #..other settings
    model = Model(forward_rad=rad, ...)

    times = np.logspace(2, 8, 200)  

    bands = np.array([1e9, 1e14, 1e17])  

    results = model.specific_flux(times, bands)

    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band 
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results['syn'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')#synchrotron
        plt.loglog(times, results['IC'][i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')#SSC

.. note::
    (IC_cooling = False, KN = False, SSC = True): The IC radiation is calculated based on synchrotron spectrum without IC cooling.

    (IC_cooling = True, KN = False, SSC = True): The IC radiation is calculated based on synchrotron spectrum with IC cooling, but without Klein-Nishina correction.

    (IC_cooling = True, KN = True, SSC = True): The IC radiation is calculated based on synchrotron spectrum with both IC cooling and Klein-Nishina correction.
  
Advanced Features
-----------------


MCMC Parameter Fitting
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import ObsData, Fitter, ParamDef, Scale

    # Create observation data object
    data = ObsData()

    # Add some observational data (light curves)
    t_data = np.array([1e3, 2e3, 5e3, 1e4, 2e4])  # Time in seconds
    flux_data = np.array([1e-26, 8e-27, 5e-27, 3e-27, 2e-27])  # Specific flux
    flux_err = np.array([1e-28, 8e-28, 5e-28, 3e-28, 2e-28])  # Flux error
    
    # Add a light curve at optical frequency (5e14 Hz)
    data.add_light_curve(nu=5e14, t=t_data, flux=flux_data, flux_err=flux_err)
    
    # Define parameters with priors
    params = [
        ParamDef("E_iso",      1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
        ParamDef("Gamma0",        5,  1000,  Scale.LOG),       # Lorentz factor at the core
        ParamDef("theta_c",     0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
        ParamDef("theta_v",     0.0,   0.0,  Scale.FIXED),     # Viewing angle [rad]
        ParamDef("p",             2,     3,  Scale.LINEAR),    # Shocked electron power law index
        ParamDef("eps_e",      1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
        ParamDef("eps_B",      1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
        ParamDef("A_star",     1e-3,     1,  Scale.LOG),       # Wind parameter
        ParamDef("xi_e",       1e-3,     1,  Scale.LOG),       # Electron acceleration fraction
    ]
    
    # Create the fitter with default model setup
    fitter = Fitter(data=data, params=params)
    
    # Run MCMC
    samples, log_probs = fitter.run_mcmc(
        n_walkers=32,  # Number of walkers
        n_steps=1000,  # Number of steps per walker
        n_burn=200,    # Number of burn-in steps to discard
        progress=True  # Show progress bar
    )
    
    # Plot the posterior distributions
    fitter.plot_corner()

Parameter Study
^^^^^^^^^^^^^^^

.. code-block:: python

    # Study the effect of electron energy index p
    p_values = np.linspace(2.0, 3.0, 5)
    
    plt.figure(figsize=(10, 6))
    
    # Fix a frequency to study (optical)
    nu_index = 1  # Optical band
    
    for p in p_values:
        # Update the radiation model
        model.radiation.p = p
        
        # Calculate new light curve
        results_p = model.calculate_light_curves(times, frequencies)
        
        # Plot
        plt.loglog(times, results_p[:, nu_index], label=f'p = {p:.1f}')
    
    plt.xlabel('Time (s)')
    plt.ylabel('Flux Density (erg/cm²/s/Hz)')
    plt.legend()
    plt.title('Effect of Electron Energy Index (p) on Optical Light Curves')
    plt.grid(True, which='both', linestyle='--', alpha=0.5)
    plt.show()
