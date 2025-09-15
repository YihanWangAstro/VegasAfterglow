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
    model = Model(jet=jet, medium=medium, observer=obs, fwd_rad=rad, resolutions=(0.3,1,10))

    # Define time range for light curve calculation
    times = np.logspace(2, 8, 200)

    # Define observing frequencies (radio, optical, X-ray bands in Hz)
    bands = np.array([1e9, 1e14, 1e17])

    # Calculate the afterglow emission at each time and frequency
    # NOTE that the times array needs to be in ascending order
    results = model.specific_flux(times, bands)

    # Visualize the multi-wavelength light curves
    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results.total[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz')

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
        plt.loglog(frequencies, results.total[:,i], color=colors[i], label=fr'${base:.1f} \times 10^{{{exp}}}$ s')

    # Add vertical lines marking the bands from the light curve plot
    for i, band in enumerate(bands):
        exp = int(np.floor(np.log10(band)))
        base = band / 10**exp
        plt.axvline(band,ls='--',color='C'+str(i))

    plt.xlabel('frequency (Hz)')
    plt.ylabel('flux density (erg/cm²/s/Hz)')
    plt.legend(ncol=2)
    plt.title('Synchrotron Spectra')

Calculate flux on time-frequency pairs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Suppose you want to calculate the flux at specific time-frequency pairs (t_i, nu_i) instead of a grid (t_i, nu_j), you can use the following method:

.. code-block:: python

    # Define time range for light curve calculation
    times = np.logspace(2, 8, 200)

    # Define observing frequencies (must be the same length as times)
    bands = np.logspace(9,17, 200)

    results = model.specific_flux_series(times, bands) #times array could be random order

    # the returned results is a FluxDict object with arrays of the same shape as the input times and bands.

Calculate flux with exposure time averaging
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For observations with finite exposure times, you can calculate time-averaged flux by sampling multiple points within each exposure:

.. code-block:: python

    # Define observation times (start of exposure)
    times = np.logspace(2, 8, 50)

    # Define observing frequencies (must be the same length as times)
    bands = np.logspace(9, 17, 50)

    # Define exposure times for each observation (in seconds)
    expo_time = np.ones_like(times) * 100  # 100-second exposures

    # Calculate time-averaged flux with 20 sample points per exposure
    results = model.specific_flux_series_with_expo(times, bands, expo_time, num_points=20)

    # The returned results is a FluxDict object with arrays of the same shape as input times and bands
    # Each flux value represents the average over the corresponding exposure time

.. note::
    The function samples `num_points` evenly spaced within each exposure time and averages the computed flux. Higher `num_points` gives more accurate time averaging but increases computation time. The minimum value is 2.


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

Stratified Medium
^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Wind

    # Create a stratified stellar wind medium;
    # smooth transited stratified medium. Inner region, n(r) = n_0, middle region n(r) \propto 1/r^2, outer region n(r)=n_ism
    # A = 0 (default): fallback to n = n_ism
    # n_0 = inf (default): wind bubble, from wind profile to ism profile
    # A = 0 & n_0 = inf: pure wind;
    wind = Wind(A_star=0.1, n_ism = 1, n_0 = 1e-3)

    #..other settings
    model = Model(medium=wind, ...)


User-Defined Medium
^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Medium

    mp = 1.67e-24 # proton mass in gram

    # Define a custom density profile function
    def density(phi, theta, r):# r in cm, phi and theta in radians
        return mp # n_ism =  1 cm^-3
        #return whatever density profile (cm^-3) you want as a function of phi, theta, and r

    # Create a user-defined medium
    medium = Medium(rho=density)

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
        k_e=2.0,              # Power-law index for energy angular dependence
        k_g=2.0               # Power-law index for Lorentz factor angular dependence
    )

    #..other settings
    model = Model(jet=jet, ...)

Two-Component Jet
^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import TwoComponentJet

    # Create a two-component jet
    jet = TwoComponentJet(
        theta_c=0.05,        # Narrow component angular size (radians)
        E_iso_c=1e53,        # Isotropic-equivalent energy of the narrow component (ergs)
        Gamma0_c=300,        # Initial Lorentz factor of the narrow component
        theta_w=0.1,         # Wide component angular size (radians)
        E_iso_w=1e52,        # Isotropic-equivalent energy of the wide component (ergs)
        Gamma0_w=100         # Initial Lorentz factor of the wide component
    )

    #..other settings
    model = Model(jet=jet, ...)

Step Power-Law Jet
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import StepPowerLawJet

    # Create a step power-law structured jet (uniform core with sharp transition)
    jet = StepPowerLawJet(
        theta_c=0.05,        # Core angular size (radians)
        E_iso_c=1e53,        # Isotropic-equivalent energy of the core component (ergs)
        Gamma0_c=300,        # Initial Lorentz factor of the core component
        E_iso_w=1e52,        # Isotropic-equivalent energy of the wide component (ergs)
        Gamma0_w=100,        # Initial Lorentz factor of the wide component
        k_e=2.0,             # Power-law index for energy angular dependence
        k_g=2.0              # Power-law index for Lorentz factor angular dependence
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
    jet = TophatJet(theta_c=0.05, E_iso=1e53, Gamma0=300, magnetar=Magnetar(L0=1e46, t0=100, q=2))

.. note::
    The magnetar spin-down injection is implemented in the default form L0*(1+t/t0)^(-q) for theta < theta_c. You can pass the `magnetar` argument to the power-law and Gaussian jet as well.


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
    model = Model(jet=jet, ..., axisymmetric=False, resolutions=(0.3, 1, 10))

    # the user-defined jet structure could be spiky, the default resolution may not resolve the jet structure. if that is the case, you can try a finer resolution (phi_ppd, theta_ppd, t_ppd)
    # where phi_ppd is the number of points per degree in the phi direction, theta_ppd is the number of points per degree in the theta direction, and t_ppd is the number of points per decade in the time direction    .

.. note::
    Setting user-defined structured jet in the Python level is OK for light curve and spectrum calculation. However, it is not recommended for MCMC parameter fitting if you do care about the performance.
    The reason is that setting user-defined profiles in the Python level leads to a large overhead due to the Python-C++ inter-process communication.
    Users are recommended to set up the user-defined jet structure in the C++ level for MCMC parameter fitting for better performance, if you want the best performance.


Radiation Processes
-------------------

Reverse Shock Emission
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Radiation

    #set the jet duration to be 100 seconds, the default is 1 second. The jet duration affects the reverse shock thickness (thin shell or thick shell).
    jet = TophatJet(theta_c=0.1, E_iso=1e52, Gamma0=300, duration = 100)

    # Create a radiation model with both forward and reverse shock synchrotron radiation
    fwd_rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3)
    rvs_rad = Radiation(eps_e=1e-2, eps_B=1e-4, p=2.4)

    #..other settings
    model = Model(fwd_rad=fwd_rad, rvs_rad=rvs_rad, resolutions=(0.5, 1, 10),...)

    times = np.logspace(2, 8, 200)

    bands = np.array([1e9, 1e14, 1e17])

    results = model.specific_flux(times, bands)

    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results.fwd.sync[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz (fwd)')
        plt.loglog(times, results.rvs.sync[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz (rvs)')#reverse shock synchrotron

.. note::
    You may increase the resolution of the grid to improve the accuracy of the reverse shock synchrotron radiation if you see spiky features.


Inverse Compton Cooling
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Radiation

    # Create a radiation model with inverse Compton cooling (without Klein-Nishina correction) on synchrotron radiation
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, ssc_cooling=True, kn=False)

    #..other settings
    model = Model(fwd_rad=rad, ...)

Self-Synchrotron Compton Radiation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from VegasAfterglow import Radiation

    # Create a radiation model with self-Compton radiation
    rad = Radiation(eps_e=1e-1, eps_B=1e-3, p=2.3, ssc=True, kn=True, ssc_cooling=True)

    #..other settings
    model = Model(fwd_rad=rad, ...)

    times = np.logspace(2, 8, 200)

    bands = np.array([1e9, 1e14, 1e17])

    results = model.specific_flux(times, bands)

    plt.figure(figsize=(4.8, 3.6),dpi=200)

    # Plot each frequency band
    for i, nu in enumerate(bands):
        exp = int(np.floor(np.log10(nu)))
        base = nu / 10**exp
        plt.loglog(times, results.fwd.sync[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz (sync)')#synchrotron
        plt.loglog(times, results.fwd.ssc[i,:], label=fr'${base:.1f} \times 10^{{{exp}}}$ Hz (SSC)')#SSC

.. note::
    (ssc_cooling = False, kn = False, ssc = True): The IC radiation is calculated based on synchrotron spectrum without IC cooling.

    (ssc_cooling = True, kn = False, ssc = True): The IC radiation is calculated based on synchrotron spectrum with IC cooling, but without Klein-Nishina correction.

    (ssc_cooling = True, kn = True, ssc = True): The IC radiation is calculated based on synchrotron spectrum with both IC cooling and Klein-Nishina correction.
