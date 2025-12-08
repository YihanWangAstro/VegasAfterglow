Python API Reference
====================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

The Python API provides a user-friendly interface to VegasAfterglow's C++ core, enabling easy model setup, simulation execution, and result analysis. All major C++ components are exposed through Python bindings with an intuitive interface.

Key Components
--------------

The Python API is organized into several core components:

* **Data Representation**: Classes for handling observational data and model parameters
* **MCMC Framework**: Tools for parameter fitting and posterior exploration
* **Model Configuration**: Classes for setting up the physical model
* **Results Processing**: Tools for handling simulation outputs

Core Classes
------------

.. _api-obsdata:

ObsData
^^^^^^^

.. autoclass:: VegasAfterglow.ObsData
   :members:
   :undoc-members:
   :show-inheritance:

The `ObsData` class is used to store and manage observational data, including light curves and spectra. It provides methods to add observational data from various sources.

Example:

.. code-block:: python

    from VegasAfterglow import ObsData

    # Create an instance to store observational data
    data = ObsData()

    # Add light curve data
    data.add_flux_density(nu=4.84e14, t=time_data, f_nu=flux_data, err=flux_error)  # All quantities in CGS units

    # Add spectrum data
    data.add_spectrum(t=3000, nu=nu_data, f_nu=spectrum_data, err=spectrum_error)  # All quantities in CGS units

.. _api-setups:

Setups
^^^^^^

.. autoclass:: VegasAfterglow.Setups
   :members:
   :undoc-members:
   :show-inheritance:

The `Setups` class defines global properties and environment settings for the model. These settings remain fixed during the MCMC process.

Example:

.. code-block:: python

    from VegasAfterglow import Setups

    # Create configuration
    cfg = Setups()

    # Source properties
    cfg.lumi_dist = 3.364e28    # Luminosity distance [cm]
    cfg.z = 1.58               # Redshift

    # Physical model configuration
    cfg.medium = "wind"        # Ambient medium: "wind", "ISM", or "user"
    cfg.jet = "powerlaw"       # Jet structure: "powerlaw", "gaussian", "tophat", or "user"

.. _api-modelparams:

ModelParams
^^^^^^^^^^^

.. autoclass:: VegasAfterglow.ModelParams
   :members:
   :undoc-members:
   :show-inheritance:

The `ModelParams` class stores the physical parameters that define the GRB afterglow model. These parameters are varied during the MCMC fitting process.

.. _api-paramdef:

ParamDef and Scale
^^^^^^^^^^^^^^^^^^

.. _ParamDef:

.. autoclass:: VegasAfterglow.ParamDef
   :members:
   :undoc-members:
   :show-inheritance:

.. _api-scale:
.. _api-scale-log:
.. _api-scale-linear:
.. _api-scale-fixed:

.. autoclass:: VegasAfterglow.Scale
   :members:
   :undoc-members:
   :show-inheritance:

These classes are used to define parameters for MCMC exploration, including their name, initial value, prior range, and sampling scale:

Example:

.. code-block:: python

    from VegasAfterglow import ParamDef, Scale

    mc_params = [
        ParamDef("E_iso",   1e50,  1e54,  Scale.LOG),       # Isotropic energy [erg]
        ParamDef("Gamma0",     5,  1000,  Scale.LOG),       # Lorentz factor at the core
        ParamDef("theta_c",  0.0,   0.5,  Scale.LINEAR),    # Core half-opening angle [rad]
        ParamDef("theta_v",  0.0,   0.0,  Scale.FIXED),     # Viewing angle [rad]
        ParamDef("p",          2,     3,  Scale.LINEAR),    # Shocked electron power law index
        ParamDef("eps_e",   1e-2,   0.5,  Scale.LOG),       # Electron energy fraction
        ParamDef("eps_B",   1e-4,   0.5,  Scale.LOG),       # Magnetic field energy fraction
        ParamDef("A_star",  1e-3,     1,  Scale.LOG),       # Wind parameter
    ]

.. _A_star:

For wind medium models, you would use the A_star parameter as shown above.

.. _n_ism:

For ISM medium models, you would use the density parameter instead:

.. code-block:: python

    ParamDef("n_ism",     0.1,  1e-3,    10,  Scale.LOG),       # ISM density [cm^-3]

For a comprehensive list of all available parameters, their physical meanings, typical ranges, and usage guidelines, see the :doc:`parameter_reference` page.

.. _api-fitter:

Fitter
^^^^^^

.. autoclass:: VegasAfterglow.Fitter
   :members:
   :undoc-members:
   :show-inheritance:

The `Fitter` class provides a high-level interface for MCMC fitting of GRB afterglow models to observational data.

Example:

.. code-block:: python

    from VegasAfterglow import Fitter

    # Create the fitter object
    fitter = Fitter(data, cfg)

    # Option 1: MCMC with emcee (faster, recommended for quick fitting)
    result = fitter.fit(
        mc_params,
        resolution=(0.3, 1, 10),       # Grid resolution (phi, theta, t)
        sampler="emcee",               # MCMC sampler
        nsteps=10000,                  # Number of steps per walker
        nburn=1000,                    # Burn-in steps to discard
        npool=8,                       # Number of parallel processes
        top_k=10,                      # Number of best-fit parameters to return
    )

    # Option 2: Nested sampling with dynesty (slower but computes Bayesian evidence)
    result = fitter.fit(
        mc_params,
        resolution=(0.3, 1, 10),       # Grid resolution (phi, theta, t)
        sampler="dynesty",             # Nested sampling algorithm
        nlive=500,                     # Number of live points
        dlogz=0.1,                     # Stopping criterion (evidence tolerance)
        sample="rwalk",                # Sampling method
        npool=8,                       # Number of parallel processes
    )

    # Generate light curves with best-fit parameters
    lc_best = fitter.flux_density_grid(result.top_k_params[0], t_out, bands)

    # Generate spectra with best-fit parameters
    spec_best = fitter.flux_density_grid(result.top_k_params[0], times, nu_out)

.. _api-fitresult:

FitResult
^^^^^^^^^

.. autoclass:: VegasAfterglow.FitResult
   :members:
   :undoc-members:
   :show-inheritance:

The `FitResult` class stores the results of an MCMC fit, including the posterior samples, log probabilities, top-k best-fit parameters, and the full bilby Result object for diagnostics.

.. _api-vegasmc:

VegasMC
^^^^^^^

.. autoclass:: VegasAfterglow.VegasMC
   :members:
   :undoc-members:
   :show-inheritance:

The `VegasMC` class is the core calculator for MCMC sampling, providing efficient computation of model likelihood based on the specified parameters.

Documenting Python Code
-----------------------

When contributing to the Python codebase, please follow these documentation guidelines:

Class and Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use NumPy-style docstrings for all classes and functions:

.. code-block:: python

    def function(param1, param2):
        """
        Brief description of the function.

        Detailed description of the function's behavior, expected inputs,
        outputs, and any other relevant information.

        Parameters
        ----------
        param1 : type
            Description of param1
        param2 : type
            Description of param2

        Returns
        -------
        type
            Description of the return value

        Examples
        --------
        >>> function(1, 2)
        3
        """

Example Class
^^^^^^^^^^^^^

Here's an example of a well-documented class:

.. code-block:: python

    class ParamDef:
        """
        Single-parameter definition for MCMC.

        This class defines a parameter to be used in MCMC fitting, including
        its name, initial value, prior range, and sampling scale.

        Parameters
        ----------
        name : str
            The parameter name
        init : float
            Initial value for the parameter
        lower : float, optional
            Lower bound for the parameter (not needed for FIXED scale)
        upper : float, optional
            Upper bound for the parameter (not needed for FIXED scale)
        scale : Scale, optional
            Sampling scale (LINEAR, LOG, or FIXED), default is LINEAR

        Notes
        -----
        When scale=LOG, we sample log10(x), then transform via 10**v.
        When scale=FIXED, the parameter never appears in the sampler.
        """

For more details on NumPy-style docstrings, see the :doc:`contributing` page.
