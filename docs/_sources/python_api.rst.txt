Python API Reference
===================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
-------

The Python API provides a user-friendly interface to VegasAfterglow's C++ core, enabling easy model setup, simulation execution, and result analysis. All major C++ components are exposed through Python bindings with an intuitive interface.

Key Components
-------------

The Python API is organized into several core components:

* **Model Creation**: Classes for defining GRB afterglow models
* **Jet Models**: Classes for different jet structure profiles (TophatJet, GaussianJet, PowerLawJet)
* **Ambient Media**: Classes for defining the circumburst environment (ISM, Wind)
* **Radiation Processes**: Components for calculating emission mechanisms (Synchrotron, SynchrotronSelfCompton)
* **Shock Dynamics**: Forward and reverse shock physics
* **Data Fitting**: Tools for comparing models with observational data using MCMC methods
* **Results Processing**: Classes for handling and visualizing simulation results

Core Classes
-----

.. autoclass:: VegasAfterglow.Model
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Observer
   :members:
   :undoc-members:
   :show-inheritance:

Jet Models
---------

.. autoclass:: VegasAfterglow.TophatJet
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.GaussianJet
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.PowerLawJet
   :members:
   :undoc-members:
   :show-inheritance:

Ambient Media
------------

.. autoclass:: VegasAfterglow.ISM
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Wind
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.UserDefinedMedium
   :members:
   :undoc-members:
   :show-inheritance:

Radiation Processes
-----------------

.. autoclass:: VegasAfterglow.Synchrotron
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.SynchrotronSelfCompton
   :members:
   :undoc-members:
   :show-inheritance:

Shock Dynamics
------------

.. autoclass:: VegasAfterglow.ForwardShock
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.ReverseShock
   :members:
   :undoc-members:
   :show-inheritance:

Data Fitting
-----------

.. autoclass:: VegasAfterglow.ObsData
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Fitter
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.ParamDef
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Scale
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Setups
   :members:
   :undoc-members:
   :show-inheritance:

Results
------

.. autoclass:: VegasAfterglow.SimulationResults
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.LightCurve
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Spectrum
   :members:
   :undoc-members:
   :show-inheritance:

Utilities
-------

.. autoclass:: VegasAfterglow.Constants
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: VegasAfterglow.Units
   :members:
   :undoc-members:
   :show-inheritance:

Documenting Python Code
---------------------

When contributing to the Python codebase, please follow these documentation guidelines:

Class and Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^

Here's an example of a well-documented class:

.. code-block:: python

    class GaussianJet:
        """
        A class representing a structured jet model with Gaussian profile.

        This class implements a jet with angular structure, where the energy
        and Lorentz factor follow a Gaussian profile with angle from the jet axis.

        Parameters
        ----------
        core_energy : float
            The isotropic-equivalent energy in the jet core (ergs)
        core_lorentz : float
            The initial bulk Lorentz factor in the jet core
        theta_core : float
            The angular size of the jet core (radians)
        viewer_angle : float, optional
            The viewing angle relative to the jet axis (radians)

        Attributes
        ----------
        theta_core : float
            Core angle of the jet (radians)
        """

        def get_energy_at_angle(self, theta):
            """
            Get the energy at a specific angle.

            Parameters
            ----------
            theta : float
                Angle from the jet axis (radians)

            Returns
            -------
            float
                Energy at the specified angle (ergs)
            """

For more details on NumPy-style docstrings, see the :doc:`contributing` page. 