Test Class Documentation
========================

This is a test page to verify that Breathe can properly process the Doxygen XML files.

GaussianJet Class
----------------

.. doxygenclass:: GaussianJet
   :project: VegasAfterglow
   :members:
   :private-members:
   :protected-members:
   :undoc-members:
   :no-link:
   :outline:
   :allow-dot-graphs:

Details for GaussianJet Functions
--------------------------------

Constructor
^^^^^^^^^^

.. doxygenfunction:: GaussianJet::GaussianJet
   :project: VegasAfterglow
   :no-link:

GaussianJet is a class that implements a Gaussian-structured relativistic jet model, where the energy and Lorentz factor distribution follow a Gaussian profile with angle. The constructor initializes a new GaussianJet with the following parameters:

* **theta_c**: Core angle (radians) - defines the angular scale of the Gaussian distribution 
* **E_iso**: Isotropic equivalent energy (ergs) - the total energy if viewed along the jet axis
* **Gamma0**: Initial Lorentz factor along the jet axis

This model is commonly used for GRB afterglows where observations suggest structured jets rather than uniform top-hat jets.

Energy Distribution
^^^^^^^^^^^^^^^^^

.. doxygenfunction:: GaussianJet::eps_k
   :project: VegasAfterglow
   :no-link:

The eps_k function calculates the isotropic-equivalent kinetic energy per solid angle as a function of the polar angle theta.

For a Gaussian jet, the energy distribution follows:

.. math::
   \epsilon_k(\theta) = \epsilon_{k,0} \exp\left(-\frac{\theta^2}{2\theta_c^2}\right)

where:
* :math:`\epsilon_{k,0}` is the on-axis energy per solid angle
* :math:`\theta_c` is the core angle of the jet
* :math:`\theta` is the polar angle from the jet axis

This function helps model how the energy is distributed across different viewing angles of the jet.

Lorentz Factor Profile
^^^^^^^^^^^^^^^^^^^^

.. doxygenfunction:: GaussianJet::Gamma0
   :project: VegasAfterglow
   :no-link:

The Gamma0 function returns the initial Lorentz factor as a function of polar angle theta.

For a Gaussian jet, the Lorentz factor follows a similar profile to the energy distribution:

.. math::
   \Gamma_0(\theta) = 1 + (\Gamma_{0,c} - 1) \exp\left(-\frac{\theta^2}{2\theta_c^2}\right)

where:
* :math:`\Gamma_{0,c}` is the on-axis maximum Lorentz factor
* :math:`\theta_c` is the core angle of the jet
* :math:`\theta` is the polar angle from the jet axis

The Lorentz factor determines the initial velocity of the ejecta at different angles, which affects the observed timescales and brightness of the afterglow.

Class Members
^^^^^^^^^^

.. cpp:member:: Real GaussianJet::T0

The deceleration time of the jet (in seconds). This is calculated based on the initial kinetic energy and circumburst medium density.

.. cpp:member:: bool GaussianJet::spreading

Boolean flag indicating whether the jet undergoes lateral spreading. When true, the jet expands sideways as it decelerates, which affects the observed light curve, particularly at late times when the jet has decelerated significantly.

Private Members
^^^^^^^^^^^^^

.. cpp:member:: Real const GaussianJet::norm_

Normalization factor for the energy distribution, calculated based on the total isotropic energy and core angle.

.. cpp:member:: Real const GaussianJet::eps_k_

The on-axis energy per solid angle value, used in the Gaussian energy distribution formula.

.. cpp:member:: Real const GaussianJet::Gamma0_

The on-axis initial Lorentz factor, representing the maximum Lorentz factor at the jet center. 