Parameter Reference
===================

This page provides a comprehensive reference for all parameters used in VegasAfterglow, including their physical meanings, typical ranges, and units. All parameters listed here are available in the code and can be set via Python interfaces.

Physical Parameters
-------------------

Observer Parameters
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``lumi_dist``
     - :math:`d_L`
     - cm
     - :math:`10^{26} - 10^{29}`
     - Luminosity distance to the source
   * - ``z``
     - :math:`z`
     - dimensionless
     - :math:`0.01 - 10`
     - Cosmological redshift
   * - ``theta_v``
     - :math:`\theta_v`
     - radians
     - :math:`0 - \pi/2`
     - Viewing angle (angle between jet axis and line of sight)

Jet Structure Parameters
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``E_iso``
     - :math:`E_{\rm iso}`
     - erg
     - :math:`10^{50} - 10^{54}`
     - Isotropic-equivalent kinetic energy of the jet
   * - ``Gamma0``
     - :math:`\Gamma_0`
     - dimensionless
     - :math:`10 - 1000`
     - Initial bulk Lorentz factor of the jet
   * - ``theta_c``
     - :math:`\theta_c`
     - radians
     - :math:`0.01 - 0.5`
     - Half-opening angle of the jet core
   * - ``duration``
     - :math:`T_{\rm dur}`
     - seconds
     - :math:`0.1 - 1000`
     - Duration of energy injection (affects reverse shock)
   * - ``k_e``
     - :math:`k_e`
     - dimensionless
     - :math:`1 - 10`
     - Energy power-law index for structured jets (PowerLawJet only)
   * - ``k_g``
     - :math:`k_g`
     - dimensionless
     - :math:`1 - 10`
     - Lorentz factor power-law index for structured jets (PowerLawJet only)

Two-Component Jet Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``theta_w``
     - :math:`\theta_w`
     - radians
     - :math:`0.1 - 0.5`
     - Half-opening angle of wide component
   * - ``E_iso_w``
     - :math:`E_{\rm iso,w}`
     - erg
     - :math:`10^{50} - 10^{53}`
     - Isotropic energy of wide component
   * - ``Gamma0_w``
     - :math:`\Gamma_{0,w}`
     - dimensionless
     - :math:`10 - 300`
     - Initial Lorentz factor of wide component

Ambient Medium Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``n_ism``
     - :math:`n_{\rm ISM}`
     - cm⁻³
     - :math:`10^{-4} - 10^{3}`
     - Number density of uniform ISM
   * - ``n0``
     - :math:`n0`
     - cm⁻³
     - :math:`10^{-4} - 10^{6}`
     - Inner region number density for wind medium
   * - ``A_star``
     - :math:`A_*`
     - dimensionless
     - :math:`10^{-3} - 10`
     - Wind parameter: :math:`\rho = A_* \times 5 \times 10^{11} r^{-2}` g/cm³

Forward Shock Radiation Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``eps_e``
     - :math:`\epsilon_e`
     - dimensionless
     - :math:`10^{-3} - 0.5`
     - Fraction of shock energy in relativistic electrons
   * - ``eps_B``
     - :math:`\epsilon_B`
     - dimensionless
     - :math:`10^{-6} - 0.5`
     - Fraction of shock energy in magnetic field
   * - ``p``
     - :math:`p`
     - dimensionless
     - :math:`2.01 - 3.5`
     - Power-law index of electron energy distribution
   * - ``xi_e``
     - :math:`\xi_e`
     - dimensionless
     - :math:`10^{-3} - 1`
     - Fraction of electrons that are accelerated

Reverse Shock Radiation Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``eps_e_r``
     - :math:`\epsilon_{e,r}`
     - dimensionless
     - :math:`10^{-3} - 0.5`
     - Reverse shock fraction of energy in electrons
   * - ``eps_B_r``
     - :math:`\epsilon_{B,r}`
     - dimensionless
     - :math:`10^{-6} - 0.5`
     - Reverse shock fraction of energy in magnetic field
   * - ``p_r``
     - :math:`p_r`
     - dimensionless
     - :math:`2.01 - 3.5`
     - Reverse shock electron energy distribution index
   * - ``xi_e_r``
     - :math:`\xi_{e,r}`
     - dimensionless
     - :math:`10^{-3} - 1`
     - Reverse shock electron acceleration fraction

Energy Injection Parameters (Magnetar)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 20 40

   * - Parameter
     - Symbol
     - Units
     - Typical Range
     - Description
   * - ``L0``
     - :math:`L_0`
     - erg/s
     - :math:`10^{44} - 10^{48}`
     - Magnetar luminosity at time t₀
   * - ``t0``
     - :math:`t_0`
     - seconds
     - :math:`10 - 10^4`
     - Characteristic magnetar spin-down timescale
   * - ``q``
     - :math:`q`
     - dimensionless
     - :math:`1 - 6`
     - Power-law index of spin-down: :math:`L(t) = L_0(1+t/t_0)^{-q}`

Model Configuration
-------------------

Jet Types
^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Jet Type
     - Description
   * - ``tophat``
     - Uniform energy and Lorentz factor within opening angle
   * - ``gaussian``
     - Gaussian angular profile for energy and Lorentz factor
   * - ``powerlaw``
     - Power-law angular dependence with indices k_e and k_g
   * - ``two_component``
     - Two-component jet with narrow core and wide wing components
   * - ``step_powerlaw``
     - Uniform core with sharp transition to power-law wing
   * - ``ejecta``
     - Generic ejecta with arbitrary angular profiles

Medium Types
^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Medium Type
     - Description
   * - ``ism``
     - Uniform interstellar medium with constant density n_ism
   * - ``wind``
     - Stellar wind medium with :math:`\rho \propto r^{-2}` profile

Physics Switches
^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Parameter
     - Default
     - Description
   * - ``rvs_shock``
     - false
     - Include reverse shock emission
   * - ``fwd_ssc``
     - false
     - Include forward shock synchrotron self-Compton
   * - ``rvs_ssc``
     - false
     - Include reverse shock synchrotron self-Compton
   * - ``ssc_cooling``
     - false
     - Include inverse Compton cooling
   * - ``kn``
     - false
     - Use Klein-Nishina cross-section for IC scattering
   * - ``magnetar``
     - false
     - Include magnetar energy injection

Computational Parameters
------------------------

Grid Resolution
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 10 60

   * - Parameter
     - Default
     - Units
     - Description
   * - ``phi_resol``
     - 0.3
     - points/degree
     - Angular resolution in azimuthal direction
   * - ``theta_resol``
     - 1.0
     - points/degree
     - Angular resolution in polar direction
   * - ``t_resol``
     - 10.0
     - points/decade
     - Temporal resolution (logarithmic spacing)

Numerical Parameters
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - Parameter
     - Default
     - Description
   * - ``rtol``
     - 1e-6
     - Relative tolerance for numerical integration

MCMC Parameters
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 20 60

   * - Parameter
     - Typical Value
     - Description
   * - ``total_steps``
     - 1000-50000
     - Total number of MCMC steps per walker
   * - ``burn_frac``
     - 0.2-0.5
     - Fraction of steps to discard as burn-in
   * - ``thin``
     - 1-10
     - Thinning factor (keep every nth sample)
   * - ``n_walkers``
     - 2×n_params to 10×n_params
     - Number of ensemble walkers

Parameter Scaling Types
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 20 80

   * - Scale Type
     - Description and Usage
   * - ``Scale.LOG``
     - Sample in log₁₀ space. Use for parameters spanning multiple orders of magnitude (energies, densities, microphysics parameters)
   * - ``Scale.LINEAR``
     - Sample in linear space. Use for parameters with limited ranges (angles, power-law indices)
   * - ``Scale.FIXED``
     - Keep parameter fixed at initial value. Use when you don't want to vary a parameter

Parameter Relationships and Constraints
---------------------------------------

Physical Constraints
^^^^^^^^^^^^^^^^^^^^

**Energy Conservation:**

- :math:`E_{\rm iso}` should be consistent with the kinetic energy available from the central engine

**Causality:**

- Light travel time sets minimum variability timescale: :math:`\delta t \geq R/c\Gamma^2`
- Jet opening angle and Lorentz factor: :math:`\theta_c \gtrsim 1/\Gamma_0` for causal contact

**Microphysics:**

- Energy fractions: :math:`\epsilon_e + \epsilon_B \leq 1` (though often :math:`\ll 1`)
- Electron power-law index: :math:`p > 2` for finite energy in fast-cooling regime

Unit System and Physical Constants
----------------------------------

VegasAfterglow uses a normalized unit system defined in ``macros.h``:

**Base Units:**
- Length: :math:`l_0 = 1.5 \times 10^{13}` cm
- Time: :math:`t_0 = l_0/c = 500` s
- Mass: :math:`m_0 = 2 \times 10^{33}` g

**Physical Constants (code units):**
- Speed of light: :math:`c = 1`
- Proton mass: :math:`m_p = 1.67 \times 10^{-24}` g
- Electron mass: :math:`m_e = m_p/1836`
- Thomson cross-section: :math:`\sigma_T = 6.65 \times 10^{-25}` cm²

**Cosmological Parameters:**
- :math:`\Omega_m = 0.27` (matter density)
- :math:`\Omega_\Lambda = 0.73` (dark energy density)
- :math:`H_0 = 67.66` km/s/Mpc (Hubble constant)

Common Unit Conversions
^^^^^^^^^^^^^^^^^^^^^^^

**Distance:**
- 1 Mpc = 3.086 × 10²⁴ cm
- 1 kpc = 3.086 × 10²¹ cm
- 1 AU = 1.5 × 10¹³ cm

**Energy:**
- 1 erg = 1 g⋅cm²/s²
- 1 keV = 1.602 × 10⁻⁹ erg
- 1 GeV = 1.602 × 10⁻³ erg

**Angles:**
- 1 degree = π/180 ≈ 0.01745 radians
- 1 arcminute = π/10800 ≈ 2.91 × 10⁻⁴ radians

Parameter Degeneracies and Fitting Strategies
---------------------------------------------

Understanding parameter correlations helps in MCMC fitting:

**Strong Correlations:**

- :math:`E_{\rm iso}` ↔ :math:`n_{\rm ISM}`: Higher energy can compensate for lower density
- :math:`\epsilon_e` ↔ :math:`\epsilon_B`: Microphysics parameters are often correlated
- :math:`\theta_c` ↔ :math:`\theta_v`: Jet geometry parameters affect observed flux similarly

**Frequency-dependent Constraints:**

- **Radio data**: Most sensitive to :math:`\epsilon_B`, :math:`n_{\rm ISM}`
- **Optical data**: Constrains :math:`\epsilon_e`, :math:`p`, :math:`E_{\rm iso}`
- **X-ray data**: Sensitive to :math:`\Gamma_0`, high-frequency cutoffs

**Time-dependent Constraints:**

- **Early times (< 1 day)**: Constrain :math:`\Gamma_0`, :math:`\epsilon_e`
- **Jet break time**: Determines :math:`\theta_c`, :math:`E_{\rm iso}`
- **Late times (> 100 days)**: Sensitive to :math:`n_{\rm ISM}`, :math:`p`

For more detailed information on parameter estimation strategies and examples of using these parameters in practice, see the :doc:`examples` and :doc:`mcmc_fitting` pages.
