Parameter Reference
===================

This page provides a comprehensive reference for all parameters used in VegasAfterglow, including their physical meanings, typical ranges, and units.

Physical Parameters
-------------------

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
   * - ``theta_v``
     - :math:`\theta_v`
     - radians
     - :math:`0 - \pi/2`
     - Viewing angle (angle between jet axis and line of sight)
   * - ``duration``
     - :math:`T_{\rm dur}`
     - seconds
     - :math:`0.1 - 1000`
     - Duration of energy injection (affects reverse shock)
   * - ``k``
     - :math:`k`
     - dimensionless
     - :math:`1 - 10`
     - Power-law index for structured jets (PowerLawJet only)
   * - ``sigma0``
     - :math:`\sigma_0`
     - dimensionless
     - :math:`0.001 - 10`
     - Initial magnetization parameter

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
   * - ``A_star``
     - :math:`A_*`
     - dimensionless
     - :math:`10^{-3} - 10`
     - Wind parameter: :math:`\rho = A_* \times 5 \times 10^{11} r^{-2}` g/cm³

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

Radiation Microphysics Parameters
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
     - Electron acceleration efficiency

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
   * - ``L_0``
     - :math:`L_0`
     - erg/s
     - :math:`10^{44} - 10^{48}`
     - Initial luminosity of magnetar spin-down
   * - ``t_0``
     - :math:`t_0`
     - seconds
     - :math:`10 - 10^4`
     - Characteristic spin-down timescale
   * - ``q``
     - :math:`q`
     - dimensionless
     - :math:`1 - 6`
     - Power-law index of spin-down: :math:`L(t) = L_0(1+t/t_0)^{-q}`

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
   * - ``theta_n``
     - :math:`\theta_n`
     - radians
     - :math:`0.01 - 0.2`
     - Half-opening angle of narrow component
   * - ``E_iso_n``
     - :math:`E_{\rm iso,n}`
     - erg
     - :math:`10^{51} - 10^{54}`
     - Isotropic energy of narrow component
   * - ``Gamma0_n``
     - :math:`\Gamma_{0,n}`
     - dimensionless
     - :math:`100 - 1000`
     - Initial Lorentz factor of narrow component
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

Computational Parameters
------------------------

Model Resolution
^^^^^^^^^^^^^^^^

.. list-table:: 
   :header-rows: 1
   :widths: 15 15 40

   * - Parameter
     - Units
     - Description
   * - ``phi_ppd``
     - points/degree
     - Angular resolution in azimuthal direction
   * - ``theta_ppd``
     - points/degree
     - Angular resolution in polar direction
   * - ``t_ppd``
     - points/decade
     - Temporal resolution (logarithmic spacing)

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
- For structured jets: :math:`E_{\rm iso} = \int E(\theta) d\Omega` over the jet solid angle

**Causality:**

- Light travel time sets minimum variability timescale: :math:`\delta t \geq R/c\Gamma^2`
- Jet opening angle and Lorentz factor: :math:`\theta_c \gtrsim 1/\Gamma_0` for causal contact

**Microphysics:**

- Energy fractions: :math:`\epsilon_e + \epsilon_B \leq 1` (though often :math:`\ll 1`)
- Electron power-law index: :math:`p > 2` for finite energy in fast-cooling regime

Typical Parameter Combinations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Short GRB (Neutron Star Merger):**

.. code-block:: python

    ParamDef("E_iso",     5e50,  5e52,  Scale.LOG),
    ParamDef("Gamma0",    50,    500,   Scale.LOG),
    ParamDef("theta_c",   0.05,  0.3,   Scale.LINEAR),
    ParamDef("theta_v",   0.1,   0.5,   Scale.LINEAR),
    ParamDef("eps_e",     0.01,  0.3,   Scale.LOG),
    ParamDef("eps_B",     1e-4,  0.1,   Scale.LOG),
    ParamDef("p",         2.1,   2.8,   Scale.LINEAR),

**Long GRB (Collapsar):**

.. code-block:: python

    ParamDef("E_iso",     1e52,  1e54,  Scale.LOG),
    ParamDef("Gamma0",    100,   1000,  Scale.LOG),
    ParamDef("theta_c",   0.02,  0.2,   Scale.LINEAR),
    ParamDef("theta_v",   0.0,   0.1,   Scale.LINEAR),
    ParamDef("eps_e",     0.1,   0.5,   Scale.LOG),
    ParamDef("eps_B",     1e-3,  0.1,   Scale.LOG),
    ParamDef("p",         2.2,   2.6,   Scale.LINEAR),

**Kilonova-associated GRB:**

.. code-block:: python

    ParamDef("E_iso",     1e49,  1e52,  Scale.LOG),
    ParamDef("Gamma0",    10,    300,   Scale.LOG),
    ParamDef("theta_c",   0.1,   0.5,   Scale.LINEAR),
    ParamDef("theta_v",   0.2,   0.8,   Scale.LINEAR),
    ParamDef("eps_e",     0.01,  0.3,   Scale.LOG),
    ParamDef("eps_B",     1e-5,  1e-2,  Scale.LOG),
    ParamDef("p",         2.1,   2.8,   Scale.LINEAR),

Unit Conversions
----------------

Common unit conversions for convenience:

**Distance:**

- 1 Mpc = 3.086 × 10²⁴ cm
- 1 kpc = 3.086 × 10²¹ cm
- Luminosity distance: :math:`d_L = (1+z) \times d_A` (angular diameter distance)

**Energy:**

- 1 BeV = 1.602 × 10⁻³ erg
- 1 keV = 1.602 × 10⁻⁹ erg
- Solar rest mass energy: :math:`M_\odot c^2 = 1.8 \times 10^{54}` erg

**Angles:**

- 1 degree = π/180 ≈ 0.01745 radians
- 1 arcminute = π/10800 ≈ 2.91 × 10⁻⁴ radians
- 1 arcsecond = π/648000 ≈ 4.85 × 10⁻⁶ radians

**Frequencies:**

- X-ray (1 keV): ν ≈ 2.4 × 10¹⁷ Hz
- Optical (V-band): ν ≈ 5.5 × 10¹⁴ Hz  
- Radio (1 GHz): ν = 10⁹ Hz

Parameter Degeneracies
----------------------

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

For more detailed information on parameter estimation strategies, see the :doc:`examples` page.