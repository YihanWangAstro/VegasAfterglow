===========================
GRB Afterglow Physics
===========================

Introduction
------------

Gamma-ray bursts (GRBs) are among the most energetic astrophysical transients in the universe. They are often associated with:

* Mergers of compact objects (e.g., binary neutron stars, neutron star-black hole systems)
* Massive star collapses (collapsar model)

The multi-wavelength afterglow emission provides crucial insights into:

* The structure of relativistic jets
* The surrounding environment
* The underlying shock physics

VegasAfterglow offers a comprehensive modeling framework for GRB afterglows with unprecedented flexibility and computational efficiency. This document outlines the key physical processes implemented in the framework, starting with jet and medium properties, followed by shock dynamics, and finally radiation mechanisms.

Jets and Medium
----------------

The properties of the relativistic outflow and the ambient medium are foundational to any afterglow model. VegasAfterglow employs a flexible approach to represent diverse jet structures and environmental conditions.

VegasAfterglow solves the afterglow emission on a 3D grid of spherical coordinates :math:`(t, \theta, \phi)`, where:

* :math:`t` is the observation time for an on-axis observer
* :math:`\theta` is the polar angle
* :math:`\phi` is the azimuthal angle

The jet is modeled as a relativistic outflow with profiles for:

* Lorentz factor :math:`\Gamma_{0}(\phi, \theta)`
* Energy :math:`\frac{dE}{d\Omega}(\phi, \theta)`
* Magnetization :math:`\sigma(\phi,\theta)`

In the afterglow calculation, the jet is decomposed into a grid of elements on :math:`\theta` and :math:`\phi`, with each jet element evolving independently.

Jet Profiles
~~~~~~~~~~~~

VegasAfterglow supports arbitrary user-defined jet profiles, with several common models built-in:

**Top-hat Jet**
    The simplest jet profile, with constant parameters within a certain opening angle :math:`\theta_{\rm c}` and zero outside:

    .. math::

        \Gamma_{0}(\phi, \theta) = \begin{cases}
            \Gamma_{\rm 0}, & \theta\leq\theta_{\rm c} \\
            1, & \theta>\theta_{\rm c}
        \end{cases}

    .. math::

        \frac{dE}{d\Omega}(\phi, \theta) = \begin{cases}
            \frac{E_{\rm iso}}{4\pi}, & \theta\leq\theta_{\rm c} \\
            0, & \theta>\theta_{\rm c}
        \end{cases}

    .. math::

        \sigma(\phi, \theta) = \begin{cases}
            \sigma_{\rm 0}, & \theta\leq\theta_{\rm c} \\
            0, & \theta>\theta_{\rm c}
        \end{cases}

**Gaussian Jet**
    A jet profile with Gaussian distribution:

    .. math::

        \Gamma_{0}(\phi, \theta) = (\Gamma_{\rm 0}-1)\exp\left(-\frac{\theta^{2}}{2\theta_{\rm c}^{2}}\right)+1

    .. math::

        \frac{dE}{d\Omega}(\phi, \theta) = \frac{E_{\rm iso}}{4\pi}\exp\left(-\frac{\theta^{2}}{2\theta_{\rm c}^{2}}\right)

    .. math::

        \sigma(\phi, \theta) = \sigma_{\rm 0}

**Power-law Jet**
    A jet profile with power-law distribution:

    .. math::

        \Gamma_{0}(\phi, \theta) = (\Gamma_{\rm 0}-1)\left(1+\frac{\theta}{\theta_{\rm c}}\right)^{-k}+1

    .. math::

        \frac{dE}{d\Omega}(\phi, \theta) = \frac{E_{\rm iso}}{4\pi}\left(1+\frac{\theta}{\theta_{\rm c}}\right)^{-k}

    .. math::

        \sigma(\phi, \theta) = \sigma_{\rm 0}

Energy Injection
~~~~~~~~~~~~~~~~

Energy injection is crucial for modeling afterglow emission of GRBs that show plateau features. VegasAfterglow provides arbitrary user-defined energy injection profiles :math:`L_{\rm inj}(\phi, \theta, t)`.

A common model is the magnetar spin-down model, approximated as an isotropic power-law decay function:

.. math::

    \frac{d L_{\rm inj}}{d\Omega}(\phi, \theta, t) = \frac{L_{\rm 0}}{4\pi}\left(1+\frac{t}{t_{\rm 0}}\right)^{-q}

.. math::

    \sigma(\phi, \theta) = \infty

where :math:`\sigma = \infty` assumes an ideal Poynting-flux dominated wind with no matter injection.

Ambient Medium
~~~~~~~~~~~~~~

The properties of the circumburst environment significantly influence afterglow evolution. The ambient medium can be user-defined through functions :math:`m(\phi,\theta, r)` and :math:`\rho(\phi, \theta, r)`. Two standard models are:

**Homogeneous Medium (ISM)**

.. math::

    \rho(\phi, \theta, r) = n_{0} m_{\rm p}

where :math:`n_{0}` is the number density of the interstellar medium.

**Wind-like Medium**

.. math::

    \rho(\phi, \theta, r) = \frac{A}{r^k}

where :math:`A` is related to the mass-loss rate and :math:`k` is the wind profile index (typically 2 for a stellar wind).

Shock Dynamics
--------------

As the relativistic jet propagates through the ambient medium, it drives a shock wave that accelerates particles and amplifies magnetic fields. Understanding the dynamics of these shocks is crucial for accurate afterglow modeling.

Shock Jump Conditions
~~~~~~~~~~~~~~~~~~~~~

The fundamental physical principles governing the behavior of relativistic shocks are formulated through the shock jump conditions, which relate quantities across the shock discontinuity.

In the rest frame of the shock, the continuity equation is:

.. math::

    n_{\rm u} u_{\rm us} = n_{\rm d} u_{\rm ds}

where:

* :math:`n_{\rm u}` and :math:`n_{\rm d}` are the proton number density in the upstream and downstream
* :math:`u_{\rm us}` and :math:`u_{\rm ds}` are the four-velocities

Combined with magnetic flux conservation, energy conservation, and enthalpy conservation across the shock front, one can derive:

.. math::

    u_{\rm us} = u_{\rm ds}\Gamma_{\rm ud}+\sqrt{(u_{\rm ds}^{2}+1)(\Gamma_{\rm ud}^{2}-1)}

.. math::

    Au_{\rm ds}^6 + Bu_{\rm ds}^4 +Cu_{\rm ds}^2+D = 0

where:

.. math::

    A &= \hat\gamma_{\rm ud}(2-\hat\gamma_{\rm ud})(\Gamma_{\rm ud}-1)+2\\
    B &= -(\Gamma_{\rm ud}+1)[(2-\hat\gamma_{\rm ud})(\hat\gamma_{\rm ud}\Gamma_{\rm ud}^2+1)+\hat\gamma_{\rm ud}(\hat\gamma_{\rm ud}-1)\Gamma_{\rm ud}]\sigma_{\rm u}\\
    &-(\Gamma_{\rm ud}-1)[\hat\gamma_{\rm ud}(2-\hat\gamma_{\rm ud})(\Gamma_{\rm ud}^2-2)+(2\Gamma_{\rm ud}+3)]\\
    C &= (\Gamma_{\rm ud}+1)[\hat\gamma_{\rm ud}(1-\hat\gamma_{\rm ud}/4)(\Gamma_{\rm ud}^2-1)+1]\sigma_{\rm u}^2\\
    &+(\Gamma_{\rm ud}^2-1)[2\hat\gamma_{\rm ud}-(2-\hat\gamma_{\rm ud})(\hat\gamma_{\rm ud}\Gamma_{\rm ud}-1)]\sigma_{\rm u}\\
    &+(\Gamma_{\rm ud}+1)(\Gamma_{\rm ud}-1)^2(\hat\gamma_{\rm ud}-1)^2\\
    D &= -(\Gamma_{\rm ud}-1)(\Gamma_{\rm ud}+1)^2(2-\hat\gamma_{\rm ud})^2\sigma_{\rm u}^2/4\\
    \hat\gamma_{\rm ud} &= \frac{4\Gamma_{\rm ud}+1}{3\Gamma_{\rm ud}}

The solution for :math:`u_{\rm ds}^2` can be derived as:

.. math::

    u_{\rm ds}^2 = 2\sqrt{\frac{-P}{3}}\cos\left(\frac{1}{3}\arccos\left(\frac{3Q}{2P}\sqrt{\frac{-3}{P}}\right)-\frac{2\pi}{3}\right)-\frac{\mathcal{B}}{3}

where:

.. math::

    \mathcal{B} &= \frac{B}{A}, \mathcal{C} = \frac{C}{A}, \mathcal{D} = \frac{D}{A}\\
    P &= \mathcal{C} - \frac{\mathcal{B}^{2}}{3}\\
    Q &= \frac{2\mathcal{B}^{3}}{27}-\frac{\mathcal{BC}}{3}+\mathcal{D}

For unmagnetized upstream where :math:`\sigma_{\rm u}=0`, the solution simplifies to:

.. math::

    u_{\rm ds}^2 = \frac{(\Gamma_{\rm ud}-1)(\hat\gamma_{\rm ud}-1)^{2}}{\hat\gamma_{\rm ud}(2-\hat\gamma_{\rm ud})(\Gamma_{\rm ud}-1)+2}

.. math::

    u_{\rm us}^2 = \frac{(\Gamma_{\rm ud}-1)(\hat\gamma_{\rm ud}\Gamma_{\rm ud}+1)^{2}}{\hat\gamma_{\rm ud}(2-\hat\gamma_{\rm ud})(\Gamma_{\rm ud}-1)+2}

and the density compression ratio becomes:

.. math::

    \frac{n_{\rm d}}{n_{\rm u}} = \frac{\hat\gamma_{\rm ud}\Gamma_{\rm ud}+1}{\hat\gamma_{\rm ud}-1} = 4\Gamma_{\rm ud}

The internal energy density generated in the downstream due to shock heating is:

.. math::

    e_{\rm d} = (\Gamma_{\rm ud}-1)n_{\rm d} m_{\rm p} c^2

The downstream co-moving magnetic field is given by:

.. math::

    \vec B_{\rm d} = \vec B_{\rm d, o} + \vec B_{\rm d, w} =
    \frac{u_{\rm us}}{u_{\rm ds}}\vec B_{\rm u,o} +
    \sqrt{8\pi\epsilon_{\rm B}e_{\rm d}}\frac{\vec B_{\rm d, w}}{B_{\rm d, w}}

where :math:`B_{\rm d, o}` is the ordered magnetic field and :math:`B_{\rm d, w}` is the unordered magnetic field generated by Weibel instability, with :math:`\epsilon_{\rm B}` being the magnetic energy fraction.

The general magnetic energy fraction can be approximated as:

.. math::

    \bar{\epsilon}_{\rm B} = \frac{B_{\rm d}^{2}}{8\pi e_{\rm d}} \sim
    \frac{\sigma_{\rm u}}{3(1+\sigma_{\rm u})} + \epsilon_{\rm B}

Shock Lorentz Factor
~~~~~~~~~~~~~~~~~~~~

The evolution of the shock's Lorentz factor determines the energy dissipation rate and consequently the afterglow light curve properties. For a given relative Lorentz factor :math:`\Gamma_{\rm ud}` between upstream and downstream, the shock jump conditions provide the downstream proton number density :math:`n_{\rm d}` and co-moving magnetic field :math:`B_{\rm d}` required for radiation calculations.

The relative Lorentz factor is given by:

.. math::

    \Gamma_{\rm ud} = \Gamma_{\rm u}\Gamma_{\rm d} - \sqrt{(\Gamma_{\rm u}^2-1)(\Gamma_{\rm d}^2-1)}

where :math:`\Gamma_{\rm u}` and :math:`\Gamma_{\rm d}` are the Lorentz factors of the upstream and downstream regions.

Reverse Shock Crossing
~~~~~~~~~~~~~~~~~~~~~~

When the relativistic jet encounters the ambient medium, a two-shock structure forms: a forward shock propagating into the medium and a reverse shock moving back into the jet material. This phase is crucial for early afterglow emission.

As the jet collides with the ambient medium, a forward shock and reverse shock pair is generated (if reverse shock conditions are met). These shocks divide the system into four regions:

1. Unshocked medium (upstream of forward shock)
2. Shocked medium (downstream of forward shock)
3. Shocked jet (downstream of reverse shock)
4. Unshocked jet (upstream of reverse shock)

During the reverse shock crossing phase, we need to combine the shock jump conditions of forward and reverse shocks. For the forward shock:

.. math::

    n_2 &= 4\Gamma_{12}n_1 \\
    e_2 &= (\Gamma_{12}-1)n_2m_{\rm p}c^2

where :math:`n_1` is the ambient medium density and :math:`\Gamma_{12}` is the relative Lorentz factor between regions 1 and 2.

For the reverse shock:

.. math::

    n_3 &= \frac{u_{4s}}{u_{3s}}n_4 = \left(\Gamma_{43}+\frac{\sqrt{(u_{\rm 3s}^{2}+1)(\Gamma_{\rm 43}^{2}-1)}}{u_{3s}}\right)n_4 \\
    e_3 &= (\Gamma_{43}-1)n_3m_{\rm p} c^2 \\
    {B_{3}^{2}} &= \frac{u_{4s}^{2}}{u_{3s}^{2}}{B_{4}^{2}}

where:

.. math::

    n_4 &= \frac{E_{\rm jet}}{4\pi r^2\Delta^\prime\Gamma_{4}m_{p}c^{2}(1+\sigma_{4})} \\
    {B_{4}^{2}} &= 4\pi\sigma_4 n_4 m_{\rm p}c^2

with :math:`\Delta^\prime` being the shock width in the co-moving frame.

The force balance at the discontinuity gives:

.. math::

    p_2 = (\hat\gamma_{12}-1)e_2=(\hat\gamma_{43}-1)e_3+ \frac{B_{3}^{2}}{8\pi} = p_3

Combining these equations, we get:

.. math::

    1 = \frac{p_2}{p_3} \sim \frac{n_1(\Gamma_{12}^2-1)}{n_4(1+\sigma_4)(\Gamma_{43}^2-1)}

Since :math:`\Gamma_1`, :math:`\Gamma_4`, :math:`n_1` and :math:`n_4` are known, and :math:`\Gamma_2=\Gamma_3`, we can solve for :math:`\Gamma_2` during the reverse shock crossing phase.

If the pressure in region 4 (unshocked jet) is too strong, the reverse shock cannot form. The reverse shock generation condition is:

.. math::

    \sigma_4 &\lesssim 8(\hat\gamma_{12}-1)\Gamma_{12}(\Gamma_{12}-1)\frac{n_{1}}{n_{4}} \\
    &= \frac{8}{3}(\Gamma_{12}^2-1)\frac{n_{1}}{n_{4}} \\
    &\sim \frac{8}{3}(\Gamma_{4}^2-1)\frac{n_{1}}{n_{4}}

For energy conservation during long-lasting reverse shock crossings, VegasAfterglow implements an effective mechanical model. The total energy as the blast wave propagates to radius :math:`r` should be:

.. math::

    E_{\rm jet}+E_{\rm medium}=N_{4,0}\Gamma_4m_pc^2(1+\sigma_4)+N_2m_pc^2

The total energy in regions 2, 3, and 4 is:

.. math::

    E_2+E_3+E_4 &= N_2[\Gamma_2+\Gamma_{\rm eff,2}(\Gamma_{21}-1)]m_pc^2\\
    &+N_3[\Gamma_3+\Gamma_{\rm eff,3}(\Gamma_{43}-1)]m_pc^2(1+\sigma_4)\\
    &+ N_4\Gamma_4m_pc^2(1+\sigma_4)

where:

.. math::

    \Gamma_{\rm eff,2} &= \frac{\hat\gamma_{21}\Gamma_2^2-\hat\gamma_{21}+1}{\Gamma_2}\\
    \Gamma_{\rm eff,3} &= \frac{\hat\gamma_{43}\Gamma_3^2-\hat\gamma_{43}+1}{\Gamma_3}

Using energy conservation and particle number conservation:

.. math::

    E_{\rm jet}+E_{\rm medium} &= E_2+E_3+E_4\\
    N_{4,0} &= N_3+N_4

And with :math:`\Gamma_3=\Gamma_2=\Gamma_{21}`, we get:

.. math::

    0 &= N_2[\Gamma_3-1+\Gamma_{\rm eff,2}(\Gamma_{3}-1)]\\
    &+N_3(1+\sigma_4)[\Gamma_3-\Gamma_4+\Gamma_{\rm eff,3}(\Gamma_{43}-1)]

The shocked proton number :math:`N_3` can be expressed as:

.. math::

    N_{3}(r) &= \int_{0}^{r} 4\pi r^{\prime 2}n_3\frac{d\Delta_3^\prime}{dr^{\prime}}dr^\prime\\
    \frac{d\Delta_3^\prime}{dr} &= \frac{\Gamma_3}{\Gamma_{4}}\left(\frac{n_{1}}{n_{4}(1+\sigma_{4})}\right)^{1/2}\left(1-\frac{\Gamma_{4}n_{4}}{\Gamma_{3}n_{3}}\right)^{-1}

The reverse shock completely crosses the jet when:

.. math::

    N_{3}(r_x) = N_{4} \equiv \frac{E_{\rm jet}}{\Gamma_{4}m_{p}c^{2}(1+\sigma_{4})}

where :math:`r_x` is the crossing radius.

Post-Reverse Shock Crossing Phase
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once the reverse shock completely traverses the jet material, the dynamics enter a new regime characterized by adiabatic expansion and gradual deceleration. This phase is essential for modeling the mid-to-late afterglow emission.

After the reverse shock has fully crossed the jet (:math:`N_3 = N_{4,0}`), the shocked jet material enters the post-crossing phase. In this phase, the dynamics follow a self-similar solution that smoothly transitions from the Newtonian regime to the relativistic Blandford-McKee solution.

The dynamics for both the forward and reverse shocks need to be treated separately after crossing.

**Forward Shock Dynamics**

The forward shock continues to propagate into the ambient medium and is governed by energy conservation in the blast wave. VegasAfterglow carefully accounts for pressure, radiation, and adiabatic expansion:

.. math::

    {d}[\Gamma_2M_{\rm tot}c^2+\Gamma_{\rm eff}U] &= c^2{dm}+\Gamma_{\rm eff}dU_{\rm rad}+L_{\rm inj}dt_{\rm eng}\\
    {dU} &= {dU_{\rm sh}}+dU_{\rm ad}+{dU_{\rm rad}}\\
    M_{\rm tot} &= M_0 + m + m_{\rm inj}

where :math:`E_{\rm jet}=\Gamma_{2}M_{0}c^{2}` with :math:`M_{0}` the equivalent initial rest mass of the jet, :math:`m` is the total collected mass, :math:`m_{\rm inj}` is the injected mass, :math:`U` is the internal energy in the blast wave, :math:`dU_{\rm sh}` is increased internal energy due to shock heating, :math:`dU_{\rm ad}` is the lost internal energy due to adiabatic expansion, :math:`dU_{\rm rad}` is lost internal energy due to radiation, and :math:`L_{\rm inj}` is the injection luminosity.

These quantities are defined as:

.. math::

    \Gamma_{\rm eff} &= \frac{\hat{\gamma}_{2}\Gamma_{2}^{2}-\hat{\gamma}_{2}+1}{\Gamma_{2}}\\
    dU_{\rm sh} &= c^2 (\Gamma_2-1)dm\\
    dU_{\rm ad} &= -(\hat{\gamma}_2-1)\bigg(\frac{3}{r}dr-\frac{1}{\Gamma_{2}}d\Gamma_2\bigg)U\\
    dU_{\rm rad} &= -\epsilon(\Gamma_2-1)c^2dm\\
    \epsilon &= \epsilon_{\rm rad}\epsilon_{e}

The blast wave equation for :math:`\Gamma_2` evolution can then be written as:

.. math::

    \dot\Gamma_{2} = \frac{-(\Gamma_{\rm eff}+1)(\Gamma_2-1)c^2\dot{m}-\Gamma_{\rm eff}\dot U_{\rm ad}+L_{\rm inj} - c^2\Gamma_{2}\dot m_{\rm inj}}{(M_{0}+ m +m_{\rm inj})c^{2}+U\frac{d\Gamma_{\rm eff}}{d\Gamma_2}}

**Reverse Shock Dynamics**

For the reverse shock, the four-velocity evolution in the post-crossing phase follows a power-law profile:

.. math::

    u = u_x \left(\frac{r}{r_x}\right)^{-g}

where :math:`u_x` is the four-velocity at the crossing radius :math:`r_x`, and :math:`g` is the power-law index that depends on the relative Lorentz factor :math:`\Gamma_{\rm rel}`.

The scaling laws after reverse shock crossing vary between thin and thick shell regimes:

1. **Thin Shell Regime** (Newtonian reverse shock):

   .. math::

       n_3 &\propto r^{-\frac{2(3+g)}{\hat\gamma_{43}+1}}\\
       p_3 &\propto r^{-\frac{2(3+g)\hat\gamma_{43}}{\hat\gamma_{43}+1}}\\
       \frac{p_3}{n_3} &\propto r^{-\frac{2(3+g)(\hat\gamma_{43}-1)}{\hat\gamma_{43}+1}}

2. **Thick Shell Regime** (Relativistic reverse shock):

   Before shell expansion:

   .. math::

       n_3 &\propto r^{-2}\\
       p_3 &\propto r^{-2\hat\gamma_{43}}\\
       \frac{p_3}{n_3} &\propto r^{2(1-\hat\gamma_{43})}

   After shell expansion:

   .. math::

       n_3 &\propto r^{-(3+g)}\\
       p_3 &\propto r^{-(3+g)\hat\gamma_{43}}\\
       \frac{p_3}{n_3} &\propto r^{(3+g)(1-\hat\gamma_{43})}

The power-law index :math:`g` transitions smoothly between different regimes:

.. math::

    g = g_{\rm low} + (g_{\rm high} - g_{\rm low})\frac{p}{1+p}

with:

.. math::

    p = \sqrt{\sqrt{\Gamma_{\rm rel} - 1}}

For a homogeneous medium (ISM), :math:`g_{\rm low} = 1.5` in the Newtonian regime, and :math:`g_{\rm high} = 3.5` in the ultra-relativistic regime (Blandford-McKee limit).

The shell width evolution is also important for determining the dynamics:

.. math::

    \Delta^\prime = \Delta^\prime_0+ c_s t^\prime

where :math:`c_s` is the sound speed in the co-moving frame:

.. math::

    c_s = c\sqrt{\frac{\hat\gamma_{43}(\hat\gamma_{43}-1)(\Gamma_{43}-1)}{1+\hat\gamma_{43}(\Gamma_{43}-1)}}

The post-crossing evolution preserves the total number of electrons (particle number conservation):

.. math::

    n_3 r^2 \Delta^\prime \propto r^0

which leads to:

.. math::

    n_3(r) = n_{3}(r_x)\left(\frac{r}{r_x}\right)^{-2}\left(\frac{\Delta^\prime}{\Delta_x^\prime}\right)^{-1}

The pressure evolves according to adiabatic expansion:

.. math::

    p_3 \propto n_3^{\hat{\gamma}_{43}}

which gives:

.. math::

    p_3(r) = p_{3}(r_x)\left(\frac{r}{r_x}\right)^{-2\hat{\gamma}_{43}}\left(\frac{\Delta^\prime}{\Delta_x^\prime}\right)^{-\hat{\gamma}_{43}}

The average electron Lorentz factor in region 3 (used for synchrotron radiation calculation) is:

.. math::

    \bar\gamma_3(r) \propto \frac{p_3}{n_3} \propto \left(\frac{r}{r_x}\right)^{2(1-\hat{\gamma}_{43})}\left(\frac{\Delta^\prime}{\Delta_x^\prime}\right)^{1-\hat{\gamma}_{43}}

This detailed treatment of post-crossing dynamics allows VegasAfterglow to accurately model the afterglow emission from the earliest phases through the transition to the self-similar regime and into the non-relativistic phase.

Jet Spreading
~~~~~~~~~~~~~

As the jet propagates and decelerates, it also expands laterally. This lateral spreading significantly affects the observed light curve, particularly for off-axis observers.

VegasAfterglow models lateral expansion of the jet through:

.. math::

    \frac{d\theta}{dt} = F(u)\sqrt{\frac{2u^2+3}{4u^2+3}}\frac{1}{2\Gamma}\frac{dr}{dt}

.. math::

    F(u) = \left[1+\left(\frac{u}{Q\theta_s}\right)^2\right]^{-1}

where :math:`u = \Gamma\beta` and :math:`\frac{dr}{dt} = \beta c/(1-\beta)`.

Radiation Processes
-------------------

The observed afterglow emission arises from various radiation mechanisms as relativistic particles interact with magnetic fields and ambient photons. VegasAfterglow implements detailed treatments of these processes to generate accurate multi-wavelength light curves.

Electron Population
~~~~~~~~~~~~~~~~~~~

Particle acceleration at relativistic shocks produces a non-thermal distribution of electrons, which are the primary emitters in the afterglow.

In relativistic shocks, the dissipation of kinetic energy accelerates particles to relativistic energies. The comoving internal energy density behind the shock is:

.. math::

    \mathcal{U} = (\Gamma-1)n_pm_pc^2

where :math:`n_p` is the downstream proton number density.

A fraction of this energy goes to magnetic fields:

.. math::

    \frac{B^{\prime2}}{8\pi} = \epsilon_B\mathcal{U}

And to non-thermal electrons:

.. math::

    (\bar\gamma-1)n_em_ec^2 = \epsilon_e\mathcal{U} = \epsilon_en_p(\Gamma-1)m_pc^2

where :math:`\bar\gamma` is the average Lorentz factor of the accelerated electrons.

The electron energy distribution typically follows a power law:

.. math::

    N(\gamma) \propto \gamma^{-p}, \quad \gamma_m \leq \gamma \leq \gamma_M

where :math:`\gamma_m` and :math:`\gamma_M` are the minimum and maximum Lorentz factors.

The minimum Lorentz factor :math:`\gamma_m` can be determined from the energy equations:

.. math::

    \gamma_m= \begin{cases}
    \frac{p-2}{p-1}\frac{\epsilon_e}{\xi}(\Gamma-1)\frac{m_p}{m_e}+1, & p>2\\
    \ln^{-1}(\frac{\gamma_M}{\gamma_m})\frac{\epsilon_e}{\xi}(\Gamma-1)\frac{m_p}{m_e}+1, & p=2\\
    \left(\frac{2-p}{p-1}\frac{\epsilon_e}{\xi}(\Gamma-1)\frac{m_p}{m_e}\gamma_M^{p-2}\right)^{1/(p-1)}+1, & p<2
    \end{cases}

where :math:`\xi = \frac{n_e}{n_p}` is the electron-to-proton number ratio.

The cooling of electrons due to synchrotron radiation and inverse Compton scattering is described by:

.. math::

    \frac{d\gamma}{dt^\prime}m_ec^2 = -P_{\rm syn} = -\frac{4}{3}\sigma_Tc\gamma^2\beta^2\frac{B^{\prime2}}{8\pi}(1+\tilde{Y})

where :math:`\tilde{Y}` is the Compton Y-parameter that accounts for the inverse Compton cooling.

The cooling time is:

.. math::

    t_c^\prime = \left|\frac{\gamma_c}{\dot{\gamma_c}}\right| = \frac{6\pi m_ec^2\gamma_c}{\sigma_Tc(\gamma_c^2-1)B^{\prime2}(1+\tilde{Y})}

The cooling Lorentz factor :math:`\gamma_c` is then:

.. math::

    \gamma_c &= \frac{1}{2}\left(\bar{\gamma}_c+\sqrt{\bar{\gamma}^2_c+4}\right)\\
    \bar{\gamma}_c &= \frac{6\pi m_ec}{\sigma_TB^{\prime2}(1+\tilde{Y})t^\prime}

The self-absorption frequency :math:`\nu_a` defines where synchrotron photons are self-absorbed. It can be derived as the intersection point between synchrotron and blackbody spectra:

.. math::

    I_{\rm \nu}^{\rm bb}(\nu_a) &= I_{\nu}^{\rm syn}(\nu_a) \sim 2kT\frac{\nu_a^2}{c^2}\\
    kT &= (\hat\gamma_{p}-1)(\gamma_{p}-1)m_ec^2

The absorption Lorentz factor :math:`\gamma_a` is given by:

.. math::

    \gamma_a = \begin{cases}
    \left(\frac{I_{\nu,p}c^2}{2kT\nu_{p}^{1/3}}\right)^{3/5}, & \nu_{p}<\nu_a\\
    \left(\frac{I_{\nu,p}}{2m_e(\hat\gamma_{p}-1)}\sqrt{\frac{3eB^\prime}{4\pi m_ec}}\right)^{2/5}, & \nu_{p}=\nu_a
    \end{cases}

where:

.. math::

    I_{\rm\nu,p}^{\rm syn} = P_{\rm syn, p}\frac{f_{\rm syn}N_e\xi}{4\pi \delta^\prime} = \frac{\sqrt{3}}{2}\frac{(p-1)B^\prime e^3}{m_ec^2}\frac{f_{\rm syn}\Sigma_e\xi}{4\pi}

The electron distribution function depends on the ordering of characteristic Lorentz factors. For slow cooling with weak absorption (:math:`\gamma_a<\gamma_m<\gamma_c`):

.. math::

    \frac{dN_e^{\rm syn}}{d\gamma}= \begin{cases}
    n_0(p-1)\gamma_m^{p-1}\gamma^{-p}, & \gamma_m<\gamma<\gamma_c\\
    n_0(p-1)\gamma_m^{p-1}\gamma_c\gamma^{-p-1}, & \gamma_c<\gamma
    \end{cases}

For fast cooling with weak absorption (:math:`\gamma_a<\gamma_c<\gamma_m`):

.. math::

    \frac{dN_e^{\rm syn}}{d\gamma}= \begin{cases}
    n_0\gamma_c\gamma^{-2}, & \gamma_c<\gamma<\gamma_m\\
    n_0\gamma_m^{p-1}\gamma_c\gamma^{-p-1}, & \gamma_m<\gamma
    \end{cases}

Synchrotron Radiation
~~~~~~~~~~~~~~~~~~~~~

The characteristic synchrotron emission frequency is:

.. math::

    \nu = \frac{3eB^\prime\gamma^2}{4\pi m_e c}

The synchrotron spectrum is approximated as a multi-segment broken power law, with different forms depending on the ordering of key frequencies:

* :math:`\nu_a` - self-absorption frequency
* :math:`\nu_m` - characteristic frequency of minimum energy electrons
* :math:`\nu_c` - cooling frequency

These define different regimes (e.g., slow cooling with :math:`\nu_m < \nu_c` or fast cooling with :math:`\nu_c < \nu_m`).

Below are the functional forms for various spectral regimes:

**(I) Slow cooling with weak absorption** (:math:`\nu_a<\nu_m<\nu_c`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu_a}{\nu_m}\right)^{1/3}\left(\frac{\nu}{\nu_a}\right)^{2}, & \nu<\nu_a\\
    \left(\frac{\nu}{\nu_m}\right)^{1/3}, & \nu_a<\nu<\nu_m\\
    \left(\frac{\nu}{\nu_m}\right)^{-(p-1)/2}, & \nu_m<\nu<\nu_c\\
    \left(\frac{\nu_c}{\nu_m}\right)^{-(p-1)/2}\left(\frac{\nu}{\nu_c}\right)^{-p/2}, & \nu_c<\nu<\nu_M
    \end{cases}

**(II) Slow cooling with weak absorption** (:math:`\nu_m<\nu_a<\nu_c`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu_m}{\nu_a}\right)^{(p+4)/2}\left(\frac{\nu}{\nu_m}\right)^{2}, & \nu<\nu_m\\
    \left(\frac{\nu_a}{\nu_m}\right)^{-(p-1)/2}\left(\frac{\nu}{\nu_a}\right)^{5/2}, & \nu_m<\nu<\nu_a\\
    \left(\frac{\nu}{\nu_m}\right)^{-(p-1)/2}, & \nu_a<\nu<\nu_c\\
    \left(\frac{\nu_c}{\nu_m}\right)^{-(p-1)/2}\left(\frac{\nu}{\nu_c}\right)^{-p/2}, & \nu_c<\nu<\nu_M
    \end{cases}

**(III) Fast cooling with weak absorption** (:math:`\nu_a<\nu_c<\nu_m`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu_a}{\nu_c}\right)^{1/3}\left(\frac{\nu}{\nu_a}\right)^{2}, & \nu<\nu_a\\
    \left(\frac{\nu}{\nu_c}\right)^{1/3}, & \nu_a<\nu<\nu_c\\
    \left(\frac{\nu}{\nu_c}\right)^{-1/2}, & \nu_c<\nu<\nu_m\\
    \left(\frac{\nu_m}{\nu_c}\right)^{-1/2}\left(\frac{\nu}{\nu_m}\right)^{-p/2}, & \nu_m<\nu<\nu_M
    \end{cases}

**(IV) Fast cooling with strong absorption** (:math:`\nu_c<\nu_a<\nu_m`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu}{\nu_a}\right)^{2}, & \nu<\nu_a\\
    \mathcal{R}_4\left(\frac{\nu}{\nu_a}\right)^{-1/2}, & \nu_a<\nu<\nu_m\\
    \mathcal{R}_4\left(\frac{\nu_m}{\nu_a}\right)^{-1/2}\left(\frac{\nu}{\nu_m}\right)^{-p/2}, & \nu_m<\nu<\nu_M
    \end{cases}

where :math:`\mathcal{R}_4=\frac{\gamma_c}{3\gamma_a}`

**(V) Slow cooling with strong absorption** (:math:`\nu_m<\nu_c<\nu_a`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu}{\nu_a}\right)^{2}, & \nu<\nu_a\\
    \mathcal{R}_5\left(\frac{\nu}{\nu_a}\right)^{-p/2}, & \nu_a<\nu<\nu_M
    \end{cases}

where :math:`\mathcal{R}_5=(p-1)\frac{\gamma_c}{3\gamma_a}\left(\frac{\gamma_m}{\gamma_a}\right)^{p-1}`

**(VI) Fast cooling with strong absorption** (:math:`\nu_c<\nu_m<\nu_a`):

.. math::

    F_\nu=F_{\rm\nu,max}\begin{cases}
    \left(\frac{\nu}{\nu_a}\right)^{2}, & \nu<\nu_a\\
    \mathcal{R}_6\left(\frac{\nu}{\nu_a}\right)^{-p/2}, & \nu_a<\nu<\nu_M
    \end{cases}

where :math:`\mathcal{R}_6=\frac{\gamma_c}{3\gamma_a}\left(\frac{\gamma_m}{\gamma_a}\right)^{p-1}`

Inverse Compton Process
~~~~~~~~~~~~~~~~~~~~~~~

Inverse Compton (IC) scattering is modeled through:

.. math::

    I_\nu^{\rm IC} = \delta^\prime\int_{\gamma_m}^\infty\frac{dN}{d\gamma}(\gamma)d\gamma\int_0^{x_0}\sigma_c(x)I_{\nu}^{\rm syn}(x)dx

where :math:`\sigma_c` is the Klein-Nishina cross-section, accounting for the quantum effects at high energies:

.. math::

    \sigma_c &= \frac{3}{4}\sigma_T\left[\frac{1+x}{x^3}\left\{\frac{2x(1+x)}{1+2x}-\ln(1+2x)\right\}\right.\\
    &\left.+\frac{1}{2x}\ln(1+2x)-\frac{1+3x}{(1+2x)^2} \right]\\
    x &= \frac{h\nu}{m_ec^2}

The IC cooling affects the electron distribution. Defining :math:`\hat\gamma_i=\frac{m_ec^2}{h\nu_i}`, the modified electron distribution due to IC cooling is:

.. math::

    \frac{dN^{\rm IC}_e}{d\gamma}(\gamma)=\frac{dN^{\rm syn}_e}{d\gamma}(\gamma)\begin{cases}
    1, & \gamma<\gamma_c\\
    \frac{1+\tilde{Y}(\gamma_c)}{1+\tilde{Y}(\gamma)}, & \gamma>\gamma_c
    \end{cases}

For :math:`\hat\gamma_m<\hat\gamma_c` (corresponding to :math:`\nu_m > \nu_c`):

.. math::

    \tilde{Y}(\gamma)=\begin{cases}
    Y_T, & \gamma<\hat\gamma_m\\
    Y_T\left(\frac{\gamma}{\hat\gamma_m}\right)^{-1/2}, & \hat\gamma_m<\gamma<\hat\gamma_c\\
    Y_T\left(\frac{\gamma}{\gamma_c}\right)^{-4/3}\left(\frac{\hat\gamma_m}{\hat\gamma_c}\right)^{1/2}, & \hat\gamma_c<\gamma
    \end{cases}

For :math:`\hat\gamma_c<\hat\gamma_m` (corresponding to :math:`\nu_c > \nu_m`):

.. math::

    \tilde{Y}(\gamma)=\begin{cases}
    Y_T, & \gamma<\hat\gamma_c\\
    Y_T\left(\frac{\gamma}{\hat\gamma_c}\right)^{(p-3)/2}, & \hat\gamma_c<\gamma<\hat\gamma_m\\
    Y_T\left(\frac{\gamma}{\gamma_m}\right)^{-4/3}\left(\frac{\hat\gamma_m}{\hat\gamma_c}\right)^{(p-3)/2}, & \hat\gamma_c<\gamma
    \end{cases}

The synchrotron spectrum is also modified by IC cooling:

.. math::

    I_\nu^{\rm syn, IC}(\nu)=I_\nu^{\rm syn}(\nu)\begin{cases}
    1, & \nu<\nu_c\\
    \frac{1+\tilde{Y}(\nu_c)}{1+\tilde{Y}(\nu)}, & \nu>\nu_c
    \end{cases}

For :math:`\hat\nu_m<\hat\nu_c`:

.. math::

    \tilde{Y}(\nu)=\begin{cases}
    Y_T, & \nu<\hat\nu_m\\
    Y_T\left(\frac{\nu}{\hat\nu_m}\right)^{-1/4}, & \hat\nu_m<\nu<\hat\nu_c\\
    Y_T\left(\frac{\nu}{\nu_c}\right)^{-2/3}\left(\frac{\hat\nu_m}{\hat\nu_c}\right)^{1/4}, & \hat\nu_c<\nu
    \end{cases}

For :math:`\hat\nu_c<\hat\nu_m`:

.. math::

    \tilde{Y}(\nu)=\begin{cases}
    Y_T, & \nu<\hat\nu_c\\
    Y_T\left(\frac{\nu}{\hat\nu_c}\right)^{(p-3)/4}, & \hat\nu_c<\nu<\hat\nu_m\\
    Y_T\left(\frac{\nu}{\nu_m}\right)^{-2/3}\left(\frac{\hat\nu_m}{\hat\nu_c}\right)^{(p-3)/4}, & \hat\nu_c<\nu
    \end{cases}

Observed Emission
-----------------

The final step in afterglow modeling is translating the emission properties in the comoving frame to the observer frame, accounting for relativistic effects and geometry.

The observed emission depends on geometry and relativistic effects. For a given observer angle, the Doppler factor is:

.. math::

    \mathcal{D} = \frac{1}{\Gamma(1-\beta\cos w)}

where :math:`\cos w = \sin\theta\cos\phi\sin\theta_v+\cos\theta\cos\theta_v`.

The observed time is:

.. math::

    t_{\rm obs} = (1+z)\left[t+\frac{r}{c}\cos w\right]

And the observed flux:

.. math::

    F_\nu(\nu,\tilde{t}) = \frac{1+z}{4\pi D_L^2}\iint_{t_{\rm obs}=\tilde{t}} 4\pi r^2 I_{\nu^\prime}(\nu^\prime)\mathcal{D}^3d\Omega dr

where :math:`\nu^\prime = \frac{1+z}{\mathcal{D}}\nu`.

References
----------

For detailed exploration of GRB afterglow physics, see:
