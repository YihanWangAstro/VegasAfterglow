C++ API Reference
================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
-------

VegasAfterglow's C++ core provides high-performance computation capabilities for GRB afterglow modeling. This section documents the C++ API for advanced users and developers who want to work directly with the C++ library.

Key Components
-------------

The C++ API is organized into several core components:

* **Jet Models**: Implementations of different jet structure models (TophatJet, GaussianJet, PowerLawJet)
* **Ambient Medium**: Classes for modeling the circumburst environment (ISM, Wind)
* **Radiation Processes**: Components for calculating synchrotron and inverse Compton emission
* **Dynamics**: Classes for evolving the blast wave and calculating shock parameters
* **Utilities**: Helper functions and mathematical tools

C++ API Documentation
-------------------

Namespaces
---------

.. doxygennamespace:: con
   :project: VegasAfterglow
   :members:
   :outline:

.. doxygennamespace:: math
   :project: VegasAfterglow
   :members:
   :outline:

.. doxygennamespace:: evn
   :project: VegasAfterglow
   :members:
   :outline:

Core Classes
---------

.. doxygenclass:: Medium
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: Ejecta
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

Jet Models
---------

.. doxygenclass:: TophatJet
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: GaussianJet
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: PowerLawJet
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:
   :allow-dot-graphs:

Ambient Medium
------------

.. doxygenclass:: ISM
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: Wind
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

Radiation Processes
-----------------

.. doxygenstruct:: SynPhotons
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenstruct:: SynElectrons
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenstruct:: InverseComptonY
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

Shock Dynamics
------------

.. doxygenclass:: Shock
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: SimpleShockEqn
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: ForwardShockEqn
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: FRShockEqn
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

Physics and Utilities
------------------

.. doxygenclass:: Observer
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

.. doxygenclass:: Coord
   :project: VegasAfterglow
   :members:
   :protected-members:
   :private-members:
   :undoc-members:

Documenting C++ Code
------------------

When contributing to the C++ codebase, please follow these documentation guidelines:

Class and Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Use Doxygen-style comments for all classes and functions:

.. code-block:: cpp

    /********************************************************************************************************************
     * @brief Brief description of the function/class
     * @details Detailed description that provides more information
     *          about what this function/class does, how it works,
     *          and any important details users should know.
     *
     * @param param1 Description of first parameter
     * @param param2 Description of second parameter
     * @return Description of return value
     * @throws Description of exceptions that might be thrown
     * @see RelatedClass, related_function()
     ********************************************************************************************************************/

Member Variable Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For member variables, use inline Doxygen comments with `///<`:

.. code-block:: cpp

    double energy; ///< Isotropic-equivalent energy in ergs
    double gamma0; ///< Initial bulk Lorentz factor

Template Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For template functions, make sure to document both the template parameters and the function parameters:

.. code-block:: cpp

    /********************************************************************************************************************
     * @brief Brief description of the template function
     * @details Detailed description of what the template function does.
     *
     * @tparam T The type of elements in the vector
     * @tparam Comparator The comparison function type
     * @param values Vector of values to be sorted
     * @param comparator Comparator function to determine sorting order
     * @return Sorted vector of values
     ********************************************************************************************************************/
    template<typename T, typename Comparator = std::less<T>>
    std::vector<T> sort_values(const std::vector<T>& values, Comparator comparator = Comparator()) {
        // Implementation details
    }

Inline Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

For inline functions, maintain the same documentation standard as regular functions:

.. code-block:: cpp

    /********************************************************************************************************************
     * @brief Compute the square of a value
     * @details This inline function efficiently computes the square of any numeric value.
     *
     * @param x The value to square
     * @return The squared value
     ********************************************************************************************************************/
    inline double square(double x) {
        return x * x;
    }

Example Class
^^^^^^^^^^^

Here's an example of a well-documented class:

.. code-block:: cpp

    /********************************************************************************************************************
     * @class GaussianJet
     * @brief Implements a Gaussian jet profile where properties follow a Gaussian distribution with angle.
     * @details This class provides a smooth model for GRB jets, characterized by core angle theta_c,
     *          isotropic equivalent energy E_iso, and initial Lorentz factor Gamma0 at the center.
     ********************************************************************************************************************/
    class GaussianJet {
    public:
        /********************************************************************************************************************
         * @brief Constructor: Initialize with core angle, isotropic energy, and initial Lorentz factor
         * @param theta_c Core angle of the jet
         * @param E_iso Isotropic equivalent energy
         * @param Gamma0 Initial Lorentz factor
         ********************************************************************************************************************/
        GaussianJet(Real theta_c, Real E_iso, Real Gamma0) noexcept;
        
        /********************************************************************************************************************
         * @brief Energy per solid angle as a function of phi and theta, with Gaussian falloff
         * @param phi Azimuthal angle (unused)
         * @param theta Polar angle
         * @return Energy per solid angle with Gaussian angular dependence
         ********************************************************************************************************************/
        Real eps_k(Real phi, Real theta) const noexcept;
        
    private:
        Real const norm_{0};    ///< Normalization factor for Gaussian distribution
        Real const eps_k_{0};   ///< Peak energy per solid angle at center
        Real const Gamma0_{1};  ///< Peak Lorentz factor at center
    };

For more details on Doxygen commands, see the :doc:`contributing` page. 