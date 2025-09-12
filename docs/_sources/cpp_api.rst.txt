C++ API Reference
=================

.. contents:: Table of Contents
   :local:
   :depth: 2

Overview
--------

VegasAfterglow's C++ core provides high-performance computation capabilities for GRB afterglow modeling. This section documents the C++ API for advanced users and developers who want to work directly with the C++ library. The library requires a C++20 compatible compiler.

Key Components
--------------

The C++ API is organized into several core components:

* **Jet Models**: Implementations of different jet structure models (top-hat, Gaussian, power-law, and user-defined)
* **Ambient Medium**: Classes for modeling the circumburst environment (ISM, wind, and user-defined)
* **Radiation Processes**: Components for synchrotron emission, synchrotron self-absorption, and inverse Compton (with Klein-Nishina corrections)
* **Shock Dynamics**: Forward and reverse shock physics with accurate evolution across relativistic and non-relativistic regimes
* **Observer Effects**: Classes for handling off-axis viewing angles and relativistic beaming
* **Utilities**: Helper functions and mathematical tools for GRB afterglow calculations

Using the C++ API
-----------------

After compiling the library, you can create custom applications that use VegasAfterglow's core functionality:

Including Headers
^^^^^^^^^^^^^^^^^

To use VegasAfterglow directly from C++, include the main header:

.. code-block:: cpp

    #include "afterglow.h"              // Main header for afterglow models

Building Custom Applications
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Compile your C++ application, linking against the VegasAfterglow library:

.. code-block:: bash

    g++ -std=c++20 -I/path/to/VegasAfterglow/include -L/path/to/VegasAfterglow/lib -o my_program my_program.cpp -lvegasafterglow

.. note:: Replace ``/path/to/VegasAfterglow/`` with the actual path to your VegasAfterglow installation.

C++ API Documentation
---------------------

Jet Models
----------

.. doxygenclass:: TophatJet
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: GaussianJet
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: PowerLawJet
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: Ejecta
   :members:
   :undoc-members:
   :allow-dot-graphs:

Ambient Medium
--------------

.. doxygenclass:: ISM
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: Wind
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: Medium
   :members:
   :undoc-members:
   :allow-dot-graphs:

Radiation Processes
-------------------

.. doxygenstruct:: SynPhotons
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenstruct:: SynElectrons
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenstruct:: InverseComptonY
   :members:
   :undoc-members:
   :allow-dot-graphs:

Shock Dynamics
--------------

.. doxygenclass:: Shock
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: SimpleShockEqn
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: ForwardShockEqn
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: FRShockEqn
   :members:
   :undoc-members:
   :allow-dot-graphs:

Physics and Utilities
---------------------

.. doxygenclass:: Observer
   :members:
   :undoc-members:
   :allow-dot-graphs:

.. doxygenclass:: Coord
   :members:
   :undoc-members:
   :allow-dot-graphs:

Performance Considerations
--------------------------

VegasAfterglow's C++ core is designed for exceptional computational performance:

* **Memory Access Patterns**: Carefully optimized to minimize cache misses
* **SIMD Optimizations**: Takes advantage of vectorization where possible
* **Multi-threading**: Core algorithms designed for parallel execution
* **Avoiding Allocations**: Minimal heap allocations in critical computation paths
* **Computational Approximations**: Efficient numerical approximations for complex computations

These optimizations enable the generation of a 30-point single-frequency light curve in approximately 0.6 milliseconds on an Apple M2 chip with a single core, and full MCMC parameter estimation with 10,000 steps in seconds to minutes on standard laptop hardware.

Documenting C++ Code
--------------------

When contributing to the C++ codebase, please follow these documentation guidelines:

Class and Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For member variables, use inline Doxygen comments with the triple-slash syntax:

.. code-block:: cpp

    double energy; ///< Isotropic-equivalent energy in ergs
    double gamma0; ///< Initial bulk Lorentz factor

Template Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For inline functions, use specialized documentation to explain why the function is inline and include important implementation details:

.. code-block:: cpp

    /**
     * @brief Compute the square of a value
     * @inlinefunc Performance-critical function used in inner loops
     *
     * @param x The value to square
     * @return The squared value
     *
     * @inline_details
     * Uses direct multiplication instead of std::pow for better performance.
     * Handles both positive and negative inputs correctly.
     */
    inline double square(double x) {
        return x * x;
    }

C++20 Attribute Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For functions with C++20 attributes, use the specialized tags:

.. code-block:: cpp

    /**
     * @brief Calculate the inverse of a value
     * @nodiscard
     * @constexpr
     *
     * @param value The input value (must not be zero)
     * @return The inverse of the input value (1/value)
     * @throws std::invalid_argument if value is zero
     */
    [[nodiscard]] constexpr double inverse(double value) {
        if (value == 0) throw std::invalid_argument("Cannot take inverse of zero");
        return 1.0 / value;
    }

Special Member Function Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For special member functions, use the dedicated aliases:

.. code-block:: cpp

    /**
     * @defaultctor
     * Initializes with default empty state.
     */
    JetModel();

    /**
     * @copyctor
     * @param other The jet model to copy
     */
    JetModel(const JetModel& other);

    /**
     * @moveassign
     * @param other The jet model to move from
     * @return Reference to this object
     */
    JetModel& operator=(JetModel&& other) noexcept;

Private Member Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Even though private members won't appear in the public API documentation, they should be properly documented in the code for maintainability:

.. code-block:: cpp

    private:
        /**
         * @brief Calculate internal jet dynamics
         *
         * @param time Current simulation time
         * @return Energy distribution at current time
         */
        double calculateDynamics(double time);

        double energy_; ///< Internal energy storage

Example Class
^^^^^^^^^^^^^

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

        /**
         * @brief Get the core angle of the jet
         * @nodiscard
         * @return Core angle in radians
         */
        [[nodiscard]] inline Real getTheta_c() const noexcept;

        /**
         * @brief Get the isotropic equivalent energy
         * @nodiscard
         * @return Energy in ergs
         */
        [[nodiscard]] inline Real getE_iso() const noexcept;

        /**
         * @brief Get the initial Lorentz factor
         * @nodiscard
         * @return Lorentz factor at jet core
         */
        [[nodiscard]] inline Real getGamma0() const noexcept;
    };

    // Implementation of inline methods would be in the .cpp file or in a separate
    // inline header file, and should not appear in the API documentation.

For more details on Doxygen commands, see the :doc:`contributing` page.
