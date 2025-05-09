Contributing
============

.. contents:: Table of Contents
   :local:
   :depth: 2

Getting Started
---------------

Thank you for considering contributing to VegasAfterglow! This document provides guidelines for contributing to the project.

Installation for Development
----------------------------

1. Fork the repository on GitHub
2. Clone your fork locally:

.. code-block:: bash

    git clone https://github.com/your-username/VegasAfterglow.git
    cd VegasAfterglow

3. Install the package in development mode:

.. code-block:: bash

    pip install -e .

4. Install development dependencies:

.. code-block:: bash

    pip install pytest black flake8 sphinx sphinx_rtd_theme breathe

Coding Standards
----------------

C++ Code
^^^^^^^^

* Follow modern C++20 practices
* Use appropriate type safety features
* Write clear, well-commented code
* Avoid raw pointers, prefer smart pointers
* Keep functions small and focused on a single task
* Use meaningful variable and function names
* Leverage C++20 features when appropriate (concepts, ranges, etc.)

Python Code
^^^^^^^^^^^

* Follow PEP 8 style guidelines
* Use type hints for function signatures
* Keep functions small and focused
* Document all functions, classes, and modules with docstrings
* Support Python 3.8 and higher

Documentation Guidelines
------------------------

Documentation for VegasAfterglow consists of both C++ API documentation (using Doxygen) and Python API documentation (using Sphinx with Napoleon). Proper documentation is essential for maintaining code quality and improving user experience.

C++ Documentation
^^^^^^^^^^^^^^^^^

All C++ code should be documented using Doxygen-style comments:

* **Classes and Functions**: Use ``/**`` style comments with appropriate tags:

  .. code-block:: cpp

      /**
       * @brief Brief description of the function or class
       * 
       * Detailed description that can span multiple lines and provide
       * more in-depth information about the functionality.
       *
       * @param paramName Description of the parameter
       * @return Description of the return value
       * @throws Description of exceptions that might be thrown
       * @see RelatedClass, anotherFunction()
       */
      ReturnType functionName(ParamType paramName);

* **Member Variables**: Use ``///<`` for inline comments:

  .. code-block:: cpp

      double mass; ///< Mass of the object in solar masses

* **Formulas**: Use LaTeX notation for mathematical expressions:

  .. code-block:: cpp

      /**
       * @brief Calculates the synchrotron spectral power
       * 
       * The formula used is:
       * \f[
       *    P(\nu) = \frac{4\pi}{3} r_e m_e c^2 \gamma^2 B^2 F(\nu/\nu_c)
       * \f]
       * where \f$r_e\f$ is the classical electron radius and \f$F\f$ is the synchrotron function.
       */

Python Documentation
^^^^^^^^^^^^^^^^^^^^

For Python code, use NumPy-style docstrings:

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

Building Documentation
----------------------

To build the documentation locally:

1. Ensure you have installed all documentation dependencies:

.. code-block:: bash

    pip install sphinx sphinx_rtd_theme breathe
    # Install Doxygen and Graphviz using your package manager
    # For example, on Ubuntu:
    # sudo apt-get install doxygen graphviz
    # On macOS:
    # brew install doxygen graphviz

2. Build the documentation:

.. code-block:: bash

    cd docs
    make all

3. View the generated documentation:

.. code-block:: bash

    # Open docs/build/html/index.html in your browser

Testing
-------

Before submitting a pull request, please run the test suite:

.. code-block:: bash

    python -m pytest tests/

For C++ tests, you can use:

.. code-block:: bash

    make tests

Pull Request Process
--------------------

1. Update the README.md and documentation with details of changes, if applicable
2. Update the test suite to include tests for your changes
3. Make sure all tests pass
4. Submit a pull request to the main repository
5. Your pull request will be reviewed by the maintainers

Issues and Feedback
-------------------

If you encounter any issues, have questions about the code, or want to request new features:

1. **GitHub Issues** - The most straightforward and fastest way to get help:
   - Open an issue at https://github.com/YihanWangAstro/VegasAfterglow/issues
   - You can report bugs, suggest features, or ask questions
   - This allows other users to see the problem/solution as well

Code of Conduct
---------------

Please be respectful and considerate of others when contributing to this project. Harassment and abusive behavior will not be tolerated.

License
-------

By contributing to this project, you agree that your contributions will be licensed under the project's BSD-3-Clause License. 