Troubleshooting
===============

This page addresses common issues and questions that users encounter when working with VegasAfterglow.

Common Issues
-------------

Model Setup and Calculation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q: My light curve has strange spikes or discontinuities**

A: This is usually caused by one of the following:

- **Time array ordering**: Ensure your time array is in ascending order when using ``specific_flux()``.
- **Resolution too low**: Increase the resolution parameters in ``Model(resolutions=(phi_ppd, theta_ppd, t_ppd))``
- **User-defined profiles**: For custom jet/medium profiles, ensure they are smooth and well-behaved

**Q: My calculation is extremely slow**

A: Performance can be improved by:

- **Reducing resolution**: Use ``resolutions=(0.1, 1, 10)`` for speed, ``(0.3, 5, 20)`` for accuracy
- **Limiting frequency/time ranges**: Calculate only the bands and times you need
- **Using built-in profiles**: Built-in jet structures are faster than user-defined Python functions
- **For MCMC**: Consider using fewer parameters or coarser resolution

**Q: I get numerical errors or NaN values**

A: Check for:

- **Extreme parameter values**: Ensure parameters are within physically reasonable ranges
- **Zero or negative values**: Some parameters (like energies, densities) must be positive
- **Resolution mismatch**: Very fine resolution with short time arrays can cause issues

**Q: My model doesn't match basic expectations (wrong slope, normalization)**

A: Verify:

- **Units**: All inputs should be in CGS units (see parameter reference table)
- **Observer frame vs source frame**: Times and frequencies are in observer frame
- **Jet opening angle**: ``theta_c`` is in radians, not degrees

MCMC Fitting Issues
^^^^^^^^^^^^^^^^^^^

**Q: MCMC chains don't converge or get stuck**

A: Try the following:

- **Check parameter ranges**: Ensure prior ranges include the true values
- **Use better initial guesses**: Set ``ParamDef`` initial values closer to expected results
- **Increase burn-in**: Use ``burn_frac=0.5`` or higher for difficult problems
- **More walkers**: Use ``n_walkers > 2 * n_parameters`` (emcee recommendation)
- **Check data quality**: Ensure observational uncertainties are realistic

**Q: MCMC is too slow for practical use**

A: Optimization strategies:

- **Reduce resolution**: Use ``resolution=(0.3, 1, 10)`` for initial exploration
- **Fewer parameters**: Fix some parameters with ``Scale.FIXED``
- **Coarser time/frequency grids**: Use fewer data points for initial fits
- **Parallel processing**: Ensure you're using multiple cores

**Q: Parameter constraints seem unrealistic**

A: Check:

- **Parameter scaling**: Use ``Scale.LOG`` for parameters spanning orders of magnitude
- **Prior ranges**: Ensure they're physically motivated and not too restrictive
- **Model degeneracies**: Some parameters may be strongly correlated
- **Data coverage**: Limited frequency/time coverage can lead to poor constraints

Data and File Issues
^^^^^^^^^^^^^^^^^^^^

**Q: I can't load my observational data**

A: Common data loading issues:

- **File format**: Ensure CSV files have the expected column names (``t``, ``Fv_obs``, ``Fv_err``, ``nu``)
- **Units**: Data should be in CGS units (times in seconds, frequencies in Hz, fluxes in erg/cm²/s/Hz)
- **Missing values**: Remove or interpolate NaN/infinite values
- **File paths**: Use absolute paths or ensure files are in the correct directory

**Q: Error messages about missing dependencies**

A: Install required packages:

.. code-block:: bash

    pip install numpy scipy matplotlib pandas corner emcee

For specific features:

.. code-block:: bash

    pip install jupyter  # For notebook examples
    pip install h5py     # For saving large datasets

Installation and Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Q: Installation fails on my system**

A: Platform-specific solutions:

- **macOS Apple Silicon**: Try ``pip install --no-deps VegasAfterglow`` then install dependencies separately
- **Windows**: Ensure Visual Studio Build Tools are installed
- **Linux**: May need development packages (``python3-dev``, ``build-essential``)
- **Conda environments**: Use ``pip`` within conda, not ``conda install``

**Q: ImportError when importing VegasAfterglow**

A: Check:

- **Python version**: VegasAfterglow requires Python 3.8+
- **Virtual environment**: Ensure you're in the correct environment
- **Installation location**: Try ``pip show VegasAfterglow`` to verify installation
- **Conflicting packages**: Try installing in a clean environment

Performance Guidelines
----------------------

Resolution Parameters
^^^^^^^^^^^^^^^^^^^^^

The ``resolutions`` parameter in ``Model()`` controls computational accuracy vs speed:

.. list-table:: Resolution Guidelines
   :header-rows: 1
   :widths: 20 25 25 30

   * - Use Case
     - Resolution
     - Speed
     - Accuracy
   * - Initial exploration
     - ``(0.2, 1, 5)``
     - Very Fast
     - Low
   * - Standard calculations
     - ``(0.3, 1, 10)``
     - Fast
     - Good
   * - MCMC fitting
     - ``(0.3, 2, 15)``
     - Moderate
     - High
   * - Publication quality
     - ``(0.3, 5, 20)``
     - Slow
     - Very High

Where ``resolutions=(phi_ppd, theta_ppd, t_ppd)``:

- ``phi_ppd``: Points per degree in azimuthal direction
- ``theta_ppd``: Points per degree in polar direction. The code sets a minimum of 56 points across the jet profile.
- ``t_ppd``: Points per decade in time direction. The code sets a minimum of 24 time points.

Memory Usage
^^^^^^^^^^^^

For large parameter studies or high-resolution calculations:

- **Limit output arrays**: Calculate only needed times/frequencies
- **Use generators**: Process results in chunks rather than storing everything
- **Clear variables**: Use ``del`` to free memory between calculations
- **Monitor usage**: Use ``htop`` or Task Manager to monitor memory consumption

Getting Help
------------

If you encounter issues not covered here:

1. **Check the examples**: The :doc:`examples` page covers many common use cases
2. **Search existing issues**: Visit our `GitHub Issues <https://github.com/YihanWangAstro/VegasAfterglow/issues>`_
3. **Create a new issue**: Include:

   - VegasAfterglow version: ``import VegasAfterglow; print(VegasAfterglow.__version__)``
   - Python version and platform
   - Minimal code example that reproduces the problem
   - Full error traceback

4. **Discussion forum**: For general questions about GRB physics or methodology

Best Practices
--------------

Model Development Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^

1. **Start simple**: Begin with built-in jet types and standard parameters
2. **Validate physics**: Check that results match analytical expectations for simple cases
3. **Parameter exploration**: Use direct model calculations before MCMC
4. **Incremental complexity**: Add features (reverse shock, IC, etc.) one at a time
5. **Resolution testing**: Verify results are converged by increasing resolution
