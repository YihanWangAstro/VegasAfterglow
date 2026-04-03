MCMC Best Practices
===================

Start Simple, Add Complexity Gradually
--------------------------------------

The most common mistake in GRB afterglow fitting is starting with an overly complex model. **Always begin with the simplest physically motivated model** and only add complexity when the data clearly demands it.

**Recommended Progression:**

1. **Start with TopHat + ISM + Forward Shock Only**

   .. code-block:: python

       fitter = Fitter(z=1.58, lumi_dist=3.364e28, jet="tophat", medium="ism")
       # All other physics options default to False

   This gives you ~7-8 free parameters. Run MCMC and examine the residuals.

2. **Check if residuals suggest additional physics**

   Examine the fit quality and residual structure. Systematic deviations at specific times or frequencies may indicate missing physics.

3. **Add ONE component at a time**

   Never jump from a simple model to enabling everything. Each addition should be justified by improved fit statistics (e.g., Bayesian evidence comparison).

Why Complex Models Are Problematic
-----------------------------------

**1. Parameter Degeneracies**

More parameters create more degeneracies. For example:

- ``E_iso`` and ``n_ism`` are degenerate in flux normalization
- ``eps_e`` and ``eps_B`` trade off in spectral shape
- ``theta_c`` and ``theta_v`` correlate for off-axis observers
- Reverse shock parameters can mimic forward shock with different microphysics

With a complex model, the MCMC may find multiple solutions that fit equally well but have very different physical interpretations.

**2. Computational Cost**

Each additional physics module increases computation time. A model with all physics enabled can be significantly slower than a basic forward-shock-only model.

**3. Overfitting Risk**

With enough free parameters, you can fit noise. A model that fits your data perfectly but has 15+ parameters may not be physically meaningful. Use Bayesian evidence (from dynesty) to compare models:

.. code-block:: python

    # Compare two models using log evidence
    result_simple = fitter_simple.fit(params_simple, sampler="dynesty", ...)
    result_complex = fitter_complex.fit(params_complex, sampler="dynesty", ...)

    # Bayes factor
    log_BF = result_complex.bilby_result.log_evidence - result_simple.bilby_result.log_evidence

    # Interpretation:
    # log_BF < 1: No preference (stick with simple model)
    # 1 < log_BF < 3: Weak preference for complex
    # 3 < log_BF < 5: Moderate preference
    # log_BF > 5: Strong preference for complex model

**4. Non-Physical Solutions**

Complex models can converge to non-physical parameter combinations. Always check:

- Is ``eps_e + eps_B < 1``? (energy conservation)
- Is ``p > 2``? (required for finite electron energy)
- Are microphysics parameters consistent between forward/reverse shocks?

When to Use Complex Models
---------------------------

Complex models are justified when:

- **Clear observational signatures** that simple models cannot explain after careful analysis
- **Comparison with similar GRBs** where the additional physics was robustly established

.. important::
    **Golden Rule**: If you cannot clearly explain WHY each physics component is needed based on your data, you probably don't need it.
