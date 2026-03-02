Basic Setup
===========

Fitter Constructor
------------------

The ``Fitter`` constructor accepts keyword arguments for source properties, model selection, physics options, and numerical parameters:

.. code-block:: python

    fitter = Fitter(
        # Source properties
        z=1.58,                   # Redshift
        lumi_dist=3.364e28,       # Luminosity distance [cm]

        # Model selection (see sections below for all options)
        jet="powerlaw",           # Jet structure type
        medium="wind",            # Ambient medium type

        # Physics options
        rvs_shock=True,           # Include reverse shock
        fwd_ssc=True,             # Forward shock inverse Compton
        rvs_ssc=False,            # Reverse shock inverse Compton
        kn=True,                  # Klein-Nishina corrections
        magnetar=True,            # Magnetar energy injection

        # Numerical parameters
        rtol=1e-5,                # Numerical tolerance
        resolution=(0.1, 0.5, 5),   # Grid resolution (phi, theta, t)
    )

Setting up Data and the Fitter
-------------------------------

The ``Fitter`` class handles both model configuration and observational data. Add light curves, spectra, and broadband flux directly to the fitter:

.. code-block:: python

    import numpy as np
    import pandas as pd
    from VegasAfterglow import Fitter, ParamDef, Scale

    # Create the fitter with model configuration
    fitter = Fitter(
        z=1.58,
        lumi_dist=3.364e28,
        jet="tophat",
        medium="ism",
    )

    # Method 1: Add data directly
    t_data = [1e3, 2e3, 5e3, 1e4, 2e4]  # Time in seconds
    flux_data = [1e-26, 8e-27, 5e-27, 3e-27, 2e-27]  # erg/cm²/s/Hz
    flux_err = [1e-28, 8e-28, 5e-28, 3e-28, 2e-28]   # Error bars

    # Add light curve at R-band frequency
    fitter.add_flux_density(nu=4.84e14, t=t_data,
                            f_nu=flux_data, err=flux_err)  # All quantities in CGS units

    # Optional: Add weights for systematic uncertainties, normalization handled internally
    weights = np.ones_like(t_data)  # Equal weights
    fitter.add_flux_density(nu=2.4e17, t=t_data,
                            f_nu=flux_data, err=flux_err,
                            weights=weights)  # All quantities in CGS units

    # Method 2: Add frequency-integrated light curve (broadband flux)
    # For instruments with wide energy bands (e.g., Swift XRT, Fermi LAT)
    from VegasAfterglow.units import band, keV

    fitter.add_flux(band=band("XRT"), t=t_data,       # named band
                    flux=flux_data, err=flux_err,
                    weights=weights)  # flux in erg/cm²/s

    # Or specify a custom band as (nu_min, nu_max) tuple
    fitter.add_flux(band=(0.3*keV, 10*keV), t=t_data,
                    flux=flux_data, err=flux_err)

    # Method 3: Load from CSV files
    bands = [2.4e17, 4.84e14, 1.4e14]  # X-ray, optical, near-IR
    lc_files = ["data/xray.csv", "data/optical.csv", "data/nir.csv"]

    for nu, fname in zip(bands, lc_files):
        df = pd.read_csv(fname)
        fitter.add_flux_density(nu=nu, t=df["t"],
                                f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

    # Add spectra at specific times
    spec_times = [1000, 10000]  # seconds
    spec_files = ["data/spec_1000s.csv", "data/spec_10000s.csv"]

    for t, fname in zip(spec_times, spec_files):
        df = pd.read_csv(fname)
        fitter.add_spectrum(t=t, nu=df["nu"],
                            f_nu=df["Fv_obs"], err=df["Fv_err"])  # All quantities in CGS units

Unit Conversion
----------------

All API inputs use CGS base units (seconds, Hz, erg/cm²/s/Hz, cm, radians). The ``VegasAfterglow.units`` module provides multiplicative constants for convenient conversion from common astronomical units:

.. code-block:: python

    from VegasAfterglow.units import day, GHz, keV, mJy, Mpc, deg

    # Fitter setup with natural units
    fitter = Fitter(z=1.58, lumi_dist=100*Mpc, jet="tophat", medium="ism")

    # Load data in astronomer-friendly units
    df = pd.read_csv("data/radio.csv")
    fitter.add_flux_density(
        nu=5*GHz,
        t=df["t_days"] * day,
        f_nu=df["flux_mJy"] * mJy,
        err=df["err_mJy"] * mJy,
    )

    # X-ray band specified in keV (converted to Hz via E=hν)
    fitter.add_flux_density(
        nu=1*keV,
        t=df_xray["t_days"] * day,
        f_nu=df_xray["flux"] * mJy,
        err=df_xray["err"] * mJy,
    )

    # Convert model output back to convenient units
    result_flux = model.flux_density_grid(times, nus)
    flux_in_mJy = result_flux.total / mJy

Available constants: ``sec``, ``ms``, ``minute``, ``hr``, ``day``, ``yr`` (time); ``Hz``, ``kHz``, ``MHz``, ``GHz`` (frequency); ``eV``, ``keV``, ``MeV``, ``GeV`` (photon energy → frequency); ``Jy``, ``mJy``, ``uJy`` (flux density); ``cm``, ``m``, ``km``, ``pc``, ``kpc``, ``Mpc``, ``Gpc`` (distance); ``rad``, ``deg``, ``arcmin``, ``arcsec`` (angle).

Magnitude Conversion
^^^^^^^^^^^^^^^^^^^^^

If your data is in magnitudes rather than flux density, use the built-in magnitude converters:

.. code-block:: python

    from VegasAfterglow.units import ABmag_to_cgs, Vegamag_to_cgs, filter

    # AB magnitudes (e.g., SDSS photometry)
    f_nu = ABmag_to_cgs(df["mag_AB"])
    f_nu_err = ABmag_to_cgs(df["mag_AB"] - df["mag_err"]) - f_nu  # asymmetric!

    # Vega magnitudes with filter name (e.g., Johnson R band)
    f_nu = Vegamag_to_cgs(df["mag_R"], "R")

    # Get the effective frequency for a filter (for the fitter)
    fitter.add_flux_density(nu=filter("R"), t=t_data, f_nu=f_nu, err=f_nu_err)

Supported magnitude systems and filters:

- **Vega magnitudes** via ``Vegamag_to_cgs(mag, filter_name)``: Johnson-Cousins (``U``, ``B``, ``V``, ``R``, ``I``), 2MASS (``J``, ``H``, ``Ks``), Swift UVOT (``v``, ``b``, ``u``, ``uvw1``, ``uvm2``, ``uvw2``), SDSS (``g``, ``r``, ``i``, ``z``)
- **AB magnitudes** via ``ABmag_to_cgs(mag)``: universal zero-point, works for any filter without specifying a name
- **ST magnitudes** via ``STmag_to_cgs(mag, filter_name)``: HST WFC3/UVIS (``F225W`` – ``F850LP``), WFC3/IR (``F105W`` – ``F160W``)

The ``filter()`` function returns the effective frequency [Hz] for any supported filter name, including survey filters (``VT_B``, ``VT_R``, ``w``) that only support AB magnitude conversion.

Data Selection and Optimization
--------------------------------

**Smart Data Subsampling with logscale_screen**

For large datasets or densely sampled observations, using all available data points can lead to computational inefficiency and biased parameter estimation. The ``logscale_screen`` utility provides intelligent data subsampling that maintains the essential information content while reducing computational overhead.

.. code-block:: python

    from VegasAfterglow import logscale_screen

    # Example: Large dense dataset
    t_dense = np.logspace(2, 7, 1000)  # 1000 time points
    flux_dense = np.random.lognormal(-60, 0.5, 1000)  # Dense flux measurements
    flux_err_dense = 0.1 * flux_dense

    # Subsample using logarithmic screening
    # This selects ~5*5=25 representative points across 5 decades in time
    indices = logscale_screen(t_dense, data_density=5)

    # Add only the selected subset
    fitter.add_flux_density(nu=5e14,
                            t=t_dense[indices],
                            f_nu=flux_dense[indices],
                            err=flux_err_dense[indices])

**Why logscale_screen is Important:**

1. **Prevents Oversampling Bias**: Dense data clusters can dominate the chi-squared calculation, causing the MCMC to over-fit specific frequency bands or time periods.

2. **Computational Efficiency**: Reduces the number of model evaluations needed during MCMC sampling, significantly improving performance.

3. **Preserves Information**: Unlike uniform thinning, logarithmic sampling maintains representation across all temporal/spectral scales.

4. **Balanced Multi-band Fitting**: Ensures each frequency band contributes proportionally to the parameter constraints.

**Data Selection Guidelines:**

- **Target 10-30 points per frequency band** for balanced constraints
- **Avoid >100 points in any single band** unless scientifically justified
- **Maintain temporal coverage** across all evolutionary phases
- **Weight systematic uncertainties** appropriately using the weights parameter

.. warning::
    **Common Data Selection Pitfalls:**

    - **Optical-heavy datasets**: Dense optical coverage can bias parameters toward optical-dominant solutions
    - **Late-time clustering**: Too many late-time points can over-constrain decay slopes at the expense of early physics
    - **Single-epoch spectra**: Broadband spectra at one time can dominate multi-epoch light curves in chi-squared space

    **Solution**: Use ``logscale_screen()`` for manual temporal reduction of over-sampled bands.

**Handling Measurement Uncertainties**

Error bars strongly influence MCMC fitting. Improperly characterized uncertainties can lead to biased parameter estimates or misleading confidence intervals.

**The Problem with Small Error Bars**

Data points with very small uncertainties dominate the chi-squared calculation:

.. math::

    \chi^2 = \sum_i \frac{(F_{\rm obs,i} - F_{\rm model,i})^2}{\sigma_i^2}

A single point with sigma = 0.01 mJy contributes 100x more to chi-squared than a point with sigma = 0.1 mJy. This causes the MCMC to:

- Heavily weight high-precision data at the expense of other observations
- Potentially fit noise or systematic effects in the "precise" data
- Produce artificially tight parameter constraints that don't reflect true uncertainties

**Common Scenarios:**

1. **Heterogeneous data quality**: X-ray data from different telescopes may have vastly different calibration uncertainties
2. **Underestimated errors**: Published error bars sometimes only include statistical uncertainty, missing systematics
3. **Calibration offsets**: Different instruments may have 5-20% flux calibration offsets not captured in quoted errors

**Solutions:**

1. **Add systematic error floors**

   For data where quoted errors seem unrealistically small:

   .. code-block:: python

       # Add 10% systematic floor to error bars
       systematic_floor = 0.1  # 10% of flux
       err_with_floor = np.sqrt(err**2 + (systematic_floor * flux)**2)

       fitter.add_flux_density(nu=nu, t=t, f_nu=flux, err=err_with_floor)

2. **Use weights to balance bands**

   Down-weight bands with suspiciously small errors:

   .. code-block:: python

       # Reduce contribution of high-precision optical data
       optical_weight = 0.5  # Effectively doubles the error contribution
       fitter.add_flux_density(nu=5e14, t=t_optical, f_nu=flux_optical,
                              err=err_optical, weights=optical_weight * np.ones_like(t_optical))

3. **Inspect residuals by band**

   After fitting, check if residuals are consistent across bands. If one band shows systematic deviations while others don't, the errors for that band may be mischaracterized.

.. tip::
    **Rule of thumb**: If your best-fit chi-squared/dof >> 1, either the model is inadequate OR the error bars are underestimated. If chi-squared/dof << 1, the errors may be overestimated.
