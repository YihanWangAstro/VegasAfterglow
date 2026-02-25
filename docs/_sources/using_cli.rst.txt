Command-Line Interface (``vegasgen``)
======================================

VegasAfterglow includes a command-line tool, ``vegasgen``, that generates multi-band afterglow light curves directly from the terminal. It provides quick access to the full physics engine without writing any Python code.

When to use the CLI vs the Python API:

- **CLI** — Quick exploration, parameter sweeps from shell scripts, generating CSV/JSON data or publication-quality plots with a single command.
- **Python API** — Custom analysis pipelines, MCMC fitting, accessing internal quantities (shock dynamics, spectra), or integrating with other Python libraries.


Installation
------------

``vegasgen`` is installed automatically with the package:

.. code-block:: bash

    pip install VegasAfterglow

Verify the installation:

.. code-block:: bash

    vegasgen --help


Basic Usage
-----------

With no arguments, ``vegasgen`` computes a light curve using sensible defaults (top-hat jet, ISM medium, on-axis observer) and prints CSV data to stdout:

.. code-block:: bash

    vegasgen

Pipe to a file or redirect output:

.. code-block:: bash

    vegasgen > lightcurve.csv
    vegasgen -o lightcurve.csv

Generate a plot instead of data:

.. code-block:: bash

    vegasgen --plot                    # interactive window
    vegasgen --plot -o lightcurve.png  # save to file


Parameter Reference
-------------------

All parameters have defaults, so you only need to specify what you want to change.

Jet Parameters
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--jet``
     - choice
     - ``tophat``
     - Jet structure: ``tophat``, ``gaussian``, or ``powerlaw``
   * - ``--theta_c``
     - float
     - ``0.1``
     - Half-opening angle [rad]
   * - ``--E_iso``
     - float
     - ``1e52``
     - Isotropic-equivalent energy [erg]
   * - ``--Gamma0``
     - float
     - ``300``
     - Initial Lorentz factor
   * - ``--k_e``
     - float
     - ``2``
     - Energy power-law index (``powerlaw`` jet only)
   * - ``--k_g``
     - float
     - ``2``
     - Lorentz factor power-law index (``powerlaw`` jet only)
   * - ``--spreading``
     - flag
     - off
     - Enable lateral spreading

Medium Parameters
^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--medium``
     - choice
     - ``ism``
     - Circumburst medium: ``ism`` or ``wind``
   * - ``--n_ism``
     - float
     - ``1.0``
     - ISM number density [cm\ :sup:`-3`]
   * - ``--A_star``
     - float
     - ``0.1``
     - Wind parameter A\ :sub:`*`

Observer Parameters
^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--z``
     - float
     - ``0.01``
     - Redshift
   * - ``--lumi_dist``
     - float
     - auto
     - Luminosity distance [cm]. If omitted, estimated from redshift using a Hubble-law approximation
   * - ``--theta_obs``
     - float
     - ``0``
     - Observer viewing angle [rad]

Radiation Parameters
^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--eps_e``
     - float
     - ``0.1``
     - Fraction of shock energy in electrons
   * - ``--eps_B``
     - float
     - ``1e-3``
     - Fraction of shock energy in magnetic field
   * - ``--p``
     - float
     - ``2.3``
     - Electron spectral index
   * - ``--xi_e``
     - float
     - ``1``
     - Fraction of electrons accelerated
   * - ``--ssc``
     - flag
     - off
     - Enable synchrotron self-Compton
   * - ``--kn``
     - flag
     - off
     - Enable Klein-Nishina corrections (requires ``--ssc``)
   * - ``--cmb_cooling``
     - flag
     - off
     - Enable CMB inverse Compton cooling

Frequency Specification
^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--nu``
     - list
     - ``1e9 5e14 1e18``
     - One or more frequencies. Accepts numeric Hz values or filter names (see :ref:`filter-names`)

Time Grid
^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--t_min``
     - float
     - ``100``
     - Start time [s]
   * - ``--t_max``
     - float
     - ``1e8``
     - End time [s]
   * - ``--num_t``
     - int
     - ``200``
     - Number of time points (logarithmically spaced)

Resolution
^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``--res PHI THETA T``
     - 3 floats
     - ``0.15 0.5 10``
     - Grid resolution: phi points per degree, theta points per degree, time points per decade

Output Options
^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1
   :widths: 20 15 15 50

   * - Argument
     - Type
     - Default
     - Description
   * - ``-o``, ``--output``
     - path
     - stdout
     - Output file path
   * - ``--format``
     - choice
     - ``csv``
     - Data format: ``csv`` or ``json``
   * - ``--flux_unit``
     - choice
     - ``mJy``
     - Flux unit: ``mJy``, ``Jy``, ``uJy``, or ``cgs``
   * - ``--time_unit``
     - choice
     - ``s``
     - Time unit: ``s``, ``day``, ``hr``, or ``min``
   * - ``--plot``
     - flag
     - off
     - Generate a plot instead of data output
   * - ``--font``
     - string
     - Times New Roman
     - Plot font family (e.g. ``Helvetica``, ``Palatino``). Sans-serif fonts are auto-detected


.. _filter-names:

Frequency Filters
-----------------

The ``--nu`` argument accepts both numeric Hz values and standard photometric filter names. You can freely mix them:

.. code-block:: bash

    vegasgen --nu 1e9 R F606W

Available filter names:

**Vega system (Johnson-Cousins, 2MASS, Swift UVOT)**

.. list-table::
   :header-rows: 1
   :widths: 30 35 35

   * - System
     - Filters
     - Reference
   * - Johnson-Cousins
     - U, B, V, R, I
     - Bessell & Murphy (2012)
   * - 2MASS
     - J, H, Ks
     - Cohen et al. (2003)
   * - Swift UVOT
     - v, b, u, uvw1, uvm2, uvw2
     - Breeveld et al. (2011)

**ST system (HST WFC3)**

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Instrument
     - Filters
   * - WFC3/UVIS
     - F225W, F275W, F336W, F438W, F475W, F555W, F606W, F625W, F775W, F814W, F850LP
   * - WFC3/IR
     - F105W, F110W, F125W, F140W, F160W

Each filter name is converted to an effective frequency using its tabulated wavelength.


Output Formats
--------------

CSV (default)
^^^^^^^^^^^^^

.. code-block:: bash

    vegasgen -o lightcurve.csv

Produces:

.. code-block:: text

    # VegasAfterglow light curve
    # jet=tophat theta_c=0.1 E_iso=1.0e+52 Gamma0=300.0 medium=ism z=0.01 theta_obs=0
    # eps_e=0.1 eps_B=0.001 p=2.3
    # t_unit=s flux_unit=mJy
    t(s),F_total(radio),F_total(optical),F_total(xray)
    1.000000e+02,1.234567e-02,4.567890e+00,7.890123e-01
    ...

JSON
^^^^

.. code-block:: bash

    vegasgen --format json -o lightcurve.json

Produces a structured JSON object:

.. code-block:: json

    {
      "parameters": {
        "jet": "tophat",
        "theta_c": 0.1,
        "E_iso": 1e52,
        "Gamma0": 300.0,
        "medium": "ism",
        "z": 0.01,
        "theta_obs": 0,
        "eps_e": 0.1,
        "eps_B": 0.001,
        "p": 2.3
      },
      "units": {"time": "s", "flux": "mJy"},
      "frequencies_Hz": [1e9, 5e14, 1e18],
      "times": [100.0, "..."],
      "flux": {
        "radio": ["..."],
        "optical": ["..."],
        "xray": ["..."]
      }
    }


Plotting
--------

The ``--plot`` flag generates publication-quality figures using matplotlib.

.. code-block:: bash

    # Interactive display
    vegasgen --plot

    # Save to file (PNG at 300 DPI, or PDF/SVG as vector)
    vegasgen --plot -o lightcurve.png
    vegasgen --plot -o lightcurve.pdf

Plot features:

- **VegasAfterglow signature palette** — a vibrant, high-contrast color scheme designed for multi-band light curves
- **Times New Roman** font by default, with LaTeX rendering. Change with ``--font``:

  .. code-block:: bash

      vegasgen --plot --font "Helvetica"
      vegasgen --plot --font "Palatino"

  Sans-serif fonts (Helvetica, Arial, etc.) are auto-detected.

- **Log-log** axes with inward-facing ticks
- **Context-aware frequency labels**: GHz for radio, nm for optical/IR, keV for X-ray, filter band names when applicable
- **Parameter summary** displayed above the plot
- **Single-column journal size** (3.5 in wide) suitable for direct inclusion in publications

Supported output formats: ``.png``, ``.pdf``, ``.jpg``, ``.svg``. If ``-o`` is not specified or does not end in a recognized image extension, an interactive matplotlib window is shown.

.. note::

   Plotting requires matplotlib. Install with ``pip install matplotlib`` or ``pip install VegasAfterglow[test]``.


Examples
--------

Default on-axis light curve
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    vegasgen

Computes a top-hat jet in ISM with default parameters at radio (1 GHz), optical (5×10\ :sup:`14` Hz), and X-ray (10\ :sup:`18` Hz).

Off-axis observer
^^^^^^^^^^^^^^^^^

.. code-block:: bash

    vegasgen --theta_obs 0.4 --plot

Sets the viewing angle to 0.4 rad (about 23°), well outside the default jet core of 0.1 rad.

Different jet structures
^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Gaussian jet
    vegasgen --jet gaussian --theta_c 0.05 --E_iso 1e53

    # Power-law jet
    vegasgen --jet powerlaw --k_e 4 --k_g 2

Wind medium
^^^^^^^^^^^

.. code-block:: bash

    vegasgen --medium wind --A_star 0.01 --plot

Synchrotron self-Compton with Klein-Nishina
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    vegasgen --ssc --kn --plot

Custom frequency bands
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Named filters
    vegasgen --nu R J Ks --plot

    # Mixed numeric and filter names
    vegasgen --nu 5e9 R 1e18 -o multiband.csv

    # HST filters
    vegasgen --nu F606W F160W --flux_unit uJy --plot

Custom time range and units
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # Early afterglow in seconds
    vegasgen --t_min 1 --t_max 1e4 --num_t 100

    # Late afterglow in days
    vegasgen --t_min 86400 --t_max 1e8 --time_unit day -o late.csv

Higher resolution
^^^^^^^^^^^^^^^^^

.. code-block:: bash

    vegasgen --res 0.1 0.3 20

Sets finer angular (phi=0.1, theta=0.3 points per degree) and temporal (20 points per decade) resolution.

Output in different formats and units
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: bash

    # JSON format in Jansky
    vegasgen --format json --flux_unit Jy -o lightcurve.json

    # CSV in microjansky and days
    vegasgen --flux_unit uJy --time_unit day -o lightcurve.csv
