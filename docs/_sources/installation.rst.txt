Installation
============

Prerequisites
-------------

VegasAfterglow requires the following to build:

.. note::
   If you install via pip (recommended), you generally do not need to install these C++ tools manually. This section is primarily for users building the C++ library directly or installing the Python package from the source code.

* **C++20 compatible compiler**:

  * **Linux**: GCC 10+ or Clang 13+
  * **macOS**: Apple Clang 13+ (with Xcode 13+) or GCC 10+ (via Homebrew)
  * **Windows**: MSVC 19.29+ (Visual Studio 2019 16.10+) or MinGW-w64 with GCC 10+

* **CMake 3.15+**:
  
  * Required for building the python package from source

* **Build tools**:

  * Make (GNU Make 4.0+ recommended) [if you want to compile & run the C++ code]

Python Installation
-------------------

From PyPI (Recommended)
^^^^^^^^^^^^^^^^^^^^^^^

The easiest way to install VegasAfterglow is via pip:

.. code-block:: bash

    pip install VegasAfterglow

VegasAfterglow requires Python 3.8 or higher.

From Source
^^^^^^^^^^^

To install from source, first clone the repository:

.. code-block:: bash

    git clone https://github.com/YihanWangAstro/VegasAfterglow.git
    cd VegasAfterglow

Then, build and install:

.. code-block:: bash

    pip install .

GitHub Repository
----------------

The source code is available on GitHub:

`https://github.com/YihanWangAstro/VegasAfterglow <https://github.com/YihanWangAstro/VegasAfterglow>`_

C++ Installation
----------------

For advanced users who want to compile and use the C++ library directly:

1. Clone the repository (if you haven't already):

   .. code-block:: bash

       git clone https://github.com/YihanWangAstro/VegasAfterglow.git
       cd VegasAfterglow

2. Compile the static library:

   .. code-block:: bash

       make lib

   This allows you to write your own C++ problem generator and use the provided VegasAfterglow interfaces.

3. (Optional) Compile and run tests:

   .. code-block:: bash

       make tests

Requirements
------------

* Python 3.8 or higher
* C++20 compatible compiler (for building from source)
* NumPy, SciPy, and other dependencies (automatically installed when using pip) 