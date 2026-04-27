============
Installation
============

System requirements
-------------------

* **Python** ≥ 3.9
* **Operating system**: Linux, macOS, or Windows
* **Memory**: ≥ 2 GB RAM recommended for populations of ≥ 50 chains

Dependencies installed automatically
--------------------------------------

.. list-table::
   :header-rows: 1
   :widths: 20 50 15

   * - Package
     - Role
     - Min version
   * - ``networkx``
     - Molecular graph construction, traversal, spring-layout
     - 2.8
   * - ``numpy``
     - Numerical arrays, random sampling, distance metrics
     - 1.22
   * - ``pandas``
     - Population statistics DataFrames and CSV export
     - 1.4
   * - ``rdkit``
     - Molecule construction, 3-D embedding (ETKDGv3), MMFF/UFF
     - 2022.09
   * - ``pysmiles``
     - SMILES generation from NetworkX molecular graphs
     - 1.0
   * - ``jsonschema``
     - Strict JSON schema validation of input files
     - 4.0
   * - ``matplotlib``
     - Distance-trajectory and metrics plots
     - 3.5
   * - ``scipy``
     - Truncated-normal size distributions
     - 1.8

Installation from source (recommended)
---------------------------------------

Clone the repository and install in **editable** mode so that local code
changes are immediately reflected::

    git clone https://github.com/daimoners/lignoforge.git
    cd lignoforge
    pip install -e .

Using a virtual environment (recommended)::

    python -m venv .venv
    source .venv/bin/activate     # Linux / macOS
    .venv\Scripts\activate        # Windows

    pip install -e .

Installation from PyPI
-----------------------

Once published on PyPI::

    pip install lignoforge

Conda environment
-----------------

::

    conda create -n lignoforge python=3.11
    conda activate lignoforge
    conda install -c conda-forge rdkit
    pip install -e .

.. note::
   RDKit is easiest to install via ``conda-forge``.  If installing with
   ``pip`` alone, use the official ``rdkit`` wheel::

       pip install rdkit

Verifying the installation
--------------------------

Run the built-in smoke test::

    python -c "from lignoforge.pipeline import LigninPipeline; print('OK')"

Or execute the demo script::

    python examples/demo_run.py --n-chains 2 --seed 42 --no-3d

A successful run prints a summary table and creates a ``demo_output/``
directory.

Upgrading
---------

Pull the latest commits and reinstall::

    git pull
    pip install -e .

Documentation dependencies
---------------------------

To build this documentation locally, install the extra requirements::

    pip install -r docs2/requirements.txt
    cd docs2
    make html

The built HTML is written to ``docs2/build/html/``.
