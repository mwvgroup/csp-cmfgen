Photometric Constraints on SNe Models with LSST
===============================================

.. |python| image:: https://img.shields.io/badge/python-3.6+-blue.svg
    :target: #

.. |license| image:: https://img.shields.io/badge/license-GPL%20v3.0-blue.svg
    :target: https://www.gnu.org/licenses/gpl-3.0.en.html

.. |travis| image:: https://www.travis-ci.com/mwvgroup/csp-cmfgen.svg?branch=master
   :target: https://www.travis-ci.com/mwvgroup/csp-cmfgen

.. |docs| image:: https://readthedocs.org/projects/csp-cmfgen/badge/?version=latest
   :target: https://csp-cmfgen.readthedocs.io/en/latest/index.html

.. rst-class:: badges

   +------------------------------------+
   | |python| |license| |travis| |docs| |
   +------------------------------------+

This documentation provides a technical overview for the comparison of CMFGEN supernova models against observations
taken by the Carnegie Supernova Project.

Using This Documentation
------------------------

The **Project Notes** section outlines the goals, technical decisions, and research notes taken by the research team.
The **Getting Started** guides outline the installation of the project's source code and how to reproduce the scientific
results. API documentation for the project's source code is broken down into two sections. The **Backend Utilities**
section documents modules that are used to manipulate data across the various stages of the analysis. The
**Science Modules** section documents Python modules used to calculate and tabulate the scientific deliverables of
this project.

Source Code and Repository Structure
------------------------------------

Source code for this project can be found online via `GitHub`_. We provide summaries below for key files in the
project repository::

    project parent/
    ├── analysis/
    ├── ascii_models/
    ├── config/
    ├── docs/
    ├── notebooks/
    ├── results/
    ├── tests/
    |
    ├── README.md
    ├── environment.yml
    ├── run_analysis.py
    ├── run_analysis.sh
    └── setup.py


+----------------------+------------------------------------------------------------------------------+
| File / Directory     | Description                                                                  |
+======================+==============================================================================+
| *analysis/*          | Python source code for the analysis performed by this project.               |
+----------------------+------------------------------------------------------------------------------+
| *ascii_models/*      | A copy of the CMFGEN models in ascii format. A Python script is also         |
|                      | provided for converting the models to `npz` format.                          |
+----------------------+------------------------------------------------------------------------------+
| *config/*            | Configuration files (and the script that generates them) that defined priors |
|                      | and arguments used when fitting light-curves.                                |
+----------------------+------------------------------------------------------------------------------+
| *docs/*              | Source code for the project's documentation.                                 |
+----------------------+------------------------------------------------------------------------------+
| *notebooks/*         | Jupyter notebooks that inspect results of the analysis.                      |
+----------------------+------------------------------------------------------------------------------+
| *results/*           | Tabulated results from the ``analysis`` package.                             |
+----------------------+------------------------------------------------------------------------------+
| *tests/*             | Tests for the ``analysis`` package.                                          |
+----------------------+------------------------------------------------------------------------------+
| *README.md*          | Github Landing document that points readers to the online documentation.     |
+----------------------+------------------------------------------------------------------------------+
| *environment.yml*    | Installation requirements for the ``analysis`` package.                      |
+----------------------+------------------------------------------------------------------------------+
| *run_analysis.py*    | A command line interface for running the ``analysis`` package.               |
+----------------------+------------------------------------------------------------------------------+
| *run_analysis.sh*    | Runs the command line interface for various combinations of arguments        |
+----------------------+------------------------------------------------------------------------------+
| *setup.py*           | Installation script for the ``analysis`` package.                            |
+----------------------+------------------------------------------------------------------------------+

.. _GitHub: https://github.com/mwvgroup/csp-cmfgen/

.. toctree::
   :hidden:
   :maxdepth: 1

   Overview<self>

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Project Notes:

   project_notes/project_outline

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Getting Started:

   getting_started/installation
   getting_started/cli


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Backend Utilities:

   backend_utils/overview
   backend_utils/models
   backend_utils/regression
   backend_utils/utils
   backend_utils/exceptions

.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: Science Modules:

   science_modules/overview
   science_modules/spectra_chisq
   science_modules/equivalent_width
   science_modules/band_fitting
   science_modules/lc_colors
