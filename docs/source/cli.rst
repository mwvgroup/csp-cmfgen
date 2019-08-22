Command Line Interface
======================

The analysis pipeline for this project provides a command line interface via
the *run_analysis.py* file. Subparsers are included for each of the analysis
steps.

.. argparse::
   :filename: ../run_analysis.py
   :func: create_parser
   :prog: run_analysis.py