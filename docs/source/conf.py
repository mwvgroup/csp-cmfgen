#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Configuration file for the Sphinx documentation builder."""

# -- Project information -----------------------------------------------------

project = 'SNe Constraints'
copyright = '2019, Daniel Perrefort'
author = 'Daniel Perrefort'
version, release = '', ''

# -- General configuration ---------------------------------------------------

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.viewcode',
    'sphinx.ext.napoleon',
    'sphinxarg.ext',
    'sphinx_copybutton'
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames. (source_suffix = ['.rst', '.md'])
source_suffix = '.rst'
exclude_patterns = []

# The master toctree document.
master_doc = 'index'

# -- Options for HTML output -------------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_theme_options = {}
html_static_path = ['../_static']

# -- Options for HTMLHelp output ---------------------------------------------

# Output file base name for HTML help builder.
htmlhelp_basename = 'nir_analysisdoc'

# -- Options for LaTeX output ------------------------------------------------

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    # 'papersize': 'letterpaper',

    # The font size ('10pt', '11pt' or '12pt').
    # 'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    # 'preamble': '',

    # Latex figure (float) alignment
    # 'figure_align': 'htbp',
}

latex_documents = [
    (master_doc, 'nir_analysis.tex', 'nir\\_analysis Documentation',
     'Daniel Perrefort', 'manual'),
]
