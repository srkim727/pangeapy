# Configuration file for the Sphinx documentation builder.

import os
import sys
sys.path.insert(0, os.path.abspath('..'))  # make pangeapy importable

project = 'pangeapy'
author = 'Seongryong Kim'
release = '0.1.0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',  # support for Google/NumPy style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.autosummary',
    'sphinx.ext.todo',
    'nbsphinx',
    'myst_parser',
    'sphinx_copybutton',
    'sphinx_design',
]

# Notebook build
nbsphinx_allow_errors = True
nbsphinx_execute = 'never'  # don't run notebooks at build time

# Theme
html_theme = 'furo'

# Options
autosummary_generate = True
autodoc_typehints = 'description'
napoleon_google_docstring = True
napoleon_numpy_docstring = True

# Paths
templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
