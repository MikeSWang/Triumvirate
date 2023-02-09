# Configuration file for the Sphinx documentation builder.
#
# For a full list of options, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import inspect
import os
import sys
from configparser import ConfigParser
from datetime import datetime
from importlib import import_module
from pathlib import Path


sys.path.insert(0, os.path.abspath("../.."))


# -- Project information -------------------------------------------------

# Parse ``setup.cfg``.
setup_cfg_parser = ConfigParser()
setup_cfg_path = Path(
    inspect.getframeinfo(inspect.currentframe()).filename
).parent.parent.parent/"setup.cfg"

setup_cfg_parser.read(setup_cfg_path)
setup_cfg = dict(setup_cfg_parser.items('metadata'))

# Extract information from package.
pkg_name = setup_cfg.get('name').lower()
pkg_author = setup_cfg.get('author')

import_module(pkg_name)
pkg = sys.modules[pkg_name]

pkf_version = pkg.__version__
pkg_date = pkg.__date__.split('-').pop(0)

# Set fields.
project = pkg_name
author = pkg_author
release = pkf_version
if datetime.now().year == int(pkg_date):
    copyright = pkg_date
else:
    copyright = u'{0}\u2013{1}'.format(pkg_date, datetime.now().year)


# -- General configuration -----------------------------------------------

exclude_patterns = ['setup', 'config', 'tests', 'examples']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.coverage',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'sphinxcontrib.bibtex',
    'myst_nb',
]

master_doc = 'index'

pygments_style = 'sphinx'

source_suffix = ['.rst', '.txt', '.md']

templates_path = ['_templates']


# -- Options for HTML output ---------------------------------------------

html_favicon = '_static/Triumvirate.ico'

html_logo = '_static/Triumvirate.png'

html_static_path = ['', '_static/']

html_title = 'Documentation Home'  # u'\u200c'

# # Uncomment for 'furo'.
# html_theme = 'furo'

# Uncomment for 'sphinx_book_theme'.
html_theme = 'sphinx_book_theme'
html_theme_options = {
    'repository_url' : 'https://github.com/MikeSWang/Triumvirate',
    'home_page_in_toc': True,
    'logo_only': True,
    'toc_title': 'On this page',
    'use_download_button': False,
    'use_fullscreen_button': False,
    'use_repository_button': True,
}


# -- Extension configuration ---------------------------------------------

autodoc_member_order = 'bysource'

bibtex_bibfiles = ['_static/refs.bib']
bibtex_reference_style = 'author_year'

intersphinx_mapping = {
    'python': ("https://docs.python.org/3", None),
    'cython': ("https://cython.readthedocs.io/en/latest/", None),
    'pytest': ("https://docs.pytest.org/en/latest/", None),
    'numpy': ("https://numpy.org/doc/latest/", None),
    'astropy': ("https://docs.astropy.org/en/latest/", None),
    'nbodykit': ("https://nbodykit.readthedocs.io/en/latest/", None),
}

napoleon_include_special_with_doc = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

nb_execution_mode = 'off'

todo_include_todos = True
