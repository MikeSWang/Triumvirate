# Configuration file for the Sphinx documentation builder.
#
# For a full list of options, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('../../'))


# -- Project information -----------------------------------------------------

project = 'Triumvirate'
copyright = ' 2022, Mike S Wang & Naonori Sugiyama'
author = 'Mike S Wang, Naonori Sugiyama'
release = 'dev0.0'


# -- General configuration ---------------------------------------------------

exclude_patterns = ['setup', 'config', 'tests', 'examples']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'nbsphinx',
]

master_doc = 'index'

pygments_style = 'sphinx'

source_suffix = ['.rst', '.txt', '.md']

templates_path = ['_templates']


# -- Options for HTML output -------------------------------------------------

html_logo = '_static/Triumvirate.png'

html_sidebars = {
    '**': ['navigation.html', 'searchbox.html'],
    'using/windows': ['windowssidebar.html', 'searchbox.html'],
}

html_static_path = ['_static']

html_theme = 'alabaster'

html_theme_options = {
    'fixed_sidebar' : True,
    'github_repo': 'Triumvirate',
    'github_user': 'MikeSWang',
    'page_width': '1150px',
    'sidebar_width': '250px',
}


# -- Extension configuration -------------------------------------------------

autodoc_member_order = 'bysource'
autosummary_generate = True

intersphinx_mapping = {
    'python': ("https://docs.python.org/3", None),
    'pytest': ("https://docs.pytest.org/en/latest/", None),
    'numpy': ("https://numpy.org/doc/stable/", None),
    'scipy': ("https://docs.scipy.org/doc/scipy/reference", None),
    'matplotlib': ("https://matplotlib.org", None),
}

napoleon_include_special_with_doc = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

todo_include_todos = True
