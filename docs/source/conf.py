"""Configure the Sphinx documentation builder.

"""
import inspect
import os
import sys
import textwrap
from configparser import ConfigParser
from datetime import datetime
from importlib import import_module
from pathlib import Path


sys.path.insert(0, os.path.abspath("../../src"))


# -- Project information -------------------------------------------------

# Directories and paths.
root_dir = Path(
    inspect.getframeinfo(inspect.currentframe()).filename
).parent.parent.parent

docs_dir = root_dir/"docs"

# Parse ``setup.cfg``.
setup_cfg_parser = ConfigParser()
setup_cfg_path = root_dir/"setup.cfg"

setup_cfg_parser.read(setup_cfg_path)
setup_cfg = dict(setup_cfg_parser.items('metadata'))

# Extract information from package.
pkg_name = setup_cfg.get('name').lower()
pkg_author = ' &'.join(setup_cfg.get('author').split(',', 1))

import_module(pkg_name)
pkg = sys.modules[pkg_name]

pkg_date = pkg.__date__.split('-').pop(0)
if os.environ.get('READTHEDOCS') == 'True':
    with open(docs_dir/"RTD_VERSION.tmp", 'r') as vers_file:
        pkg_version = vers_file.readline().strip()
else:
    pkg_version = pkg.__version__

# Set fields.
project = pkg_name
# author = pkg_author
release = pkg_version
if datetime.now().year == int(pkg_date):
    copyright = '{0}, {1}'.format(pkg_date, pkg_author)
else:
    copyright = u'{0}\u2013{1}, {2}'.format(
        pkg_date, datetime.now().year, pkg_author
    )


# -- General configuration -----------------------------------------------

exclude_patterns = ['setup', 'config', 'tests', 'examples']

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.coverage',
    'sphinx.ext.githubpages',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx.ext.graphviz',
    'sphinx.ext.todo',
    'sphinx.ext.viewcode',
    'breathe',
    'exhale',
    'myst_nb',
    'sphinxcontrib.bibtex',
    'sphinx_copybutton',
    'sphinx_tabs.tabs',
    'sphinx_togglebutton',
    # 'sphinx_inline_tabs',
]

master_doc = 'index'

pygments_style = 'sphinx'

source_suffix = ['.rst', '.md', '.ipynb']

templates_path = ['_templates']


# -- Options for HTML output ---------------------------------------------

html_favicon = '_static/Triumvirate.ico'

html_logo = '_static/Triumvirate.png'

html_static_path = ['', '_static/']

html_title = 'Triumvirate Documentation'  # u'\u200c'

html_last_updated_fmt = '%Y-%m-%d at %H:%M:%S. Version {}'.format(release)

html_theme = 'sphinx_book_theme'
html_theme_options = {
    'path_to_docs': "docs/source",
    'repository_url': 'https://github.com/MikeSWang/Triumvirate',
    'toc_title': 'On this page',
    # 'extra_footer': "<div>{}</div>".format(None),
    'home_page_in_toc': False,
    'use_download_button': False,
    'use_fullscreen_button': False,
    'use_repository_button': False,
    'use_source_button': True,
    'icon_links': [
        {
            'name': 'GitHub',
            'url': "https://github.com/MikeSWang/Triumvirate",
            'icon': 'fa-brands fa-github',
        },
        {
            'name': 'PyPI',
            'url': "https://pypi.org/project/Triumvirate",
            'icon': (
                "https://img.shields.io/pypi/v/Triumvirate"
                "?logo=PyPI&color=informational"
            ),
            'type': 'url',
        },
        {
            'name': 'conda',
            'url': "https://anaconda.org/msw/triumvirate",
            'icon': (
                "https://img.shields.io/conda/vn/msw/triumvirate"
                "?logo=Anaconda&color=informational"
            ),
            'type': 'url'
        },
    ]
}


# -- Extension configuration ---------------------------------------------

autodoc_member_order = 'bysource'

autosummary_generate = True

bibtex_bibfiles = ['_static/refs.bib']
bibtex_reference_style = 'author_year'

breathe_projects = {'Triumvirate': "apidoc_cpp/xml/"}
breathe_default_project = 'Triumvirate'
breathe_implementation_filename_extensions = ['.cpp']

exhale_args = {
    'containmentFolder': "./apidoc_cpp",
    'rootFileName': "apidoc_cpp.rst",
    'rootFileTitle': "C++ Library",
    'createTreeView': True,
    'treeViewIsBootstrap': True,
    'exhaleExecutesDoxygen': True,
    'exhaleUseDoxyfile': True,
    'doxygenStripFromPath': "..",
    'fullToctreeMaxDepth': 1,
    'afterTitleDescription': textwrap.dedent("""
        .. seealso::

            `Doxygen version <../apiref_doxy/index.html>`_
            of the C++ API reference.
    """),
}

intersphinx_mapping = {
    'python': ("https://docs.python.org/3", None),
    'numpy': ("https://numpy.org/doc/stable/", None),
    'scipy': ("https://docs.scipy.org/doc/scipy", None),
    'astropy': ("https://docs.astropy.org/en/latest/", None),
    'nbodykit': ("https://nbodykit.readthedocs.io/en/latest/", None),
}

myst_enable_extensions = [
    'amsmath',
    'colon_fence',
    'deflist',
    'dollarmath',
    'fieldlist',
    'html_admonition',
    'html_image',
    'replacements',
    'substitution',
    'tasklist',
]

myst_dmath_double_inline = True

myst_substitutions = {
    'Triumvirate':
        "<span style=\"font-variant: small-caps\">Triumvirate</span>",
}

napoleon_include_special_with_doc = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

nb_execution_mode = 'off'

todo_include_todos = True
