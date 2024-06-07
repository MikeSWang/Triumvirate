"""Configure the Sphinx documentation builder.

"""
import inspect
import os
import sys
import textwrap
from configparser import ConfigParser
from datetime import datetime
from importlib import import_module
from itertools import chain
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

# Set custom fields.
anaconda_user = 'msw'
github_user = 'MikeSWang'
github_repo = pkg_name.capitalize()
proj_name = github_repo
if 'dev' in pkg_version:
    doc_version = 'dev'
else:
    doc_version = pkg_version


# -- General configuration -----------------------------------------------

exclude_patterns = ['setup', 'config', 'tests', 'examples']

extensions = [
    'breathe',
    'exhale',
    'myst_nb',
    'hoverxref.extension',
    'notfound.extension',
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
    'sphinx_copybutton',
    'sphinx_design',
    # 'sphinx_inline_tabs',
    'sphinx_search.extension',
    # 'sphinx_tabs.tabs',
    'sphinx_togglebutton',
    'sphinxcontrib.bibtex',
]

master_doc = 'index'

pygments_style = 'sphinx'

source_suffix = ['.rst', '.md', '.ipynb']

templates_path = ['_templates']


# -- Options for HTML output ---------------------------------------------

ANACONDA_CHANNEL_URL = f"https://anaconda.org/{anaconda_user}/{pkg_name}"

WINCONV_MOD_URL = (
    "https://triumvirate.readthedocs.io/en/latest/"
    "apidoc_py/triumvirate.winconv.html"
)

html_context = {
    'github_user': github_user,
    'github_repo': github_repo,
    'github_version': 'main',
    'doc_path': "docs/source",
}  # theme: pydata_sphinx_theme

html_css_files = ["css/custom.css"]  # theme: pydata_sphinx_theme

html_favicon = f"_static/image/{proj_name}.ico"

html_js_files = ["js/custom-icons.js"]  # theme: pydata_sphinx_theme

# html_logo = f'_static/image/{proj_name}.png'  # theme: sphinx_book_theme
html_logo = \
    f'_static/image/{proj_name}-ls.png'         # theme: pydata_sphinx_theme

html_show_sourcelink = False  # theme: pydata_sphinx_theme

html_static_path = ["", "_static/"]

html_title = f'{proj_name} Documentation'  # u'\u200c'

html_last_updated_fmt = "%a, %Y-%m-%dT%H:%M:%S%z. Version {}".format(release)

html_sidebars = {
    "**": [
        'sidebar-nav-bs',
    ],
    # "background": [],
    # "installation": [],
    # "releases": [],
    # "faq": [],
}  # theme: pydata_sphinx_theme

html_theme = 'pydata_sphinx_theme'  # 'sphinx_book_theme'  #
html_theme_options = {
    'announcement': (
        "&#128064; "
        "Version 0.4 to be released soon (expected in 2024) with "
        f"<a href=\"{WINCONV_MOD_URL}\">window convolution</a> "
        "functionality!"
        "<br>"
        "You can try out the developmental release "
        f"<a href=\"{ANACONDA_CHANNEL_URL}/files?channel=dev\">"
        "<code>msw/label/dev::triumvirate</code>"
        "</a> "
        "from Anaconda."
    ),
    'back_to_top_button': True,         # theme: pydata_sphinx_theme
    # 'path_to_docs': "docs/source",    # theme: sphinx_book_theme
    # 'repository_url': f'https://github.com/{github_user}/{github_repo}',
    #                                   # theme: sphinx_book_theme
    # 'toc_title': 'On this page',      # theme: sphinx_book_theme
    # 'extra_footer': "<div>{}</div>".format(None),
    'footer_start': ['last-updated',],  # theme: pydata_sphinx_theme
    'footer_center': [],                # theme: pydata_sphinx_theme
    'footer_end': ['copyright',],       # theme: pydata_sphinx_theme
    'header_links_before_dropdown': 6,  # theme: pydata_sphinx_theme
    # 'home_page_in_toc': False,        # theme: sphinx_book_theme
    'navbar_align': 'left',             # theme: pydata_sphinx_theme
    'navbar_start': [
        'navbar-logo',
    ],                                  # theme: pydata_sphinx_theme
    'navbar_center': [
        'version-switcher',
        'navbar-nav',
    ],                                  # theme: pydata_sphinx_theme
    'navbar_persistent': [
        'search-button-field',
    ],                                  # theme: pydata_sphinx_theme
    'navbar_end': [
        'theme-switcher',
        'navbar-icon-links',
    ],                                  # theme: pydata_sphinx_theme
    'navigation_with_keys': False,      # theme: pydata_sphinx_theme
    'switcher': {
        'json_url': (
            "https://raw.githubusercontent.com/"
            "{github_user}/{github_repo}/{github_version}/"
            "docs/versions.json"
        ).format(**html_context),
        'version_match': doc_version,
    },                                  # theme: pydata_sphinx_theme
    # 'use_download_button': False,     # theme: sphinx_book_theme
    # 'use_fullscreen_button': False,   # theme: sphinx_book_theme
    # 'use_repository_button': False,   # theme: sphinx_book_theme
    # 'use_source_button': True,        # theme: sphinx_book_theme
    'icon_links': [
        {
            'name': 'GitHub',
            'url': f"https://github.com/{github_user}/{github_repo}",
            'icon': 'fa-brands fa-github',
        },
        # {
        #     'name': 'PyPI',
        #     'url': f"https://pypi.org/project/{proj_name}",
        #     'icon': (
        #         f"https://img.shields.io/pypi/v/{proj_name}"
        #         "?logo=PyPI&color=informational"
        #     ),
        #     'type': 'url',
        # },
        {
            'name': 'PyPI',
            'url': f"https://pypi.org/project/{proj_name}",
            'icon': 'fa-custom fa-pypi',
        },
        # {
        #     'name': 'Conda',
        #     'url': f"https://anaconda.org/{anaconda_user}/{pkg_name}",
        #     'icon': (
        #         "https://img.shields.io/conda/vn/"
        #         f"{anaconda_user}/{pkg_name}"
        #         "?logo=Anaconda&color=informational"
        #     ),
        #     'type': 'url'
        # },
        {
            'name': 'Conda',
            'url': f"https://anaconda.org/{anaconda_user}/{pkg_name}",
            'icon': 'fa-custom fa-anaconda',
        },
    ]
}


# -- Extension configuration ---------------------------------------------

autodoc_member_order = 'bysource'

TYPE_ALIASES_BY_MODULE = {
    'triumvirate.wincov': [
        'DegreeType', 'MultipoleLike', 'Multipole3PtLike',
        'CoeffType', 'WinConvTermLike', 'FormulaeDictType',
    ],
}

autodoc_type_aliases = {
    ta: f"{mod_name}.{ta}"
    for mod_name in TYPE_ALIASES_BY_MODULE.keys()
    for ta in TYPE_ALIASES_BY_MODULE.get(mod_name, [])
}

autosummary_generate = True

bibtex_bibfiles = ["_static/text/refs.bib"]
bibtex_reference_style = 'author_year'

breathe_projects = {proj_name: "apidoc_cpp/xml/"}
breathe_default_project = proj_name
breathe_implementation_filename_extensions = ['.cpp']

exhale_args = {
    'containmentFolder': "./apidoc_cpp",
    'rootFileName': "apidoc_cpp.rst",
    'rootFileTitle': "C++ Library",
    'createTreeView': True,
    'treeViewIsBootstrap': False,  # HACK: unknown issue with `exhale`
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

hoverxref_auto_ref = True
hoverxref_mathjax = True
hoverxref_sphinxtabs = True
hoverxref_domains = [
    'cpp',
    'py',
]
hoverxref_roles = [
    'ref',  # 'doc',
    # 'package', 'mod',
    'obj', 'exc', 'const', 'data',
    'class', 'func', 'meth', 'attr',
]
hoverxref_role_types = {role: 'tooltip' for role in hoverxref_roles}
hoverxref_intersphinx = [
    'numpy',
    'sympy',
]

intersphinx_mapping = {
    'python': ("https://docs.python.org/3", None),
    'numpy': ("https://numpy.org/doc/stable/", None),
    'scipy': ("https://docs.scipy.org/doc/scipy", None),
    'sympy': ("https://docs.sympy.org/latest/", None),
    'mpi4py': ("https://mpi4py.readthedocs.io/en/latest/", None),
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
    proj_name: f"<span style=\"font-variant: small-caps\">{proj_name}</span>",
}

napoleon_include_special_with_doc = True
napoleon_google_docstring = False
napoleon_use_param = False
napoleon_use_ivar = True

nb_execution_mode = 'off'

notfound_context = {
    'title': "Page not found",
    'body': (
        "<h1>Have you been sniffing around?</h1>\n\n"
        "Like many a scientist's dream, this page is lost&hellip; "
        "or has never existed."
    ),
}

todo_include_todos = True


# -- Application set-up --------------------------------------------------

def skip_cdef_member(app, what, name, obj, skip, options):

    SKIP_LISTS = {
        'dataobjs.Binning': [
            'bin_centres', 'bin_edges', 'bin_max', 'bin_min', 'bin_widths',
            'num_bins', 'scheme', 'space',
        ],
    }

    skip_list = list(chain.from_iterable(SKIP_LISTS.values()))

    return True if name in skip_list else None


def setup(app):

    app.connect('autodoc-skip-member', skip_cdef_member)

    return {
        'parallel_read_safe': True,
        'parallel_read_unsafe': False,
        'parallel_write_safe': True,
        'parallel_write_unsafe': False,
    }
