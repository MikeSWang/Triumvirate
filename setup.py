"""Set up Triumvirate and its Cythonised extension modules.

"""
import os
import platform
from distutils.sysconfig import get_config_vars
from distutils.util import convert_path
from setuptools import setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy


PKG_NAME = 'Triumvirate'


# -- Repository ----------------------------------------------------------

pkg_dir = pkg_name = PKG_NAME.lower()

# Extract package information.
pkg_info = {}
with open(convert_path(os.path.join(pkg_dir, "__init__.py"))) as pkg_pyfile:
    exec(pkg_pyfile.read(), pkg_info)

# Determine repository branch.
version = pkg_info.get('__version__')
if any([segment in version for segment in ['a', 'b', 'rc']]):
    branch = 'main'
elif 'dev' in version:
    branch = 'dev'
else:
    branch = 'stable'


# -- Compilation ---------------------------------------------------------

# Specify language.
language = 'c++'

# Suppress irrelevant compilation warnings.
config_vars = get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')

# Enforce compiler choice.
if platform.system() == 'Linux':
    cxx_default = 'g++'
elif platform.system() == 'Darwin':
    cxx_default = 'clang++'
else:
    cxx_default = 'g++'

os.environ['CC'] = os.environ.get('PY_CXX', cxx_default)
os.environ['CXX'] = os.environ.get('PY_CXX', cxx_default)

# Modify compilation options.
cflags = os.environ.get('PY_CFLAGS', '').split() + ['-std=c++11',]
ldflags = [
    lib for lib in os.environ.get('PY_LDFLAGS', '').split()
    if not lib.startswith('-l')
]


# -- Extensions ----------------------------------------------------------

def return_extension(module_name, extra_cpp_sources, **ext_kwargs):
    """Return an extension from given source files and options.

    Parameters
    ----------
    module_name : str
        Module name as in the built extension module
        ``<pkg_name>.<module_name>``. This sets the .pyx source file to
        ``<pkg_dir>/<module_name>.pyx`` and the .cpp source file to
        ``<pkg_dir>/src/modules/<module_name>.cpp``.
    extra_cpp_sources : list of str
        Additional .cpp source files under ``<pkg_dir>/src/modules``.
    **ext_kwargs
        Options to pass to :class:`Cython.Distutils.Extension`.

    Returns
    -------
    :class:`Cython.Distutils.Extension`
        Cython extension.

    """
    ext_module_kwargs = {
        'include_dirs': includes,
        'extra_compile_args': cflags,
        'extra_link_args': ldflags,
        'define_macros': macros,
    }

    if ext_kwargs:
        ext_module_kwargs.update(ext_kwargs)

    sources = [
        os.path.join(pkg_dir, f"{module_name}.pyx"),
        os.path.join(pkg_src, f"{module_name.lstrip('_')}.cpp"),
    ]
    for cpp_source in extra_cpp_sources:
        sources.append(os.path.join(pkg_src, cpp_source))

    return Extension(
        f'{pkg_name}.{module_name}',
        sources=sources, language=language, **ext_module_kwargs
    )


# Set source, include and library paths.
pkg_src = os.path.join(pkg_dir, "src/modules")

pkg_include = os.path.join(pkg_dir, "include")

npy_include = numpy.get_include()

ext_includes = [
    incl for incl in os.environ.get('PY_INCLUDES', "").replace("-I", "").split()
    if pkg_dir not in incl
]

ext_libraries = ['gsl', 'gslcblas', 'm', 'fftw3', 'fftw3_omp',]

includes = [pkg_include, npy_include,] + ext_includes
libraries = ext_libraries

# Set macros.
pkg_macros = [
    ('TRV_EXTCALL', None),
    # ('TRV_USE_LEGACY_CODE', None),
    # ('DBG_MODE', None),
    # ('DBG_NOAC', None),
]

npy_macros = [
    ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
]

macros = pkg_macros + npy_macros

# Define extension modules.
module_config = {
    'parameters': {'extra_cpp_sources': ["monitor.cpp",]},
    'dataobjs': {'extra_cpp_sources': ["monitor.cpp", "parameters.cpp",]},
    '_particles': {'extra_cpp_sources': ["monitor.cpp",]},
    '_twopt': {
        'extra_cpp_sources': [
            "monitor.cpp",
            "maths.cpp",
            "parameters.cpp",
            "dataobjs.cpp",
            "particles.cpp",
            "field.cpp",
        ],
        'libraries': libraries,
    },
    '_threept': {
        'extra_cpp_sources': [
            "monitor.cpp",
            "maths.cpp",
            "parameters.cpp",
            "dataobjs.cpp",
            "particles.cpp",
            "field.cpp",
            "twopt.cpp",
        ],
        'libraries': libraries,
    },
    '_fftlog': {
        'extra_cpp_sources': ["monitor.cpp", "maths.cpp", "arrayops.cpp",],
        'libraries': libraries,
    },
}

cython_directives = {
   'language_level': '3',
   'c_string_encoding': 'utf-8',
}

cython_modules = [
    return_extension(key, **val)
    for key, val in module_config.items()
]


# -- Set-up --------------------------------------------------------------

if __name__ == '__main__':
    setup(
        license=pkg_info.get('__license__'),
        author=pkg_info.get('__author__'),
        maintainer=pkg_info.get('__maintainer__'),
        maintainer_email=pkg_info.get('__maintainer_email__'),
        description=pkg_info.get('__description__'),
        cmdclass={'build_ext': build_ext},
        ext_modules=cythonize(
            cython_modules, compiler_directives=cython_directives
        )
    )
