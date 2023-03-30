"""Set up Triumvirate and its Cythonised extension modules.

"""
import logging
import os
import platform
from multiprocessing import cpu_count
from setuptools import setup
from setuptools.command.build_clib import build_clib

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension
from extension_helpers._openmp_helpers import check_openmp_support
from extension_helpers._setup_helpers import pkg_config


# ========================================================================
# Package
# ========================================================================

# -- Names ---------------------------------------------------------------

PROJECT_NAME = 'Triumvirate'

pkg_name = PROJECT_NAME.lower()


# -- Directories ---------------------------------------------------------

pkg_dir = os.path.join('src', pkg_name)

pkg_include_dir = os.path.join(pkg_dir, "include")
pkg_src_dir = os.path.join(pkg_dir, "src/modules")


# ========================================================================
# Commands
# ========================================================================

class BuildExt(build_ext):
    """Modified :class:`Cython.Distutils.build_ext`.

    """
    def finalize_options(self):
        """Modify parallelisation option in
        :meth:`Cython.Distutils.build_ext.build_ext.build_extensions`.

        """
        super().finalize_options()

        # Set unset parallel option to integer `num_procs` through
        # the environmental variable 'PY_BUILD_PARALLEL'.
        if self.parallel is None and num_procs is not None:
            self.parallel = num_procs

    def build_extensions(self):
        """Modify compiler and compilation configuration in
        :meth:`Cython.Distutils.build_ext.build_ext.build_extensions`.

        """
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')

        if isinstance(self.parallel, int):
            _num_procs = num_procs
        else:
            _num_procs = 1

        disutils_logger = logging.getLogger()
        disutils_logger.info(
            "running build_ext on %d processes with %s compiler (PY_CXX=%s)",
            _num_procs, self.compiler.compiler_so[0],
            os.environ.get('PY_CXX', '\'\'')
        )

        super().build_extensions()


class BuildClib(build_clib):
    """Modified :class:`setuptools.command.build_clib.build_clib`.

    """
    def build_libraries(self, libraries):
        """Modify compiler and compilation configuration in
        :meth:`setuptools.command.build_clib.build_clib.build_libraries`.

        """
        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')

        disutils_logger = logging.getLogger()
        disutils_logger.info(
            "running build_clib with %s compiler (PY_CXX=%s)",
            self.compiler.compiler_so[0],
            os.environ.get('PY_CXX', '\'\'')
        )

        super().build_libraries(libraries)


# ========================================================================
# Compilation
# ========================================================================

# -- Environment choices -------------------------------------------------

# Assume GNU compiler as default with matching OpenMP implementation.
COMPILER_CHOICES = {
    'default': 'g++',
    'linux': 'g++',
    'darwin': 'g++',
}
OPENMP_LIB_CHOICES = {
    'default': 'gomp',  # ''
    'linux': 'gomp',  # ''
    'darwin': '',  # 'omp'
}

build_platform = platform.system().lower()
if build_platform not in ['linux', 'darwin']:
    build_platform = 'default'


# -- Compiler ------------------------------------------------------------

# Specify language.
language = 'c++'

# Enforce platform-dependent default compiler choice.
compiler = os.environ.get('PY_CXX', COMPILER_CHOICES[build_platform])

os.environ['CC'] = compiler
os.environ['CXX'] = compiler
os.environ['CPP'] = compiler


# -- Options -------------------------------------------------------------

# Specify compilation and linker flags.
cflags = os.environ.get('PY_CFLAGS', '').split()
ldflags = [
    _ldflag for _ldflag in os.environ.get('PY_LDFLAGS', '').split()
    if not _ldflag.startswith('-l')  # set in `ext_libs` instead
]

# Add optional flags for OpenMP support (enabled by default).
disable_omp = os.environ.get('PY_NO_OMP')
if disable_omp is None and check_openmp_support():
    # Enforce platform-dependent default OpenMP library choice.
    default_libomp = OPENMP_LIB_CHOICES[build_platform]
    default_ldflags_omp = '-l' + default_libomp if default_libomp else ''

    ldflags_omp = os.environ.get('PY_LDOMP', default_ldflags_omp).split()

    # Ensure OpenMP options are enabled (by default).
    for _cflag in ['-fopenmp', '-DTRV_USE_OMP', '-DTRV_USE_FFTWOMP']:
        if _cflag not in cflags:
            cflags.append(_cflag)

    # Ensure OpenMP-related libraries are linked (by default).
    for _ldflag in ['-lfftw3_omp',] + ldflags_omp:  # noqa: E231
        if _ldflag not in ldflags:
            ldflags.append(_ldflag)

# Set macros.
pkg_macros = [
    ('TRV_EXTCALL', None),
    # ('TRV_USE_LEGACY_CODE', None),
    # ('DBG_MODE', None),
    # ('DBG_FLAG_NOAC', None),
]

npy_macros = [
    ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
]

macros = pkg_macros + npy_macros

# Set include directories.
npy_include = numpy.get_include()

ext_includes = [
    incl
    for incl in os.environ.get('PY_INCLUDES', "").replace("-I", "").split()
    if pkg_dir not in incl
]

includes = [pkg_include_dir, npy_include,] + ext_includes  # noqa: E231

# Set libraries.
libdirs = []

libs = ['gsl', 'gslcblas', 'm', 'fftw3',]  # noqa: E231

# Check for missing options using :mod:`extension_helpers`.
checked_options = pkg_config(libs, default_libraries=[])
for incl_ in checked_options['include_dirs']:
    if incl_ not in includes:
        includes.append(incl_)
for lib_ in checked_options['library_dirs']:
    if lib_ not in libdirs:
        libdirs.append(lib_)
for flag_ in checked_options['extra_compile_args']:
    if flag_ not in cflags and flag_ not in ldflags:
        cflags.append(flag_)


# ========================================================================
# Build
# ========================================================================

# Set and prepend package as (a) static library/libraries.
pkg_lib = 'trv'

pkg_library = (
    pkg_lib,
    {
        'sources': [
            os.path.join(pkg_src_dir, _cpp_source)
            for _cpp_source in os.listdir(pkg_src_dir)
        ],
        'macros': macros,
        'cflags': cflags,
        'include_dirs': includes,
    }
)

# The linking order matters, and to ensure that, `-ltrv` is duplicated
# also in `libraries` keyword argument in :meth:`setuptools.setup`.
libs = [pkg_lib,] + libs  # noqa: E231

# Add optional parallelisation for build jobs.
flag_parallel = os.environ.get('PY_BUILD_PARALLEL', '').strip()
if '-j' == flag_parallel:
    num_procs = cpu_count()
elif '-j' in flag_parallel:
    try:
        num_procs = int(flag_parallel.lstrip('-j'))
    except ValueError:
        num_procs = None
else:
    num_procs = None


# ========================================================================
# Extensions
# ========================================================================

def define_extension(module_name,
                     auto_cpp_source=False, extra_cpp_sources=None,
                     **ext_kwargs):
    """Return an extension from given source files and options.

    Parameters
    ----------
    module_name : str
        Module name as in the built extension module
        ``<pkg_name>.<module_name>``. This sets the .pyx source file to
        ``<pkg_dir>/<module_name>.pyx``.
    auto_cpp_source : bool, optional
        If `True` (default is `False`), add
        ``<pkg_src_dir>/<module_name>.cpp`` to the list of source files
        with any prefix underscores removed from ``<module_name>``.
    extra_cpp_sources : list of str, optional
        Additional .cpp source files under ``<pkg_src_dir>``
        (default is `None`).
    **ext_kwargs
        Options to pass to :class:`Cython.Distutils.Extension`.

    Returns
    -------
    :class:`Cython.Distutils.Extension`
        Cython extension.

    """
    ext_module_kwargs = {
        'define_macros': macros,
        'include_dirs': includes,
        'library_dirs': libdirs,
        'libraries': libs,
        'extra_compile_args': cflags,
        'extra_link_args': ldflags,
    }

    if extra_cpp_sources is None:
        extra_cpp_sources = []

    if ext_kwargs:
        ext_module_kwargs.update(ext_kwargs)

    sources = [os.path.join(pkg_dir, f"{module_name}.pyx"),]  # noqa: E231
    if auto_cpp_source:
        sources.append(
            os.path.join(pkg_src_dir, f"{module_name}.cpp".lstrip('_'))
        )
    for _cpp_source in extra_cpp_sources:
        sources.append(os.path.join(pkg_src_dir, _cpp_source))

    return Extension(
        f'{pkg_name}.{module_name}',
        sources=sources, language=language, **ext_module_kwargs
    )


# Cythonise extension modules.
ext_module_configs = {
    'parameters': {},
    'dataobjs': {},
    '_particles': {},
    '_twopt': {},
    '_threept': {},
    '_fftlog': {},
}

extensions = [
    define_extension(name, **cfg)
    for name, cfg in ext_module_configs.items()
]

cython_directives = {
   'language_level': '3',
   'c_string_encoding': 'utf-8',
   'embedsignature': True,
}

cython_ext_modules = cythonize(
    extensions,
    compiler_directives=cython_directives,
    nthreads=num_procs,
)


# ========================================================================
# Set-up
# ========================================================================

if __name__ == '__main__':
    setup(
        use_scm_version={
            'version_scheme': 'post-release',
            # 'local_scheme': 'no-local-version',
        },
        cmdclass={
            'build_clib': BuildClib,
            'build_ext': BuildExt,
        },
        ext_modules=cython_ext_modules,
        libraries=[pkg_library,],  # noqa: E231
    )
