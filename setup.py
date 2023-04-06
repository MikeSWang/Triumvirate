"""Set up Triumvirate and its Cythonised extension modules.

"""
import logging
import os
import platform
import sys
from copy import deepcopy
from multiprocessing import cpu_count
from setuptools import setup
from setuptools.command.build_clib import build_clib

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension
from extension_helpers._openmp_helpers import check_openmp_support
from extension_helpers._setup_helpers import pkg_config


distutils_logger = logging.getLogger()  # recycled from :mod:`distutils`


# ========================================================================
# Package
# ========================================================================

PROJECT_NAME = 'Triumvirate'


def get_pkg_name():
    """Get package name in lower-case format.

    Returns
    -------
    str
        Package name.

    """
    return PROJECT_NAME.lower()


def get_pkg_dir(topdir="src"):
    """Get package directory.

    Parameters
    ----------
    topdir : str, optional
        Top directory.

    Returns
    -------
    str
        Package directory.

    """
    return os.path.join(topdir, get_pkg_name())


def get_pkg_include_dir(subdir="include"):
    """Get package C++ header directory.

    Parameters
    ----------
    subdir : str, optional
        Subdirectory of C++ headers in the package directory.

    Returns
    -------
    str
        C++ header directory.

    """
    return os.path.join(get_pkg_dir(), subdir)


def get_pkg_src_dir(subdir="src/modules"):
    """Get package C++ source directory.

    Parameters
    ----------
    subdir : str, optional
        Subdirectory of C++ sources in the package directory.

    Returns
    -------
    str
        C++ source directory.

    """
    return os.path.join(get_pkg_dir(), subdir)


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

        _num_procs = get_build_num_procs()
        if self.parallel is None and _num_procs is not None:
            self.parallel = _num_procs

    def build_extensions(self):
        """Modify compiler and compilation configuration in
        :meth:`Cython.Distutils.build_ext.build_ext.build_extensions`.

        """
        OPTS_TO_REMOVE = ['-Wstrict-prototypes',]  # noqa: E231
        for opt in OPTS_TO_REMOVE:
            try:
                self.compiler.compiler_so.remove(opt)
            except ValueError:
                pass

        try:
            _num_procs = max(int(self.parallel), 1)
        except TypeError:
            _num_procs = 1

        distutils_logger.info(
            "running build_ext on %d process(es) with %s compiler",
            _num_procs, self.compiler.compiler_so[0],
        )

        super().build_extensions()


class BuildClib(build_clib):
    """Modified :class:`setuptools.command.build_clib.build_clib`.

    """
    def build_libraries(self, libraries):
        """Modify compiler and compilation configuration in
        :meth:`setuptools.command.build_clib.build_clib.build_libraries`.

        """
        OPTS_TO_REMOVE = ['-Wstrict-prototypes',]  # noqa: E231
        for opt in OPTS_TO_REMOVE:
            try:
                self.compiler.compiler_so.remove(opt)
            except ValueError:
                pass

        distutils_logger.info(
            "running build_clib with %s compiler",
            self.compiler.compiler_so[0],
        )

        super().build_libraries(libraries)


def get_build_num_procs():
    """Get build parallel processes.

    Returns
    -------
    int or None
        Number of parallel processes, which may be `None` in the case
        of no parallelisation.

    """
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

    return num_procs


def prioprint(*args, **kwargs):
    """`print` in high priority.

    """
    print(*args, file=sys.stderr, **kwargs)


def display_py_environs():
    """Display customised environment variables passed to setup.

    """
    PY_ENV_VARS = [
        'PY_CXX',
        'PY_CXXFLAGS',
        'PY_LDFLAGS',
        'PY_INCLUDES',
        'PY_OPTS_OMP',
        'PY_NO_OMP',
        'PY_BUILD_PARALLEL',
    ]
    for env_var in PY_ENV_VARS:
        prioprint("{}={}".format(env_var, os.environ.get(env_var)))


def display_py_options():
    """Display customised compilation options used in build.

    """
    PY_OPT_TYPES = [
        'macros',
        'cflags',
        'ldflags',
        'libs',
        'include_dirs',
        'lib_dirs',
    ]
    for opt_type in PY_OPT_TYPES:
        prioprint("{}={}".format(opt_type, globals().get(opt_type)))


# ========================================================================
# Compilation
# ========================================================================

# -- Language ------------------------------------------------------------

EXT_LANG = 'c++'


# -- Environment ---------------------------------------------------------

# Default to GCC compiler and OpenMP implementation.
COMPILERS = {
    'default': 'g++',
    'linux': 'g++',
    'darwin': 'g++',  # alternatively 'clang++'
}
OPENMP_LIBS = {
    'default': 'gomp',  # alternatively 'omp'
    'linux': '',  # set by ``-fopenmp``
    'darwin': '',  # set by ``-fopenmp``
}


def get_platform():
    """Get the build platform.

    Returns
    -------
    {'linux', 'darwin', 'default'}

    """
    build_platform = platform.system().lower()
    if build_platform not in ['linux', 'darwin']:
        build_platform = 'default'

    return build_platform


def get_compiler():
    """Get the build compiler.

    Returns
    -------
    str

    """
    compiler = os.environ.get('PY_CXX', COMPILERS[get_platform()])

    return compiler


def set_cli_compiler(compiler=None):
    """Set command-line compiler through the ``CC``, ``CXX`` and ``CPP``
    environmental variables.

    """
    compiler = compiler or get_compiler()

    os.environ['CC'] = compiler
    os.environ['CXX'] = compiler
    os.environ['CPP'] = compiler


# -- Dependencies --------------------------------------------------------

LIBS_CORE = ['gsl', 'fftw3',]  # noqa: E231
LIBS_FULL = ['gsl', 'gslcblas', 'm', 'fftw3',]  # noqa: E231

PKG_LIB_NAME = 'trv'


# -- Options -------------------------------------------------------------

def parse_cli_cflags():
    """Parse command-line ``PY_CXXFLAGS`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str
        Parsed `macros` and `cflags`.

    """
    cli_cflags = os.environ.get('PY_CXXFLAGS', '').split()

    parsed_macros, parsed_cflags = [], []
    for cflag_ in cli_cflags:
        if cflag_.startswith('-D'):
            macro_ = cflag_.lstrip('-D').split('=')
            if len(macro_) == 1:
                parsed_macros.append((macro_.pop(), None))
            elif len(macro_) == 2:
                parsed_macros.append(tuple(macro_))
            else:
                raise ValueError(
                    "Invalid macro in CXXFLAGS: {}.".format(cflag_)
                )
        else:
            parsed_cflags.append(cflag_)

    return parsed_macros, parsed_cflags


def parse_cli_ldflags():
    """Parse command-line ``PY_LDFLAGS`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str, list of str
        Parsed `ldflags`, `libs` and `libdirs`.

    """
    cli_ldflags = os.environ.get('PY_LDFLAGS', '').split()

    parsed_ldflags, parsed_libs, parsed_lib_dirs = [], [], []
    for ldflag_ in cli_ldflags:
        if ldflag_.startswith('-l'):
            parsed_libs.append(ldflag_.lstrip('-l'))
        elif ldflag_.startswith('-L'):
            parsed_lib_dirs.append(ldflag_.lstrip('-L'))
        else:
            parsed_ldflags.append(ldflag_)

    return parsed_ldflags, parsed_libs, parsed_lib_dirs


def parse_cli_includes():
    """Parse command-line ``PY_INCLUDES`` components for extension
    keyword arguments.

    Returns
    -------
    list of str
        Parsed `ldflags`.

    """
    parsed_include_dirs = \
        os.environ.get('PY_INCLUDES', '').replace('-I', '').split()

    return parsed_include_dirs


def parse_cli_opts_omp():
    """Parse command-line ``PY_OPTS_OMP`` components for extension
    keyword arguments.

    """
    cli_ldflags_omp = os.environ.get('PY_OPTS_OMP')
    if cli_ldflags_omp is None:
        return None

    parsed_ldflags_omp, parsed_libs_omp, parsed_lib_dirs_omp = [], [], []
    for ldflag_ in cli_ldflags_omp:
        if ldflag_.startswith('-l'):
            parsed_libs_omp.append(ldflag_.lstrip('-l'))
        elif ldflag_.startswith('-L'):
            parsed_lib_dirs_omp.append(ldflag_.lstrip('-L'))
        else:
            parsed_ldflags_omp.append(ldflag_)

    return parsed_ldflags_omp, parsed_libs_omp, parsed_lib_dirs_omp


def add_options_nonomp(macros, cflags, ldflags, libs, lib_dirs, include_dirs):
    """Add any missing non-OpenMP compilation options.

    Parameters
    ----------
    macros : list of tuple (str, Union[str, None])
        Macros without '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs : list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple (str, Union[str, None])
        Macros without '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs :list of str
        ``INCLUDES`` directories without the '-I' prefix.

    """
    libs_core = deepcopy(LIBS_CORE)
    libs_full = deepcopy(LIBS_FULL)

    for lib_ in libs:
        if lib_ not in libs_full:
            libs_full.append(lib_)
            libs_core.append(lib_)

    checked_options = pkg_config(libs_core, default_libraries=libs_full)

    for macro_ in checked_options['define_macros']:
        if macro_ not in macros:
            macros.append(macro_)
    for macro_ in checked_options['undef_macros']:
        if macro_ not in macros:
            macros.append(macro_)
    for flag in checked_options['extra_compile_args']:
        if flag not in cflags + ldflags:
            cflags.append(flag)
    for include_dir_ in checked_options['include_dirs']:
        if include_dir_ not in include_dirs:
            include_dirs.append(include_dir_)
    for lib_dir_ in checked_options['library_dirs']:
        if lib_dir_ not in lib_dirs:
            lib_dirs.append(lib_dir_)

    libs = checked_options['libraries']

    return macros, cflags, ldflags, libs, lib_dirs, include_dirs


def add_options_openmp(macros, cflags, ldflags, libs, lib_dirs, include_dirs):
    """Add OpenMP options, if enabled.

    Parameters
    ----------
    macros : list of tuple (str, Union[str, None])
        Macros without '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs :list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple (str, Union[str, None])
        Processed macros without '-D' prefix.
    cflags : list of str
        Processed ``CXXFLAGS`` components.
    ldflags : list of str
        Processed ``LDFLAGS`` components without '-l' prefix.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs :list of str
        ``INCLUDES`` directories without the '-I' prefix.

    """
    # Check if OpenMP is explicitly disabled.
    if os.environ.get('PY_NO_OMP') is not None:
        prioprint("OpenMP is disabled explicitly.")
        return macros, cflags, ldflags, libs, lib_dirs, include_dirs

    # Check if OpenMP is unsupported.
    if not check_openmp_support():
        prioprint(
            "OpenMP is disabled as "
            "it is not supported in this build enviromnent."
        )
        return macros, cflags, ldflags, libs, lib_dirs, include_dirs

    # Otherwise, enabled OpenMP by default.
    prioprint("OpenMP is enabled by default.")

    # Adapt macros, with the removal of duplicate options from `cflags`.
    MACROS_OMP = [
        ('TRV_USE_OMP', None),
        ('TRV_USE_FFTWOMP', None),
    ]

    for macro_ in MACROS_OMP:
        if macro_ not in macros:
            macros.append(macro_)

        macro_cflag_ = '-D' + macro_[0]
        if macro_cflag_ in cflags:
            cflags.remove(macro_cflag_)

    # Adapt `cflags`.
    CXXFLAGS_OMP = [
        '-fopenmp',
    ]

    for cflag_ in CXXFLAGS_OMP:
        if cflag_ not in cflags:
            cflags.append(cflag_)

    # Adapt `ldflags`, `libs` and `lib_dirs`.
    LIBS_OMP_BASED = [
        'fftw3_omp',
    ]  # noqa: E231
    for lib_ in LIBS_OMP_BASED:
        if lib_ and lib_ not in libs:
            libs.append(lib_)

    try:
        ldflags_omp, libs_omp, lib_dirs_omp = parse_cli_opts_omp()
    except TypeError:
        ldflags_omp = ['-fopenmp',]  # noqa: E231
        libs_omp = [OPENMP_LIBS[get_platform()],]  # noqa: E231
        lib_dirs_omp = []  # noqa: E231

    for ldflag_ in ldflags_omp:
        if ldflag_ not in ldflags:
            ldflags.append(ldflag_)
    for lib_ in libs_omp:
        if lib_ and lib_ not in libs:
            libs.append(lib_)
    for lib_dir_ in lib_dirs_omp:
        if lib_dir_ not in lib_dirs:
            lib_dirs.append(lib_dir_)

    return macros, cflags, ldflags, libs, lib_dirs, include_dirs


def add_options_pkgs(macros, cflags, ldflags, libs, lib_dirs, include_dirs):
    """Add options required by this package and external packages.

    Parameters
    ----------
    macros : list of tuple (str, Union[str, None])
        Macros without '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs : list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple (str, Union[str, None])
        Processed macros without '-D' prefix.
    cflags : list of str
        Processed ``CXXFLAGS`` components.
    ldflags : list of str
        Processed ``LDFLAGS`` components without '-l' prefix.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs : list of str
        ``INCLUDES`` directories without the '-I' prefix.

    """
    # Adapt `macros`.
    MACROS_PKG = [
        ('TRV_EXTCALL', None),
        # ('TRV_USE_LEGACY_CODE', None),
        # ('DBG_MODE', None),
        # ('DBG_FLAG_NOAC', None),
    ]
    MACROS_EXT = [
        ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
    ]

    for macro_ in MACROS_PKG + MACROS_EXT:
        macros.append(macro_)

    # Adapt `libs`. Note the linking order matters.
    libs = [PKG_LIB_NAME,] + libs  # noqa: E231

    # Adapt `include_dirs`.
    include_dirs_pkg = [
        os.path.abspath(get_pkg_include_dir()),
    ]
    include_dirs_ext = [
        numpy.get_include(),
    ]

    for include_dir_ in include_dirs_pkg + include_dirs_ext:
        if include_dir_ not in include_dirs:
            include_dirs.append(include_dir_)

    return macros, cflags, ldflags, libs, lib_dirs, include_dirs


# ========================================================================
# Build
# ========================================================================


# ========================================================================
# Extensions
# ========================================================================

CYTHON_DIRECTIVES = {
   'language_level': '3',
   'c_string_encoding': 'utf-8',
   'embedsignature': True,
}

EXT_CONFIGS = {
    'parameters': {},
    'dataobjs': {},
    '_particles': {},
    '_twopt': {},
    '_threept': {},
    '_fftlog': {},
}


def define_pkg_library(lib_name=PKG_LIB_NAME):
    """Define package library to be built and linked against.

    Parameters
    ----------
    lib_name : str
        Package libaray name.

    Returns
    -------
    tuple (str, dict)
        Package library and corresponding configurations.

    """
    pkg_src_dir = get_pkg_src_dir()

    pkg_sources = [
        os.path.join(pkg_src_dir, _cpp_source)
        for _cpp_source in os.listdir(pkg_src_dir)
    ]

    pkg_cfg = {
        cfg_opt: globals()[cfg_opt]
        for cfg_opt in ['macros', 'cflags', 'include_dirs',]  # noqa: E231
    }

    pkg_library = (
        lib_name,
        dict(sources=pkg_sources, **pkg_cfg)
    )

    return pkg_library


def define_pkg_extension(ext_name,
                         auto_cpp_source=False, extra_cpp_sources=None,
                         **ext_kwargs):
    """Define package extension from given source files and options.

    Parameters
    ----------
    ext_name : str
        Extension module name in the form ``<pkg_name>.<ext_name>``.
        This sets the .pyx source file to ``<pkg_dir>/<ext_name>.pyx``.
    auto_cpp_source : bool, optional
        If `True` (default is `False`), add
        ``<pkg_src_dir>/<ext_name>.cpp`` to the list of source files
        with any prefix underscores removed from ``<ext_name>``.
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
    EXT_CFG_MAPPING = [
        ('define_macros', 'macros'),
        ('extra_compile_args', 'cflags'),
        ('extra_link_args', 'ldflags'),
        ('include_dirs', 'include_dirs'),
        ('library_dirs', 'lib_dirs'),
        ('libraries', 'libs'),
    ]

    # Define Cython sources.
    pkg_dir = get_pkg_dir()

    sources = [os.path.join(pkg_dir, f"{ext_name}.pyx"),]  # noqa: E231

    # Add C++ sources.
    pkg_src_dir = get_pkg_src_dir()

    if auto_cpp_source:
        sources.append(
            os.path.join(pkg_src_dir, f"{ext_name}.cpp".lstrip('_'))
        )

    if extra_cpp_sources is None:
        extra_cpp_sources = []
    for _cpp_source in extra_cpp_sources:
        sources.append(os.path.join(pkg_src_dir, _cpp_source))

    # Determine extension configurations.
    ext_module_kwargs = {
        cfg: globals()[cfg_var]
        for (cfg, cfg_var) in EXT_CFG_MAPPING
    }

    if ext_kwargs:
        ext_module_kwargs.update(ext_kwargs)

    return Extension(
        '{}.{}'.format(get_pkg_name(), ext_name),
        sources=sources, language=EXT_LANG, **ext_module_kwargs
    )


# ========================================================================
# Set-up
# ========================================================================

if __name__ == '__main__':

    # Check build environment.
    display_py_environs()

    # Parse compilation options.
    set_cli_compiler()

    macros, cflags = parse_cli_cflags()
    ldflags, libs, lib_dirs = parse_cli_ldflags()
    include_dirs = parse_cli_includes()

    # Adapt compilation options.
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = add_options_nonomp(
        macros, cflags, ldflags, libs, lib_dirs, include_dirs
    )
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = add_options_openmp(
        macros, cflags, ldflags, libs, lib_dirs, include_dirs
    )
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = add_options_pkgs(
        macros, cflags, ldflags, libs, lib_dirs, include_dirs
    )

    # Check compilation options.
    display_py_options()

    # Define build targets.
    pkg_libraries = [define_pkg_library(),]  # noqa: E231

    pkg_extensions = [
        define_pkg_extension(ext, **cfg)
        for ext, cfg in EXT_CONFIGS.items()
    ]

    # Run build commands.
    cython_ext_modules = cythonize(
        pkg_extensions,
        compiler_directives=CYTHON_DIRECTIVES,
        nthreads=get_build_num_procs(),
    )

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
        libraries=pkg_libraries,
    )
