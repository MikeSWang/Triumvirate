"""Set up Triumvirate and its Cythonised extension modules.

"""
import os
import platform
import sys
from copy import deepcopy
from distutils.log import _global_log as distutils_logger
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
        Top-level (sub-)directory containing the package at the
        repository root.

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
        Subdirectory containing C++ headers in the package directory.

    Returns
    -------
    str
        C++ header directory.

    """
    return os.path.join(get_pkg_dir(), subdir)


def get_pkg_src_dir(subdir="src"):
    """Get package C++ source directory.

    Parameters
    ----------
    subdir : str, optional
        Subdirectory containing C++ sources in the package directory.

    Returns
    -------
    str
        C++ source directory.

    """
    return os.path.join(get_pkg_dir(), subdir)


def get_pkg_version_scheme(default_ver_scheme='no-guess-dev',
                           default_loc_scheme='node-and-date'):
    """Get package version scheme from the environment.

    Parameters
    ----------
    default_ver_scheme : str, optional
        Fallback default version scheme for the package
        (default is 'no-guess-dev').
    default_loc_scheme : str, optional
        Fallback default local scheme for the package
        (default is 'node-and-date').

    Returns
    -------
    dict
        Package version scheme(s).

    See Also
    --------
    :pkg:`setuptools_scm`
        For available version schemes.

    """
    ver_scheme = os.environ.get('PY_SCM_VER_SCHEME', '').strip() \
        or default_ver_scheme
    loc_scheme = os.environ.get('PY_SCM_LOC_SCHEME', '').strip() \
        or default_loc_scheme

    scheme = {
        'version_scheme': ver_scheme,
        'local_scheme': loc_scheme,
    }

    prioprint("Versioning scheme: {}".format(scheme))

    return scheme


# ========================================================================
# Build
# ========================================================================

class BuildExt(build_ext):
    """Modified :class:`Cython.Distutils.build_ext`.

    """

    def finalize_options(self):
        """Modify parallelisation option in
        :meth:`Cython.Distutils.build_ext.finalize_options`.

        """
        super().finalize_options()

        _num_procs = get_build_num_procs()
        if self.parallel is None and _num_procs is not None:
            self.parallel = _num_procs

    def build_extensions(self):
        """Modify compiler and compilation configuration in
        :meth:`Cython.Distutils.build_ext.build_extensions`.

        """
        OPTS_TO_REMOVE = ['-Wstrict-prototypes', '-Wl,-pie',]  # noqa: E231
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
        OPTS_TO_REMOVE = ['-Wstrict-prototypes', '-Wl,-pie',]  # noqa: E231
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
    """Get the number of build parallel processes.

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
    """`print` with high priority.

    """
    print(*args, file=sys.stderr, **kwargs)


def display_py_environs():
    """Display customised environment variables passed to setup.

    """
    PY_ENV_VARS = [
        'PY_CXX',
        'PY_INCLUDES',  # typically included in 'CPPFLAGS'
        'PY_CXXFLAGS',  # untypically includes macros in 'CPPFLAGS'
        'PY_LDFLAGS',  # untypically includes 'LDLIBS'
        'PY_NO_OMP',  # disable OpenMP explicitly
        'PY_OMP',  # enable OpenMP explicitly unless overriden by `PY_NO_OMP`
        'PY_CXXFLAGS_OMP',
        'PY_LDFLAGS_OMP',
        'PY_BUILD_PARALLEL',
        'PY_SCM_VER_SCHEME',
        'PY_SCM_LOC_SCHEME',
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

EXT_LANG = 'c++'


# -- OS ------------------------------------------------------------------

def get_platform():
    """Get the build platform.

    Returns
    -------
    {'linux', 'darwin', 'default'}
        Build platform.

    """
    build_platform = platform.system().lower()
    if build_platform not in ['linux', 'darwin']:
        build_platform = 'default'

    return build_platform


# -- Compiler ------------------------------------------------------------

# List default compilers (GCC assumed).
COMPILERS = {
    'default': 'g++',
    'linux': 'g++',
    'darwin': 'g++',
}


def get_compiler():
    """Get the build compiler.

    Returns
    -------
    str
        Build compiler.

    """
    compiler = os.environ.get('PY_CXX', COMPILERS[get_platform()])

    return compiler


def set_cli_compiler(compiler=None):
    """Set command-line compiler through the ``CC``and ``CXX``
    environmental variables.

    """
    _compiler = compiler or get_compiler()

    os.environ['CC'] = _compiler
    os.environ['CXX'] = _compiler


# -- Dependencies --------------------------------------------------------

PKG_LIB_NAME = 'trv'

LIBS_CORE = ['gsl', 'fftw3',]  # noqa: E231
LIBS_FULL = ['gsl', 'gslcblas', 'm', 'fftw3',]  # noqa: E231

# Default to GCC OpenMP implementation.
OPENMP_LIBS = {
    'default': 'gomp',  # or LLVM 'omp'
    'linux': '',  # set by ``-fopenmp`` or environment
    'darwin': '',  # set by ``-fopenmp`` or environment
}


# -- Options -------------------------------------------------------------

# List default compilation options.
CFLAGS = ['-std=c++17',]


def convert_macro(macro):
    """Convert a macro flag to a tuple.

    Parameters
    ----------
    macro : str
        Macro as a flag with '-D' prefix.

    Returns
    -------
    tuple[str, Union[str, None]]
        Macro as a 2-tuple.

    """
    macro_ = macro.lstrip('-D').split('=')
    if len(macro_) == 1:
        return (macro_.pop(), None)
    if len(macro_) == 2:
        return tuple(macro_)
    raise ValueError("Invalid macro to convert: {}.".format(macro_))


def parse_cli_cflags():
    """Parse command-line ``PY_CXXFLAGS`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str
        Parsed `macros` and `cflags`.

    """
    cli_cflags = os.environ.get('PY_CXXFLAGS', '').split()

    parsed_macros, parsed_cflags = [], CFLAGS.copy()
    for cflag_ in cli_cflags:
        if cflag_.startswith('-D'):
            macro_ = convert_macro(cflag_)
            parsed_macros.append(macro_)
        else:
            parsed_cflags.append(cflag_)

    return parsed_macros, parsed_cflags


def parse_cli_ldflags():
    """Parse command-line ``PY_LDFLAGS`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str, list of str
        Parsed `ldflags`, `libs` and `lib_dirs`.

    """
    cli_ldflags = os.environ.get('PY_LDFLAGS', '').split()

    parsed_ldflags, parsed_libs, parsed_lib_dirs = [], [], []
    for ldflag_ in cli_ldflags:
        if ldflag_.startswith('-l'):
            lib_ = ldflag_.lstrip('-l')
            parsed_libs.append(lib_)
        elif ldflag_.startswith('-L'):
            ldflag_L = ldflag_.lstrip('-L')
            parsed_lib_dirs.append(ldflag_L)
        else:
            parsed_ldflags.append(ldflag_)

    return parsed_ldflags, parsed_libs, parsed_lib_dirs


def parse_cli_includes():
    """Parse command-line ``PY_INCLUDES`` components for extension
    keyword arguments.

    Returns
    -------
    list of str
        Parsed `include_dirs`.

    """
    parsed_include_dirs = \
        os.environ.get('PY_INCLUDES', '').replace('-I', '').split()

    return parsed_include_dirs


def parse_cli_cflags_omp():
    """Parse command-line ``PY_CXXFLAGS_OMP`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str
        Parsed `macros` and `include_dirs` for OpenMP.

    """
    cli_cflags_omp = os.environ.get('PY_CXXFLAGS_OMP', '').split()
    if not cli_cflags_omp:
        return None

    parsed_macros_omp, parsed_cflags_omp = [], []
    for cflag_ in cli_cflags_omp:
        if cflag_.startswith('-D'):
            macro_ = cflag_.lstrip('-D').split('=')
            if len(macro_) == 1:
                parsed_macros_omp.append((macro_.pop(), None))
            elif len(macro_) == 2:
                parsed_macros_omp.append(tuple(macro_))
            else:
                raise ValueError(
                    "Invalid macro in CXXFLAG_OMP: {}.".format(cflag_)
                )
        else:
            parsed_cflags_omp.append(cflag_)

    return parsed_macros_omp, parsed_cflags_omp


def parse_cli_ldflags_omp():
    """Parse command-line ``PY_LDFLAGS_OMP`` components for extension
    keyword arguments.

    Returns
    -------
    list of str, list of str
        Parsed `ldflags`, `libs` and `lib_dirs` for OpenMP.

    """
    cli_ldflags_omp = os.environ.get('PY_LDFLAGS_OMP', '').split()
    if not cli_ldflags_omp:
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
    """Add non-OpenMP compilation options.

    Parameters
    ----------
    macros : list of tuple[str, Union[str, None]]
        Macros without the '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components with non-'-D' and non-'-I' prefixes.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' and non-'-L' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs : list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple[str, Union[str, None]]
        Extended macros without the '-D' prefix.
    cflags : list of str
        Extended ``CXXFLAGS`` components with non-'-D' and
        non-'-I' prefixes.
    ldflags : list of str
        Extended ``LDFLAGS`` components with non-'-l' and
        non-'-L' prefixes.
    libs : list of str
        Extended libraries without the '-l' prefix.
    lib_dirs : list of str
        Extended library directories without the '-L' prefix.
    include_dirs :list of str
        Extended ``INCLUDES`` directories without the '-I' prefix.

    """
    libs_core = deepcopy(LIBS_CORE)
    libs_full = deepcopy(LIBS_FULL)

    for lib_ in libs:
        if lib_ not in libs_full:
            libs_core.append(lib_)
            libs_full.append(lib_)

    checked_options = pkg_config(libs_core, default_libraries=libs_full)

    libs = checked_options['libraries']

    for macro_ in checked_options['define_macros']:
        if macro_ not in macros:
            macros.append(macro_)
    for flag in checked_options['extra_compile_args']:
        if flag not in cflags + ldflags:
            cflags.append(flag)
    for lib_dir_ in checked_options['library_dirs']:
        if lib_dir_ not in lib_dirs:
            lib_dirs.append(lib_dir_)
    for include_dir_ in checked_options['include_dirs']:
        if include_dir_ not in include_dirs:
            include_dirs.append(include_dir_)

    return macros, cflags, ldflags, libs, lib_dirs, include_dirs


def add_options_omp(macros, cflags, ldflags, libs, lib_dirs, include_dirs):
    """Add OpenMP options, if enabled.

    By default, GCC OpenMP is assumed.

    Parameters
    ----------
    macros : list of tuple[str, Union[str, None]]
        Macros without the '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components with non-'-D' and non-'-I' prefixes.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' and non-'-L' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs :list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple[str, Union[str, None]]
        Extended macros without the '-D' prefix.
    cflags : list of str
        Extended ``CXXFLAGS`` components with non-'-D' and
        non-'-I' prefixes.
    ldflags : list of str
        Extended ``LDFLAGS`` components with non-'-l' and
        non-'-L' prefixes.
    libs : list of str
        Extended libraries without the '-l' prefix.
    lib_dirs : list of str
        Extended library directories without the '-L' prefix.
    include_dirs :list of str
        Extended ``INCLUDES`` directories without the '-I' prefix.

    """
    # Check if OpenMP is explicitly disabled.
    if os.environ.get('PY_NO_OMP') is not None:
        prioprint("OpenMP is disabled explicitly.")
        return macros, cflags, ldflags, libs, lib_dirs, include_dirs

    # Check if OpenMP is unsupported.
    if os.environ.get('PY_OMP') is None and not check_openmp_support():
        prioprint(
            "OpenMP is disabled as "
            "it is not supported in this build enviromnent."
        )
        return macros, cflags, ldflags, libs, lib_dirs, include_dirs

    # Enable OpenMP by default otherwise.
    prioprint("OpenMP is enabled.")

    # Adapt `macros` and `cflags`.
    MACROS_OMP_RELATED = [
        ('TRV_USE_OMP', None),
        ('TRV_USE_FFTWOMP', None),
    ]
    for macro_ in MACROS_OMP_RELATED:
        if macro_ not in macros:
            macros.append(macro_)

    try:
        macros_omp, cflags_omp = parse_cli_cflags_omp()
    except TypeError:
        macros_omp = []
        cflags_omp = ['-fopenmp',]  # noqa: E231

    for macro_ in macros_omp:
        macro_ = convert_macro(macro_)
        if macro_ not in macros:
            macros.append(macro_)
    for cflag_ in cflags_omp:
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
        ldflags_omp, libs_omp, lib_dirs_omp = parse_cli_ldflags_omp()
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
    macros : list of tuple[str, Union[str, None]]
        Macros without the '-D' prefix.
    cflags : list of str
        ``CXXFLAGS`` components with non-'-D' and non-'-I' prefixes.
    ldflags : list of str
        ``LDFLAGS`` components with non-'-l' and non-'-L' prefixes.
    libs : list of str
        Libraries without the '-l' prefix.
    lib_dirs : list of str
        Library directories without the '-L' prefix.
    include_dirs : list of str
        ``INCLUDES`` directories without the '-I' prefix.

    Returns
    -------
    macros : list of tuple[str, Union[str, None]]
        Extended macros without the '-D' prefix.
    cflags : list of str
        Extended ``CXXFLAGS`` components with non-'-D' and
        non-'-I' prefixes.
    ldflags : list of str
        Extended ``LDFLAGS`` components with non-'-l' and
        non-'-L' prefixes.
    libs : list of str
        Extended libraries without the '-l' prefix.
    lib_dirs : list of str
        Extended library directories without the '-L' prefix.
    include_dirs : list of str
        Extended ``INCLUDES`` directories without the '-I' prefix.

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

    # Adapt `ldflags`. Note the `rpath` option is added.
    for lib_dir_ in lib_dirs:
        ldflags.append(f'-Wl,-rpath,{lib_dir_}')

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


def cleanup_options(*args):
    """Clean up compilation options by removing duplicates.

    Parameters
    ----------
    args : list of list of str
        List of compilation options.

    Returns
    -------
    list of list of str
        Cleaned-up compilation options.

    """
    ret_args = []
    for arg in args:
        ret_args.append(list(dict.fromkeys(arg)))

    return ret_args


# ========================================================================
# Targets
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
    '_field': {},
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
    tuple[str, dict]
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

    pkg_library = (lib_name, dict(sources=pkg_sources, **pkg_cfg))

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
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = add_options_omp(
        macros, cflags, ldflags, libs, lib_dirs, include_dirs
    )
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = add_options_pkgs(
        macros, cflags, ldflags, libs, lib_dirs, include_dirs
    )
    macros, cflags, ldflags, libs, lib_dirs, include_dirs = cleanup_options(
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

    setup_commands = {
        'build_clib': BuildClib,
        'build_ext': BuildExt,
    }

    setup(
        use_scm_version=get_pkg_version_scheme(),
        cmdclass=setup_commands,
        ext_modules=cython_ext_modules,
        libraries=pkg_libraries,
    )
