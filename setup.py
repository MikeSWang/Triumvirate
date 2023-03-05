"""Set up Triumvirate and its Cythonised extension modules.

"""
import os
import platform
from setuptools import setup
from setuptools.command.build_clib import build_clib

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy


PKG_NAME = 'Triumvirate'


# -- Repository ----------------------------------------------------------

pkg_name = PKG_NAME.lower()
pkg_dir = os.path.join('src', pkg_name)

# # Extract package information.
# pkg_info = {}
# with open(os.path.join(pkg_dir, "__init__.py")) as pkg_pyfile:
#     exec(pkg_pyfile.read(), pkg_info)

# # Determine repository branch.
# version = pkg_info.get('__version__')
# if any([segment in version for segment in ['a', 'b', 'rc']]):
#     branch = 'main'
# elif 'dev' in version:
#     branch = 'dev'
# else:
#     branch = 'stable'


# -- Compilation ---------------------------------------------------------

# Specify language.
language = 'c++'

# Enforce compiler choice.
if platform.system() == 'Linux':
    compiler = 'g++'
elif platform.system() == 'Darwin':
    compiler = 'clang++'
else:
    compiler = 'g++'

os.environ['CC'] = os.environ.get('PY_CXX', compiler)
os.environ['CXX'] = os.environ.get('PY_CXX', compiler)

# Modify compilation options.
cflags = os.environ.get('PY_CFLAGS', '').split()
ldflags = [
    lib
    for lib in os.environ.get('PY_LDFLAGS', '').split()
    if not lib.startswith('-l')
]

# Add optional compilation option: OpenMP.
py_trvomp = os.environ.get('PY_TRVOMP')
if py_trvomp is not None and py_trvomp.lower() in ['', '1', 'true']:
    for flag in ['-fopenmp', '-DTRV_USE_OMP', '-DTRV_USE_FFTWOMP']:
        if flag not in cflags:
            cflags.extend(flag)


# -- Extensions ----------------------------------------------------------

class BuildExt(build_ext):
    """Modified :class:`Cython.Distutils.build_ext`.

    """
    def build_extensions(self):
        """Override the compiler and compilation options.

        """
        # _compiler = os.environ.get('PY_CXX', compiler)
        # self.compiler.set_executable('compiler_cxx', _compiler)
        # self.compiler.set_executable('compiler_so', _compiler)
        # self.compiler.set_executable('linker_so', _compiler)

        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')

        super().build_extensions()


class BuildClib(build_clib):
    """Modified :class:`setuptools.command.build_clib`.

    """
    def build_libraries(self, libraries):
        """Override the compiler and compilation options.

        """
        # _compiler = os.environ.get('PY_CXX', compiler)
        # self.compiler.set_executable('compiler_cxx', _compiler)
        # self.compiler.set_executable('compiler_so', _compiler)
        # self.compiler.set_executable('linker_so', _compiler)

        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')

        super().build_libraries(libraries)


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
        'include_dirs': includes,
        'extra_compile_args': cflags,
        'extra_link_args': ldflags,
        'define_macros': macros,
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
    for cpp_source in extra_cpp_sources:
        sources.append(os.path.join(pkg_src_dir, cpp_source))

    return Extension(
        f'{pkg_name}.{module_name}',
        sources=sources, language=language, **ext_module_kwargs
    )


# Set package macros, sources, include directories and libraries.
pkg_macros = [
    ('TRV_EXTCALL', None),
    # ('TRV_USE_LEGACY_CODE', None),
    # ('DBG_MODE', None),
    # ('DBG_FLAG_NOAC', None),
]

pkg_include_dir = os.path.join(pkg_dir, "include")

pkg_src_dir = os.path.join(pkg_dir, "src/modules")

pkg_library = (
    'trv',
    {
        'sources': [
            os.path.join(pkg_src_dir, cpp_source)
            for cpp_source in os.listdir(pkg_src_dir)
        ],
        'include_dirs': [pkg_include_dir,],  # noqa: E231
        'macros': pkg_macros,
    }
)

# Set macros.
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
ext_libs = ['gsl', 'gslcblas', 'm', 'fftw3', 'fftw3_omp',]  # noqa: E231

libs = [] + ext_libs  # noqa: E231

# Define extension modules.
ext_module_names = [
    'parameters',
    'dataobjs',
    '_particles',
    '_twopt',
    '_threept',
    '_fftlog',
]

ext_module_config = {
    ext_mod_name: {'libraries': libs}
    for ext_mod_name in ext_module_names
}

cython_directives = {
   'language_level': '3',
   'c_string_encoding': 'utf-8',
   'embedsignature': True,
}

cython_ext_modules = cythonize(
    [define_extension(name, **cfg) for name, cfg in ext_module_config.items()],
    compiler_directives=cython_directives
)


# -- Set-up --------------------------------------------------------------

if __name__ == '__main__':
    setup(
        # license=pkg_info.get('__license__'),
        cmdclass={
            'build_clib': BuildClib,
            'build_ext': BuildExt,
        },
        ext_modules=cython_ext_modules,
        libraries=[pkg_library,]  # noqa: E231
    )
