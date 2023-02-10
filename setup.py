"""Set up Triumvirate and its Cythonised extension modules.

"""
import os
import platform
from setuptools import setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy


PKG_NAME = 'Triumvirate'


# -- Repository ----------------------------------------------------------

pkg_dir = pkg_name = PKG_NAME.lower()

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
    lib for lib in os.environ.get('PY_LDFLAGS', '').split()
    if not lib.startswith('-l')
]

py_trvomp = os.environ.get('PY_TRVOMP')  # reserved env. var. for ``setup.py``
if py_trvomp is not None and py_trvomp.lower() in ['', '1', 'true']:
    for flag in ['-fopenmp', '-DTRV_USE_OMP', '-DTRV_USE_FFTWOMP']:
        if flag not in cflags:
            cflags.extend(flag)


# -- Extensions ----------------------------------------------------------

class BuildExt(build_ext):
    """Modify :class:`Cython.Distutils.build_ext`.

    """
    def build_extensions(self):
        """Override the compiler and compilation options.

        """
        # self.compiler.set_executable('compiler_cxx', compiler)
        # self.compiler.set_executable('compiler_so', compiler)
        # self.compiler.set_executable('linker_so', compiler)

        if '-Wstrict-prototypes' in self.compiler.compiler_so:
            self.compiler.compiler_so.remove('-Wstrict-prototypes')

        super().build_extensions()


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
    incl for incl
    in os.environ.get('PY_INCLUDES', "").replace("-I", "").split()
    if pkg_dir not in incl
]

ext_libraries = ['gsl', 'gslcblas', 'm', 'fftw3', 'fftw3_omp',]  # noqa: E231

includes = [pkg_include, npy_include,] + ext_includes  # noqa: E231
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
    'parameters': {'extra_cpp_sources': ["monitor.cpp",]},  # noqa: E231
    'dataobjs': {
        'extra_cpp_sources': ["monitor.cpp", "parameters.cpp",]  # noqa: E231
    },
    '_particles': {'extra_cpp_sources': ["monitor.cpp",]},  # noqa: E231
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
        'extra_cpp_sources': [
            "monitor.cpp", "maths.cpp", "arrayops.cpp",  # noqa: E231
        ],
        'libraries': libraries,
    },
}

cython_directives = {
   'language_level': '3',
   'c_string_encoding': 'utf-8',
   'embedsignature': True,
}

cython_modules = [
    return_extension(key, **val)
    for key, val in module_config.items()
]


# -- Set-up --------------------------------------------------------------

if __name__ == '__main__':
    setup(
        # license=pkg_info.get('__license__'),
        cmdclass={'build_ext': BuildExt},
        ext_modules=cythonize(
            cython_modules, compiler_directives=cython_directives
        )
    )
