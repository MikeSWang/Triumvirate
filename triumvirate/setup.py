import distutils.sysconfig
import os
from distutils.core import setup

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

os.environ['CC'] = 'g++'

config_vars = distutils.sysconfig.get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')

ext_modules = [
    Extension(
        'paramset',
        sources=["paramset.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11'],
        include_dirs=[".", numpy.get_include()],
        library_dirs=["."],
        libraries=['m', 'gsl', 'fftw3', 'gslcblas'],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    ),
    Extension(
        'params',
        sources=["params.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11'],
        include_dirs=[".", numpy.get_include()],
        library_dirs=["."],
        libraries=['m', 'gsl', 'fftw3', 'gslcblas'],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    ),
    Extension(
        'threept',
        sources=["threept.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11'],
        include_dirs=[".", numpy.get_include()],
        library_dirs=["."],
        libraries=['m', 'gsl', 'fftw3', 'gslcblas'],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    ),
    Extension(
        '_morph',
        sources=["_morph.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11'],
        include_dirs=[".", numpy.get_include()],
        library_dirs=["."],
        libraries=['m', 'gsl', 'fftw3', 'gslcblas'],
        define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
    ),
]

setup(
    name='triumvirate',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(ext_modules, language_level='3')
)
