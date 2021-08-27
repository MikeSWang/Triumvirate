import distutils.sysconfig
import os
from distutils.core import setup, Extension

import numpy
from Cython.Build import cythonize
from Cython.Distutils import build_ext

os.environ['CC'] = 'g++'

config_vars = distutils.sysconfig.get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')

hankel_fftlog_ext = Extension(
    'hankel_fftlog',
    sources=["triumvirate/modelling/hankel_fftlog.pyx"],
    language='c++',
    extra_compile_args=['-std=c++11'],
    include_dirs=[".", numpy.get_include()],
    library_dirs=["."],
    libraries=['m', 'gsl'],
    define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
)

setup(
    name='hankel_fftlog',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(hankel_fftlog_ext, language_level='3')
)
