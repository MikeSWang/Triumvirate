import os
from distutils.sysconfig import get_config_vars
from distutils.util import convert_path
from setuptools import find_packages, setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy


PKG_NAME = 'Triumvirate'

# -- Repository ---------------------------------------------------------------

pkgdir = PKG_NAME.lower()
pkginfo = {}

with open(convert_path(f"{pkgdir}/_pkginfo.py")) as fpinfo:
    exec(fpinfo.read(), pkginfo)

with open("README.md", 'r') as freadme:
    readme = freadme.read()

with open("requirements.txt", 'r') as frequirements:
    requirements = [pkg.strip() for pkg in frequirements]

version_tag = pkginfo.get('_release_')
if version_tag == 'latest':
    branch = 'main'
elif version_tag == 'development':
    branch = 'dev'
else:
    branch = version_tag


# -- Compilation --------------------------------------------------------------

os.environ['CC'] = 'g++'  # change as necessary

config_vars = get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')


# -- Extensions ---------------------------------------------------------------

language = 'c++'
extra_compile_args = ['-std=c++11',]
libraries = ['m', 'gsl', 'fftw3', 'gslcblas']
self_include = f"{pkgdir}/include"
npy_include = numpy.get_include()
npy_macros = ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')

ext_modules = [
    Extension(
        f'{pkgdir}.parameters',
        sources=[f"{pkgdir}/parameters.pyx"],
        language=language,
        extra_compile_args=extra_compile_args
    ),
    Extension(
        f'{pkgdir}._catalogue',
        sources=[f"{pkgdir}/_catalogue.pyx"],
        language=language,
        extra_compile_args=extra_compile_args,
        include_dirs=[npy_include,],
        define_macros=[
            npy_macros,
            # ('DBGNZ', None),
        ],
    ),
    Extension(
        f'{pkgdir}._twopt',
        sources=[f"{pkgdir}/_twopt.pyx"],
        language=language,
        extra_compile_args=extra_compile_args,
        include_dirs=[self_include, npy_include,],
        libraries=libraries,
        define_macros=[
            npy_macros,
            # ('DBGDK', None),
        ],
    ),
    Extension(
        f'{pkgdir}._threept',
        sources=[f"{pkgdir}/_threept.pyx"],
        language=language,
        extra_compile_args=extra_compile_args,
        include_dirs=[self_include, npy_include,],
        libraries=libraries,
        define_macros=[
            npy_macros,
            # ('DBGDK', None),
        ],
    ),
    Extension(
        f'{pkgdir}._fftlog',
        sources=[f"{pkgdir}/_fftlog.pyx"],
        language=language,
        extra_compile_args=extra_compile_args,
        include_dirs=[self_include, npy_include,],
        libraries=libraries,
        define_macros=[npy_macros,],
    ),
]

setup(
    name=PKG_NAME,
    version=pkginfo.get('__version__'),
    license=pkginfo.get('__license__'),
    author=pkginfo.get('__author__'),
    author_email=pkginfo.get('__email__'),
    description=pkginfo.get('__description__'),
    long_description=readme,
    long_description_content_type="text/markdown",
    classifiers=[
        "Operating System :: OS Independent",
        "Programming Language :: C++",
        "Programming Language :: Cython",
        "Programming Language :: Python",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Information Analysis",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    project_urls={
        "Documentation": "https://mikeswang.github.io/Triumvirate",
        "Source": "https://github.com/MikeSWang/Triumvirate/",
    },
    packages=find_packages(),
    install_requires=requirements,
    python_requires='>=3.6',
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(ext_modules, language_level='3')
)
