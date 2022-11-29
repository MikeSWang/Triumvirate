"""Set-up for the Python package and Cythonised extension modules.

"""
import os
from distutils.sysconfig import get_config_vars
from distutils.util import convert_path
from setuptools import find_packages, setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy


PKG_NAME = 'Triumvirate'

# -- Repository ----------------------------------------------------------

# Extract instructions and package dependencies.
with open("README.md", 'r') as freadme:
    readme = freadme.read()

with open("requirements.txt", 'r') as requirements_file:
    requirements = [pkg.strip() for pkg in requirements_file]

# Extract package information.
pkgdir = PKG_NAME.lower()
pkginfo = {}
with open(convert_path(os.path.join(pkgdir, "_pkginfo.py"))) as pkginfo_file:
    exec(pkginfo_file.read(), pkginfo)

# Determine repository branch.
version_tag = pkginfo.get('__version__')
if any([segment in version_tag for segment in ['a', 'b', 'rc']]):
    branch = 'main'
elif 'dev' in version_tag:
    branch = 'dev'
else:
    branch = 'stable'


# -- Compilation ---------------------------------------------------------

# Specify language.
language = 'c++'

# Set compiler.
os.environ['CC'] = 'g++'

# Modify compilation options.
options = ['-std=c++11',]
links = []

if int(os.environ.get('PY_USEOMP', 0)):
    options.append('-fopenmp')
    links.append('-fopenmp')

# Suppress irrelevant compiler warnings.
config_vars = get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')


# -- Extensions ----------------------------------------------------------

# Set source, include and library paths.
self_modulesrc = os.path.join(pkgdir, "src/modules")

self_include = os.path.join(pkgdir, "include")
npy_include = numpy.get_include()

ext_includes = os.environ.get('PY_INCLUDES', '').replace("-I", "").split()
ext_includes = [incl_ for incl_ in ext_includes if pkgdir not in incl_]

ext_libraries = ['gsl', 'gslcblas', 'fftw3', 'fftw3_omp',]

includes = [self_include,] + [npy_include,] + ext_includes
libraries = ext_libraries

# Set macros.
self_macros = [
    ('TRV_EXTCALL', None),
    # ('TRV_USE_LEGACY_CODE', None),
    # ('DBG_MODE', None),
    # ('DBG_NOAC', None),
]

if int(os.environ.get('PY_USEOMP', 0)):
    self_macros.append(('TRV_USE_OMP', None))
    self_macros.append(('TRV_USE_FFTWOMP', None))
if int(os.environ.get('PY_DBGPARS', 0)):
    self_macros.append(('DBG_MODE', None))
    self_macros.append(('DBG_PARS', None))

npy_macros = [
    ('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION'),
]

macros = self_macros + npy_macros

# Define extension modules.
modules = [
    Extension(
        f'{pkgdir}.parameters',
        sources=[
            os.path.join(pkgdir, "parameters.pyx"),
            os.path.join(self_modulesrc, "parameters.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        define_macros=macros,
    ),
    Extension(
        f'{pkgdir}.dataobjs',
        sources=[
            os.path.join(pkgdir, "dataobjs.pyx"),
            os.path.join(self_modulesrc, "dataobjs.cpp"),
            os.path.join(self_modulesrc, "parameters.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        define_macros=macros,
    ),
    Extension(
        f'{pkgdir}._particles',
        sources=[
            os.path.join(pkgdir, "_particles.pyx"),
            os.path.join(self_modulesrc, "particles.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        define_macros=macros,
    ),
    Extension(
        f'{pkgdir}._twopt',
        sources=[
            os.path.join(pkgdir, "_twopt.pyx"),
            os.path.join(self_modulesrc, "twopt.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
            os.path.join(self_modulesrc, "maths.cpp"),
            os.path.join(self_modulesrc, "parameters.cpp"),
            os.path.join(self_modulesrc, "dataobjs.cpp"),
            os.path.join(self_modulesrc, "particles.cpp"),
            os.path.join(self_modulesrc, "field.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        libraries=libraries,
        define_macros=macros,
    ),
    Extension(
        f'{pkgdir}._threept',
        sources=[
            os.path.join(pkgdir, "_threept.pyx"),
            os.path.join(self_modulesrc, "threept.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
            os.path.join(self_modulesrc, "maths.cpp"),
            os.path.join(self_modulesrc, "parameters.cpp"),
            os.path.join(self_modulesrc, "dataobjs.cpp"),
            os.path.join(self_modulesrc, "particles.cpp"),
            os.path.join(self_modulesrc, "field.cpp"),
            os.path.join(self_modulesrc, "twopt.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        libraries=libraries,
        define_macros=macros,
    ),
    Extension(
        f'{pkgdir}._fftlog',
        sources=[
            os.path.join(pkgdir, "_fftlog.pyx"),
            os.path.join(self_modulesrc, "fftlog.cpp"),
            os.path.join(self_modulesrc, "monitor.cpp"),
            os.path.join(self_modulesrc, "maths.cpp"),
            os.path.join(self_modulesrc, "arrayops.cpp"),
        ],
        language=language,
        extra_compile_args=options,
        extra_link_args=links,
        include_dirs=includes,
        libraries=libraries,
        define_macros=macros,
    ),
]


# -- Set-up --------------------------------------------------------------

setup(
    name=PKG_NAME.capitalize(),
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
    python_requires='>=3.6',
    install_requires=requirements,
    packages=find_packages(),
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(
        modules,
        language_level='3',
    )
)
