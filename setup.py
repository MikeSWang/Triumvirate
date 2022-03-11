import os
from distutils.sysconfig import get_config_vars
from setuptools import find_packages, setup

from Cython.Build import cythonize
from Cython.Distutils import build_ext, Extension

import numpy

# -- Repository ---------------------------------------------------------------

with open("version.txt", 'r') as version_info:
    version_number, version_tag = [v.strip() for v in version_info]
    if version_tag == 'latest':
        branch = 'main'
    elif version_tag == 'development':
        branch = 'dev'
    else:
        branch = version_tag

with open("README.md", 'r') as readme:
    long_description = readme.read()

with open("requirements.txt", 'r') as requirements:
    dependencies = [pkg.strip() for pkg in requirements]


# -- Compilation --------------------------------------------------------------

os.environ['CC'] = 'g++'

config_vars = get_config_vars()
for key, val in config_vars.items():
    if isinstance(val, str):
        config_vars[key] = val.replace('-Wstrict-prototypes', '')


# -- Extensions ---------------------------------------------------------------

# kwarg: define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
ext_modules = [
    Extension(
        'triumvirate.parameters',
        sources=["triumvirate/parameters.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11',]
    ),
    Extension(
        'triumvirate._catalogue',
        sources=["triumvirate/_catalogue.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11',],
        include_dirs=[numpy.get_include(),],
    ),
    Extension(
        'triumvirate._twopt',
        sources=["triumvirate/_twopt.pyx"],
        language='c++',
        extra_compile_args=['-std=c++11',],
        include_dirs=["triumvirate/include", numpy.get_include(),],
        libraries=['m', 'gsl', 'fftw3', 'gslcblas'],
    ),
]

setup(
    name="Triumvirate",
    version=version_number,
    license="GPLv3",
    author="Mike S Wang",
    author_email="mikeshengbo.wang@ed.ac.uk",
    description=(
        "Measuring three-point correlators in "
        "large-scale structure clustering analysis."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Cython",
        "Programming Language :: C++",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Information Analysis",
    ],
    url="https://github.com/MikeSWang/Triumvirate/",
    project_urls={
        "Documentation": "https://mikeswang.github.io/Triumvirate",
        "Source": "https://github.com/MikeSWang/Triumvirate/",
    },
    packages=find_packages(),
    python_requires='>=3.6',
    install_requires=dependencies,
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize(ext_modules, language_level='3')
)
