<p align="center">
<img src="https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/Triumvirate.png" alt=Triumvirate-Logo width=67%>
</p>

# Three-Point Clustering Measurements in LSS

[![Build](
https://github.com/MikeSWang/Triumvirate/actions/workflows/build.yml/badge.svg?branch=main
)](https://github.com/MikeSWang/Triumvirate/actions/workflows/build.yml)
[![Tests](
https://github.com/MikeSWang/Triumvirate/actions/workflows/build.yml/badge.svg?branch=main
)](https://github.com/MikeSWang/Triumvirate/actions/workflows/tests.yml)
[![Docs](
https://readthedocs.org/projects/triumvirate/badge/?version=stable
)](https://triumvirate.readthedocs.io/en/stable)
[![Release](
https://img.shields.io/github/v/release/MikeSWang/Triumvirate?label=release
)](https://github.com/MikeSWang/Triumvirate/releases/latest)
[![Licence](
https://img.shields.io/badge/licence-GPL--3.0-informational
)](https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE)

``Triumvirate`` is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS) cosmological
analyses.


## Documentation

Comprehensive documentation including installation instructions, tutorials
and API can be found at [triumvirate.readthedocs.io](
https://triumvirate.readthedocs.io).

The scientific context is explained on the [Background](
https://triumvirate.readthedocs.io/en/stable/background.html) page
in the documentation.


## Installation

### Python package

``Triumvirate`` as a Python package is distributed through [PyPI](
https://pypi.org/project/Triumvirate) and [conda-forge](
https://anaconda.org/conda-forge/triumvirate). Instructions for installation
can be found on the [Installation](
https://triumvirate.readthedocs.io/en/stable/installation.html#python-package)
page in the documentation.

### C++ program

``Triumvirate`` as a 'black-box' C++ program is built as a binary executable
using `make`. Instructions for compilation can be found on the [Installation](
https://triumvirate.readthedocs.io/en/stable/installation.html#c-program)
page in the documentation.

> **NOTE**: In the future, the C++ code will also be released as a library.

### Development mode

Both the Python package and the C++ program can be set up in development
mode with `make`. First `git clone` this repository and `git checkout` the
branch/release you would like to edit
```bash
git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
cd Triumvirate
```
then at the repository directory root run
```bash
make clean
make [py|cpp]install [useomp=[true|1]]
```
where ``install`` builds both and ``pyinstall``/``cppinstall`` is for
Python/C++ build only. To enable OpenMP parallelisation, append
``useomp=true`` or ``useomp=1`` to the end of the second line as shown above.

The latest release is on the ``main`` branch. The default ``Makefile``
(located at the repository diretory root) suits most use cases, but you may
modify it as appropriate for your need.


## Attribution

[![JOSS](
https://img.shields.io/badge/JOSS-doi-brightgreen
)](https://joss.theoj.org/papers/?/status.svg)
[![arXiv](
https://img.shields.io/badge/arXiv-yymm.%3F-important
)](https://arxiv.org/abs/?.?)
[![MNRAS](
https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fsty3249-blue
)](https://doi.org/10.1093/mnras/sty3249)
[![MNRAS](
https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fstx2333-blue
)](https://doi.org/10.1093/mnras/stx2333)

To acknowledge the use of ``Triumvirate`` in your published research, please
cite the publications linked above which contain the relevant information
in the BibTeX format.


## Acknowledgement

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository [``hitomi``](
https://github.com/naonori/hitomi).


## Contributing

User feedback and contributions are highly valued and very welcome. Please
refer to the contribution guidelines in [``CONTRIBUTING.md``](
https://github.com/MikeSWang/Triumvirate/blob/main/CONTRIBUTING.md).


## Releases

Changes in recent releases are listed in the change log in
[``CHANGELOG.md``](
https://github.com/MikeSWang/Triumvirate/blob/main/CHANGELOG.md).


## Licence

``Triumvirate`` is made freely available under the [GPL-3.0 licence](
https://www.gnu.org/licenses/gpl-3.0.en.html). Please see [``LICENCE``](
https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE) for full
terms and conditions.

&copy; 2023 MS Wang & NS Sugiyama
