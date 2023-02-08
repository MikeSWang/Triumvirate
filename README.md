<p align="center">
<img src="https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/Triumvirate.png" alt=Triumvirate-Logo width=50%>
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

`Triumvirate` is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS) cosmological
analysis.


## Documentation

Comprehensive documentation including installation instructions, tutorials
and API can be found at [https://triumvirate.readthedocs.io](https://triumvirate.readthedocs.io).


## Installation

### Python package

Python releases are distributed through [PyPI](https://pypi.org/project/Triumvirate)
and [conda-forge](https://anaconda.org/conda-forge/triumvirate).
Instructions for installation can be found in the
[documentation](https://triumvirate.readthedocs.io/en/stable/install.html).

### C++ program

The C++ program can be custom-built using `make`. Please refer to the
[documentation](https://triumvirate.readthedocs.io/en/stable/install.html)
for detailed instructions.

### Development mode

Both the Python package and the C++ program can be built in development
mode with `make`. After `git clone`-ing this repository, run the following
commands at the root of the directory
```
make clean
make install
```
to install both, or replace `install` with `pyinstall`/`cppinstall` for
Python/C++ build only. To enable OpenMP, append `useomp=true` or `useomp=1`
to the end of the second line.

Depending on the environment, you may wish to modify [`Makefile`](
https://github.com/MikeSWang/Triumvirate/tree/main/Makefile) as appropriate.


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

To acknowledge the use of `Triumvirate` in your published research, please
cite the publications linked above which contain the relevant information
in the BibTeX format.


## Acknowledgement

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in [`hitomi`](https://github.com/naonori/hitomi).


## Licence

`Triumvirate` is made freely available under the [GPL-3.0 licence](
https://www.gnu.org/licenses/gpl-3.0.en.html). Please see [`LICENCE`](
https://github.com/MikeSWang/Triumvirate/tree/main/Makefile) for full
terms and conditions.

&copy; 2023 MS Wang & NS Sugiyama
