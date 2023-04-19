<p align="center">
<img src="docs/source/_static/Triumvirate.png" alt=Triumvirate-Logo width=67%>
</p>

# Three-Point Clustering Measurements in LSS

[![CI](https://img.shields.io/github/actions/workflow/status/MikeSWang/Triumvirate/ci.yml?label=ci&logo=GitHubActions)](https://github.com/MikeSWang/Triumvirate/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/readthedocs/triumvirate/latest?logo=ReadtheDocs)](https://triumvirate.readthedocs.io/en/latest/)
[![Release](https://img.shields.io/github/v/release/MikeSWang/Triumvirate?display_name=tag&sort=semver&logo=Git)](https://github.com/MikeSWang/Triumvirate/releases/latest)

``Triumvirate`` is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS) cosmological
analyses.


## Documentation

Comprehensive documentation including the [scientific background](
https://triumvirate.readthedocs.io/en/latest/background.html),
[installation instructions](
https://triumvirate.readthedocs.io/en/latest/installation.html),
[tutorials](https://triumvirate.readthedocs.io/en/latest/tutorials.html) and
[API reference](https://triumvirate.readthedocs.io/en/latest/apiref.html)
can be found at [triumvirate.readthedocs.io](
https://triumvirate.readthedocs.io/en/latest/).


## Installation

### Python package

[![PyPI](https://img.shields.io/pypi/v/Triumvirate?logo=PyPI&color=informational)](https://pypi.org/project/Triumvirate)
[![conda](https://img.shields.io/conda/vn/msw/triumvirate?logo=Anaconda&color=informational)](https://anaconda.org/msw/triumvirate)

``Triumvirate`` as a Python package is distributed through [PyPI](
https://pypi.org/project/Triumvirate) and [conda](
https://anaconda.org/msw/triumvirate). Instructions for installation
can be found on the [Installation](
https://triumvirate.readthedocs.io/en/latest/installation.html#python-package)
page in the documentation.

### C++ library & program

``Triumvirate`` as either a static library or a binary executable can be
built using `make`. Instructions for compilation can be found on the
[Installation](
https://triumvirate.readthedocs.io/en/latest/installation.html#c-library-program)
page in the documentation.

### Development mode

Both the Python package and the C++ library/program can be set up in
development mode with `make`, provided that dependency requirements are
satisfied (GSL and FFTW3 libraries are mandatory while an OpenMP library
is optional).

First `git clone` the desired branch/release from the GitHub repository and
change into the repository directory path:

```console
$ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
$ cd Triumvirate
```

Then, execute in terminal:

```console
$ make clean
$ make ([py|cpp]install)|(cpp[libinstall|appbuild]) [useomp=(true|1)]
```

where ``cpplibinstall`` or ``cppappbuild`` respectively builds the C++
static library or binary executable only, ``cppinstall`` builds both,
``pyinstall`` builds the Python package only, and ``install`` builds
all of the above. To enable OpenMP parallelisation, append ``useomp=true`` or ``useomp=1`` to the end of the second line as shown above.

> :bulb: The latest release is on the [``main``](
> https://github.com/MikeSWang/Triumvirate/tree/main) branch. The default
> [``Makefile``](Makefile) (located at the repository directory root) should
> work in most build environments, but may need to be modified as appropriate.

> :bulb: See the [Installation](
> https://triumvirate.readthedocs.io/en/latest/installation.html#dependencies)
> page in the documentation for more details about dependency requirements.

> :warning: If enabling OpenMP, ensure the C++ compiler used supports it and is
> configured accordingly. The default [``Makefile``](Makefile) (located at the
> repository directory root) assumes the GCC compiler and OpenMP library. See
> the [Installation](
> https://triumvirate.readthedocs.io/en/latest/installation.html#openmp-support)
> page in the documentation for more details.

> :bulb: Pass option ``-j[N] -O`` to `make` to run multiple concurrent jobs
> for parallel building (optional parameter ``N`` is the number of
> parallel jobs; see [GNU Make Manual](
> https://www.gnu.org/software/make/manual/html_node/Options-Summary.html)).


## Attribution

[![JOSS](https://joss.theoj.org/papers/a8325e3897dd726d9df42286bf72d19f/status.svg)](https://joss.theoj.org/papers/a8325e3897dd726d9df42286bf72d19f)
[![arXiv](https://img.shields.io/badge/arXiv-2304.03643-b31b1b)](https://arxiv.org/abs/2304.03643)
[![MNRAS](https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fsty3249-informational)](https://doi.org/10.1093/mnras/sty3249)
[![MNRAS](https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fstx2333-informational)](https://doi.org/10.1093/mnras/stx2333)

To acknowledge the use of ``Triumvirate`` in your published research, please
cite the publications linked above which contain the relevant information
in the BibTeX format.


## Acknowledgement

<img src="docs/source/_static/ERC-Logo-Flag.png" alt="ERC" width="40%">

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(Grant agreement No. [853291](https://doi.org/10.3030/853291)).

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository [``hitomi``](
https://github.com/naonori/hitomi).


## Contributing

User feedback and contributions are very welcome. Please refer to the
[contribution guidelines](CONTRIBUTING.md).


## Discussions

A [community forum](https://github.com/MikeSWang/Triumvirate/discussions)
for users and developers exists, where you can receive
announcements, post questions, share ideas and get updates.


## Releases

Release notes are included in the [change log](CHANGELOG.md).


## Licence

[![Licence](https://img.shields.io/github/license/MikeSWang/Triumvirate?label=licence&style=flat-square&color=informational)](https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE)

``Triumvirate`` is made freely available under the [GPL-3.0 licence](
https://www.gnu.org/licenses/gpl-3.0.en.html). Please see
[``LICENCE``](LICENCE) (located at the repository directory root) for
full terms and conditions.

&copy; 2023 Mike S Wang & Naonori S Sugiyama
