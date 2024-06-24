..
    This read-me .rst file is for PyPI project description only, and
    should be periodically compared against the official read-me .md file
    which is rendered on GitHub and included in documentation.

.. figure:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/image/Triumvirate.png
    :class: dark-light
    :align: center
    :width: 67 %

==========================================
Three-Point Clustering Measurements in LSS
==========================================

.. .. image:: https://img.shields.io/github/actions/workflow/status/MikeSWang/Triumvirate/ci.yml?label=ci&logo=GitHubActions
..     :target: https://github.com/MikeSWang/Triumvirate/actions/workflows/ci.yml
..     :alt: CI

.. .. image:: https://img.shields.io/readthedocs/triumvirate/latest?logo=ReadtheDocs
..     :target: https://readthedocs.org/projects/triumvirate/builds/
..     :alt: Docs

.. image:: https://img.shields.io/github/v/release/MikeSWang/Triumvirate?display_name=tag&sort=semver&logo=Git
    :target: https://github.com/MikeSWang/Triumvirate/releases/latest
    :alt: Release

|Triumvirate| is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS) cosmological
analyses.


Documentation
=============

.. image:: https://img.shields.io/badge/Read%20the%20Docs-latest-informational?logo=ReadtheDocs
    :target: https://triumvirate.readthedocs.io/en/latest/
    :alt: Documentation

Comprehensive documentation including the `scientific background
<https://triumvirate.readthedocs.io/en/latest/background.html>`_,
`installation instructions
<https://triumvirate.readthedocs.io/en/latest/installation.html>`_,
`tutorials
<https://triumvirate.readthedocs.io/en/latest/tutorials.html>`_ and
`API reference
<https://triumvirate.readthedocs.io/en/latest/apiref.html>`_
can be found at `triumvirate.readthedocs.io
<https://triumvirate.readthedocs.io/en/latest/>`_.


Installation
============

Python package
--------------

.. image:: https://img.shields.io/pypi/v/Triumvirate?logo=PyPI&color=informational
    :target: https://pypi.org/project/Triumvirate
    :alt: PyPI

.. image:: https://img.shields.io/conda/v/msw/triumvirate?logo=Anaconda&color=informational
    :target: https://anaconda.org/msw/triumvirate
    :alt: Conda

|Triumvirate| as a Python package is distributed through
`PyPI <https://pypi.org/project/Triumvirate>`_ and
`Conda <https://anaconda.org/msw/triumvirate>`_. Instructions for installation
can be found on the `Installation
<https://triumvirate.readthedocs.io/en/latest/installation.html#python-package>`__
page in the documentation.


C++ library & program
---------------------

|Triumvirate| as either a static library or a binary executable can be
built using `make`. Instructions for compilation can be found on the
`Installation
<https://triumvirate.readthedocs.io/en/latest/installation.html#c-program>`__
page in the documentation.


Development mode
----------------

Both the Python package and the C++ library/program can be set up in
development mode with `make`, provided that dependency requirements are
satisfied (GSL and FFTW3 libraries are mandatory while an OpenMP library
is optional).

First `git clone` the desired branch/release from the GitHub repository and
change into the repository directory path:

.. code-block:: console

    $ git clone git@github.com:MikeSWang/Triumvirate.git --branch <branch-or-release>
    $ cd Triumvirate

Then, execute in shell:

.. code-block:: console

    $ make clean
    $ make [py|cpp]install [useomp=(true|1)]

where ``cpplibinstall`` or ``cppappbuild`` respectively builds the C++
static library or binary executable only, ``cppinstall`` builds both,
``pyinstall`` builds the Python package only, and ``install`` builds
all of the above. To enable OpenMP parallelisation, append ``useomp=true``
or ``useomp=1`` to the end of the second line as shown above.

.. note::

    The latest release is on the |main|_ branch. The default |Makefile|_
    (located at the repository directory root) should work in most
    build environments, but may need to be modified as appropriate.

.. note::

    See the `Installation
    <https://triumvirate.readthedocs.io/en/latest/installation.html#dependencies>`__
    page in the documentation for more details about dependency requirements.

.. note::

    If enabling OpenMP, ensure the C++ compiler used supports it and is
    configured accordingly. The default |Makefile|_ (located at the repository
    directory root) assumes the GCC compiler and OpenMP library. See the
    `Installation
    <https://triumvirate.readthedocs.io/en/latest/installation.html#openmp-support>`__
    page in the documentation for more details.

.. note::

    Pass option ``-j[N] -O`` to `make` to run multiple concurrent jobs
    for parallel building (optional parameter ``N`` is the number of
    parallel jobs; see `GNU Make Manual
    <https://www.gnu.org/software/make/manual/html_node/Options-Summary.html>`_).


Attribution
===========

.. image:: https://joss.theoj.org/papers/10.21105/joss.05571/status.svg
    :target: https://doi.org/10.21105/joss.05571
    :alt: JOSS

.. image:: https://img.shields.io/badge/zenodo-10.5281%2Fzenodo.10072128-1682D4
    :target: https://doi.org/10.5281/zenodo.10072128
    :alt: Zenodo

.. image:: https://img.shields.io/badge/arXiv-2304.03643-b31b1b
    :target: https://arxiv.org/abs/2304.03643
    :alt: arXiv

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fsty3249-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/sty3249
    :alt: MNRAS

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fstx2333-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/stx2333
    :alt: MNRAS

To acknowledge the use of |Triumvirate| in your published research, please
cite the publications linked above; for convenience, you can refer to the
files |CitationCFF|_ and |CitationMD|_ for the relevant information in
different formats.


Acknowledgement
===============

.. figure:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/image/ERC-Logo-Flag.png
    :alt: ERC
    :align: left
    :width: 40%

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(Grant agreement ID: `853291 <https://doi.org/10.3030/853291>`_).

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository |hitomi|_.

We thank the JOSS reviewers, William Coulton
(`@wcoulton <https://github.com/wcoulton>`_) and Alfonso Veropalumbo
(`@alfonso-veropalumbo <https://github.com/alfonso-veropalumbo>`_), for
their valuable feedback and suggestions (`openjournals/joss-reviews#5571
<https://github.com/openjournals/joss-reviews/issues/5571>`_),
which have improved the functionality and documentation of the code.


Contributing
============

.. image:: https://img.shields.io/conda/pn/msw/triumvirate
    :alt: Platforms

.. image:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/image/GitHub-Codespace-badge.png
    :target: https://codespaces.new/MikeSWang/Triumvirate?hide_repo_select=true&ref=main
    :alt: Codespaces
    :width: 249px

User feedback and contributions are very welcome. Please refer to the
`contribution guidelines
<https://github.com/MikeSWang/Triumvirate/blob/main/.github/CONTRIBUTING.md>`_.


Discussions & Wiki
==================

.. image:: https://img.shields.io/github/discussions/MikeSWang/Triumvirate
    :target: https://github.com/MikeSWang/Triumvirate/discussions
    :alt: Discussions

A `community forum <https://github.com/MikeSWang/Triumvirate/discussions>`_
for users and developers is hosted on GitHub, where you can receive
announcements, post questions, share ideas and get updates.

A `wiki site <https://github.com/MikeSWang/Triumvirate/wiki>`_ collects wisdoms
for specific use cases and user environments.


Releases
========

Release notes are included in the `change log
<https://github.com/MikeSWang/Triumvirate/blob/main/CHANGELOG.md>`_.


Licence
=======

.. image:: https://img.shields.io/github/license/MikeSWang/Triumvirate?label=licence&style=flat-square&color=informational
    :target: https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE
    :alt: GPL-3.0 Licence

|Triumvirate| is made freely available under the `GPL-3.0 licence
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_ (or any later version).
Please see |Licence|_ (located at the repository directory root) for full
terms and conditions.

|Copyright| 2023 Mike S Wang & Naonori S Sugiyama


.. |Triumvirate| replace:: ``Triumvirate``

.. |hitomi| replace:: ``hitomi``
.. _hitomi: https://github.com/naonori/hitomi

.. |main| replace:: ``main``
.. _main: https://github.com/MikeSWang/Triumvirate/tree/main

.. |Makefile| replace:: ``Makefile``
.. _Makefile: https://github.com/MikeSWang/Triumvirate/blob/main/Makefile

.. |CitationCFF| replace:: ``CITATION.cff``
.. _CitationCFF: https://github.com/MikeSWang/Triumvirate/blob/main/CITATION.cff

.. |CitationMD| replace:: ``CITATION.md``
.. _CitationMD: https://github.com/MikeSWang/Triumvirate/blob/main/CITATION.md

.. |Licence| replace:: ``LICENCE``
.. _Licence: https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE

.. |Copyright| unicode:: U+000A9
