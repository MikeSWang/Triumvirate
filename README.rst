..
    This read-me .rst file is for PyPI project description only, and
    should be periodically compared against the official read-me .md file
    which is rendered on GitHub and included in documentation.

.. figure:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/image/Triumvirate.png
    :class: dark-light
    :align: center
    :width: 67 %

========================================
Three-Point Clustering Statistics in LSS
========================================

|Triumvirate| is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics and performing survey window convolution in
large-scale structure (LSS) cosmological analyses.


Documentation
=============

.. image:: https://img.shields.io/badge/Read%20the%20Docs-stable-informational?logo=ReadtheDocs
    :target: https://triumvirate.readthedocs.io/en/stable/
    :alt: Documentation

Comprehensive documentation including the `scientific background
<https://triumvirate.readthedocs.io/en/stable/background.html>`_,
`installation instructions
<https://triumvirate.readthedocs.io/en/stable/installation.html>`_,
`tutorials
<https://triumvirate.readthedocs.io/en/stable/tutorials.html>`_ and
`API reference
<https://triumvirate.readthedocs.io/en/stable/apiref.html>`_
can be found at `triumvirate.readthedocs.io
<https://triumvirate.readthedocs.io/en/stable/>`_.


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


GPU variants
------------

.. .. image:: https://img.shields.io/pypi/v/Triumvirate-CUDA?logo=PyPI&color=informational
..     :target: https://pypi.org/project/Triumvirate-CUDA
..     :alt: PyPI

.. .. image:: https://img.shields.io/conda/v/msw/triumvirate-cuda?logo=Anaconda&color=informational
..     :target: https://anaconda.org/msw/triumvirate-cuda
..     :alt: Conda

CUDA variants of the Python package are/will be made available as
``Triumvirate-CUDA`` on `PyPI <https://pypi.org/project/Triumvirate-CUDA>`_
and ``triumvirate-cuda`` through `Conda
<https://anaconda.org/msw/triumvirate-cuda>`_.


Attribution
===========

.. image:: https://joss.theoj.org/papers/10.21105/joss.05571/status.svg
    :target: https://doi.org/10.21105/joss.05571
    :alt: JOSS

.. image:: https://img.shields.io/badge/10.1088%2F1475-7516%2F2025%2F06%2F031-grey?logo=doi
    :target: https://doi.org/10.1088/1475-7516/2025/06/031
    :alt: JCAP

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fsty3249-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/sty3249
    :alt: MNRAS

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fstx2333-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/stx2333
    :alt: MNRAS

.. image:: https://img.shields.io/badge/arXiv-2304.03643-b31b1b
    :target: https://arxiv.org/abs/2304.03643
    :alt: arXiv

.. image:: https://img.shields.io/badge/arXiv-2411.14947-b31b1b
    :target: https://arxiv.org/abs/2411.14947
    :alt: arXiv

.. image:: https://img.shields.io/badge/zenodo-10.5281%2Fzenodo.10072128-1682D4
    :target: https://doi.org/10.5281/zenodo.10072128
    :alt: Zenodo

To acknowledge the use of |Triumvirate| in your published research, please
cite the relevant publications linked above; for convenience, you can refer to
the files |CitationCFF|_ and |CitationMD|_ for the relevant information in
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

.. image:: https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fgithub.com%2FMikeSWang%2FTriumvirate%2Fraw%2Fmain%2Fdeploy%2Fpkg%2Fpyproject%2F.pyproject.toml&logo=python
   :alt: Python-Version

.. image:: https://img.shields.io/badge/std-c%2B%2B17-informational?logo=cplusplus
   :alt: C++-Standard

.. .. image:: https://img.shields.io/github/release-date-pre/MikeSWang/Triumvirate
..     :target: https://github.com/MikeSWang/Triumvirate/releases/latest
..     :alt: Release-Date

.. .. image:: https://img.shields.io/github/commits-since/MikeSWang/Triumvirate/latest/main
..     :alt: Commits-Since

.. .. image:: https://img.shields.io/github/issues/MikeSWang/Triumvirate/build
..     :target: https://github.com/MikeSWang/Triumvirate/issues?q=is%3Aopen+is%3Aissue+label%3Abuild
..     :alt: Build-Issues

.. .. image:: https://img.shields.io/github/issues/MikeSWang/Triumvirate/bug
..     :target: https://github.com/MikeSWang/Triumvirate/issues?q=is%3Aopen+is%3Aissue+label%3Abug
..     :alt: Bug-Issues

.. .. image:: https://img.shields.io/github/issues/MikeSWang/Triumvirate/feature
..     :target: https://github.com/MikeSWang/Triumvirate/issues?q=is%3Aopen+is%3Aissue+label%3Afeature
..     :alt: Feature-Issues

.. .. image:: https://img.shields.io/github/issues-pr/MikeSWang/Triumvirate
..     :target: https://github.com/MikeSWang/Triumvirate/pulls
..     :alt: Pull-Requests

.. .. image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit
..     :target: https://github.com/pre-commit/pre-commit
..     :alt: pre-commit

.. .. image:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/image/GitHub-Codespace-badge.png
..     :target: https://codespaces.new/MikeSWang/Triumvirate?hide_repo_select=true&ref=main
..     :alt: Codespaces
..     :width: 249px

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
    :alt: GPLv3+ Licence

|Triumvirate| is made freely available under the `GPLv3+ licence
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_.
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
