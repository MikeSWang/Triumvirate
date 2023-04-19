..
    This read-me .rst file is for PyPI project description only, and
    should be periodically compared against the official read-me .md file
    which is rendered on GitHub and included in documentation.

.. figure:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/Triumvirate.png
    :align: center
    :width: 67 %

==========================================
Three-Point Clustering Measurements in LSS
==========================================

.. image:: https://img.shields.io/github/actions/workflow/status/MikeSWang/Triumvirate/ci.yml?label=ci&logo=GitHubActions
    :target: https://github.com/MikeSWang/Triumvirate/actions/workflows/ci.yml
    :alt: CI

.. image:: https://img.shields.io/readthedocs/triumvirate/latest?logo=ReadtheDocs
    :target: https://triumvirate.readthedocs.io/en/latest/
    :alt: Docs

.. image:: https://img.shields.io/github/v/release/MikeSWang/Triumvirate?display_name=tag&sort=semver&logo=Git
    :target: https://github.com/MikeSWang/Triumvirate/releases/latest
    :alt: Release

|Triumvirate| is a Python/C++ software package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS) cosmological
analyses.


Documentation
=============

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

.. image:: https://img.shields.io/conda/vn/msw/triumvirate?logo=Anaconda&color=informational
    :target: https://anaconda.org/msw/triumvirate
    :alt: conda

|Triumvirate| as a Python package is distributed through
`PyPI <https://pypi.org/project/Triumvirate>`_ and
`conda <https://anaconda.org/msw/triumvirate>`_. Instructions for installation
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

Then, execute in terminal:

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

.. image:: https://joss.theoj.org/papers/a8325e3897dd726d9df42286bf72d19f/status.svg
    :target: https://joss.theoj.org/papers/a8325e3897dd726d9df42286bf72d19f
    :alt: JOSS

.. image:: https://img.shields.io/badge/arXiv-2304.03643-b31b1b
    :target: https://arxiv.org/abs/2304.03643
    :alt: arXiv

.. image:: https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fsty3249-informational
    :target: https://doi.org/10.1093/mnras/sty3249
    :alt: MNRAS

.. image:: https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fstx2333-informational
    :target: https://doi.org/10.1093/mnras/stx2333
    :alt: MNRAS

To acknowledge the use of |Triumvirate| in your published research, please
cite the publications linked above which contain the relevant information
in the BibTeX format.


Acknowledgement
===============

.. figure:: https://github.com/MikeSWang/Triumvirate/raw/main/docs/source/_static/ERC-Logo-Flag.png
    :alt: ERC
    :align: left
    :width: 40%

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(Grant agreement No. `853291 <https://doi.org/10.3030/853291>`_).

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository |hitomi|_.


Contributing
============

User feedback and contributions are very welcome. Please refer to the
`contribution guidelines
<https://github.com/MikeSWang/Triumvirate/blob/main/CONTRIBUTING.md>`_.


Discussions
===========

A `community forum <https://github.com/MikeSWang/Triumvirate/discussions>`_
for users and developers is hosted on GitHub, where you can receive
announcements, post questions, share ideas and get updates.


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
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_. Please see |Licence|_
(located at the repository directory root) for full terms and conditions.

&copy; 2023 Mike S Wang & Naonori S Sugiyama


.. |Triumvirate| replace:: ``Triumvirate``

.. |hitomi| replace:: ``hitomi``
.. _hitomi: https://github.com/naonori/hitomi

.. |main| replace:: ``main``
.. _main: https://github.com/MikeSWang/Triumvirate/tree/main

.. |Makefile| replace:: ``Makefile``
.. _Makefile: https://github.com/MikeSWang/Triumvirate/blob/main/Makefile

.. |Licence| replace:: ``Licence``
.. _Licence: https://github.com/MikeSWang/Triumvirate/blob/main/Licence
