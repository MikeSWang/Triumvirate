.. title:: Triumvirate Documentation

========================================================================
**Three-Point Clustering Correlators in Large-Scale Structure Analysis**
========================================================================

.. image:: https://img.shields.io/badge/arXiv-yymm.nnnnn-important
    :target: https://arxiv.org/abs/yymm.nnnnn
    :alt: arXiv eprint
.. image:: https://img.shields.io/github/v/release/MikeSWang/Triumvirate?label=release
    :target: https://github.com/MikeSWang/Triumvirate/releases/latest
    :alt: GitHub release (latest by date)
.. image:: https://readthedocs.org/projects/triumvirate/badge/?version=latest
    :target: https://triumvirate.readthedocs.io/en/latest
    :alt: Documentation status
.. image:: https://travis-ci.com/MikeSWang/Triumvirate.svg?branch=main
    :target: https://travis-ci.com/MikeSWang/Triumvirate
    :alt: Build status
.. image:: https://img.shields.io/badge/licence-GPLv3-informational
    :target: https://github.com/MikeSWang/Triumvirate/tree/main/LICENCE
    :alt: Licence

|Triumvirate| is a Python/C++ software package for the measurement and
model comparison of three-point clustering correlators in large-scale
structure analysis.

.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


Pre-release guide
=================

.. note::

    If you can see this text, please note that the repository is currently
    in pre-release mode, with specific instructions detailed in this
    section overriding any other instructions found elsewhere.


After ``git clone``-ing this repository, at the root of its directory
run the following commands to install the code::

    make clean
    make install

If you wish to use the Python (or C++) package only, replace ``make all``
with ``make pyinstall`` (or ``make cppinstall``).  Depending on your
machine environment, certain parts of the |Makefile|_ may need replacing
(e.g. the C++ compiler and compilation flags).

See Python scripts/notebooks in |examples|_ for simple examples of how to
perform basic tasks.  For the latest API documentation, you may need to
build it yourself using ``sphinx-autodoc``.

.. todo::

    Add instructions for generating API docs here.


Please note that both testing and documentation are incomplete and there
may be unknown bugs.  Any feedback is welcome, and more importantly,
if you would like to contribute to testing and documentation please get
in touch.

.. |Makefile| replace:: ``Makefile``
.. _Makefile: https://github.com/MikeSWang/Triumvirate/tree/main/Makefile

.. |examples| replace:: ``examples``
.. _examples: https://github.com/MikeSWang/Triumvirate/tree/main/examples


Installation
============


Documentation
=============

- :doc:`tutorials` (under construction): notebooks showcasing the use of
  |Triumvirate| will be gradually added, so look out for any updates!
  For now, |examples|_ contains some scripts/notebooks that demonstrates
  how basic tasks can be performed.
- :doc:`apidoc`: comprehensive API documentation.


Attribution
===========

The underlying C++ code was originally developed by Naonori Sugiyama.
Please refer to/cite Sugiyama et al. (2019)
`[doi:10.1093/mnras/sty3249] <https://doi.org/10.1093/mnras/sty3249>`_
accordingly.


Licence
=======

|Triumvirate| is made freely available under the `GPL v3.0
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_ licence.

Copyright 2021–22.

.. toctree::
    :hidden:

    tutorials
    apidoc
