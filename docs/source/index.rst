.. title:: Triumvirate Documentation

.. figure:: _static/Triumvirate.png
    :align: center
    :width: 67 %

##########################################################################
**Three-Point Clustering Measurements in LSS**
##########################################################################

|Triumvirate| is a Python/C++ package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS)
cosmological analyses.

.. tip::

    * The scientific context is explained in :doc:`background`.
    * Installation instructions can be found in :doc:`installation`.
    * To get started quickly, read :doc:`tutorials`.
    * For full API, please refer to :doc:`apiref`.


Attribution
===========

.. image:: https://img.shields.io/badge/JOSS-doi-brightgreen
    :target: https://joss.theoj.org/papers/?/status.svg
    :alt: JOSS
.. image:: https://img.shields.io/badge/arXiv-yymm.%3F-b31b1b?style=flat-square
    :target: https://arxiv.org/abs/?.?
    :alt: arXiv
.. image:: https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fsty3249-blue?style=flat-square
    :target: https://doi.org/10.1093/mnras/sty3249
    :alt: MNRAS
.. image:: https://img.shields.io/badge/doi-10.1093%2Fmnras%2Fstx2333-blue?style=flat-square
    :target: https://doi.org/10.1093/mnras/stx2333
    :alt: MNRAS

|br| To acknowledge the use of |Triumvirate| in your published research, please
cite the publications :cite:p:`Sugiyama:2018,Sugiyama:2019` linked above which
contain the relevant information in the BibTeX format.


Acknowledgement
===============

This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(grant agreement 853291).

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository |hitomi|_.


Contributing
============

User feedback and contributions are highly valued and very welcome. Please
refer to the `contribution guidelines
<https://github.com/MikeSWang/Triumvirate/blob/main/CONTRIBUTING.md>`_.


Releases
========

Changes in recent releases are listed in the `change log
<https://github.com/MikeSWang/Triumvirate/blob/main/CHANGELOG.md>`_.


Licence
=======

.. image:: https://img.shields.io/github/license/MikeSWang/Triumvirate?color=informational&style=flat-square
    :target: https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE
    :alt: GPL-3.0 Licence

|Triumvirate| is made freely available under the `GPL-3.0 licence
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |hitomi| replace:: ``hitomi``
.. _hitomi: https://github.com/naonori/hitomi


.. |br| raw:: html

    <br/>


.. toctree::
    :hidden:

    background
    installation
    tutorials
    apiref
