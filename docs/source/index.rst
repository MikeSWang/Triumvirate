.. title:: Home

.. figure:: _static/Triumvirate.png
    :align: center
    :width: 67 %

.. rubric:: **Three-Point Clustering Measurements in LSS**

|Triumvirate| is a Python/C++ package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS)
cosmological analyses.

.. tip::

    - The scientific context is explained in :doc:`background`.
    - Installation instructions can be found in :doc:`installation`.
    - To get started quickly, read :doc:`tutorials`.
    - For full API, please refer to :doc:`apiref`.
    - Frequently asked questions are answered in :doc:`faq`.


Attribution
===========

.. image:: https://joss.theoj.org/papers/10.21105/joss.05571/status.svg
    :target: https://doi.org/10.21105/joss.05571
    :alt: JOSS

.. image:: https://img.shields.io/badge/arXiv-2304.03643-b31b1b
    :target: https://arxiv.org/abs/2304.03643
    :alt: arXiv

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fsty3249-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/sty3249
    :alt: MNRAS

.. image:: https://img.shields.io/badge/10.1093%2Fmnras%2Fstx2333-grey?logo=doi
    :target: https://doi.org/10.1093/mnras/stx2333
    :alt: MNRAS

|br| To acknowledge the use of |Triumvirate| in your published research, please
cite the companion journal article :cite:p:`Wang:2023` and related publications
:cite:p:`Sugiyama:2018,Sugiyama:2019` with links above containing the
relevant information in the BibTeX format.


Acknowledgement
===============

.. figure:: _static/ERC-Logo-Flag.png
    :alt: ERC
    :align: left
    :width: 40%

|br| This project has received funding from the European Research Council (ERC)
under the European Union's Horizon 2020 research and innovation programme
(Grant agreement ID: `853291 <https://doi.org/10.3030/853291>`_).

Key underlying numerical algorithms were originally developed by
Naonori S Sugiyama, and are available in the GitHub repository |hitomi|_.

We thank the JOSS reviewers, William Coulton
(`@wcoulton <https://github.com/wcoulton>`_) and Alfonso Veropalumbo
(`@alfonso-veropalumbo <https://github.com/alfonso-veropalumbo>`_), for
their valuable feedback and suggestions, which have improved the
functionality and documentation of the code.


Contributing
============

User feedback and contributions are very welcome. Please refer to the
`contribution guidelines
<https://github.com/MikeSWang/Triumvirate/blob/main/CONTRIBUTING.md>`_.


Discussions & Wiki
==================

.. image:: https://img.shields.io/github/discussions/MikeSWang/Triumvirate
    :target: https://github.com/MikeSWang/Triumvirate/discussions
    :alt: Discussions

|br| A `community forum <https://github.com/MikeSWang/Triumvirate/discussions>`_
for users and developers is hosted on GitHub, where you can receive
announcements, post questions, share ideas and get updates.

A `wiki site <https://github.com/MikeSWang/Triumvirate/wiki>`_ collects wisdoms
for specific use cases and user environments.


Releases
========

Release notes are included in :doc:`releases`.


Licence
=======

.. image:: https://img.shields.io/github/license/MikeSWang/Triumvirate?label=licence&style=flat-square&color=informational
    :target: https://github.com/MikeSWang/Triumvirate/blob/main/LICENCE
    :alt: GPL-3.0 Licence

|br| |Triumvirate| is made freely available under the `GPL-3.0 licence
<https://www.gnu.org/licenses/gpl-3.0.en.html>`_.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |hitomi| replace:: ``hitomi``
.. _hitomi: https://github.com/naonori/hitomi


.. |br| raw:: html

    <br/>


.. toctree::
    :caption: Documentation
    :hidden:

    background
    installation
    tutorials
    apiref
    releases
    faq
