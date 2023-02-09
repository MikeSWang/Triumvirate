**********
Background
**********

Given a catalogue of discrete particles (such as galaxies) with their spatial
coordinates, |Triumvirate| computes estimators of the multipoles of the
three-point correlation function, also known as the bispectrum in Fourier
space, in the tri-polar spherical harmonic (TripoSH) decomposition
proposed by :cite:t:`Sugiyama:2019`. The objective of |Triumvirate| is to
provide efficient end-to-end measurement of clustering statistics which can be
fed into downstream galaxy survey analyses to constrain and test cosmological
models.

For the TripoSH decomposition mentioned above, |Triumvirate| offers the
following functionalities:

* :mod:`triumvirate.threept` measures three-point clustering statistics,
  namely multipoles of the bispectrum in Fourier space and of the three-point
  correlation function (3PCF) in configuration space, as in eqs. (41) & (48)
  of :cite:t:`Sugiyama:2019`;

* :mod:`triumvirate.twopt` measures two-point clustering statistics, namely
  multipoles of the power spectrum and two-point correlation function (2PCF),
  as in eqs. (40) & (49) of :cite:t:`Sugiyama:2018`;

* both modules mentioned above include the local plane-parallel estimator for
  a pair of survey-like data and random catalogues, and the global
  plane-parallel estimator for a cubic-box simulation, as in e.g. eq. (52)
  of :cite:t:`Sugiyama:2019`;

* both modules mentioned above can measure the configuration-space window
  function multipoles, which are used to convolve theoretical models derived
  in Fourier space through the Hankel transform, as in e.g. eqs. (58) & (63)
  of :cite:p:`Sugiyama:2019` and eq. (56) of :cite:t:`Sugiyama:2018`.

.. bibliography::
    :style: plain


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>


.. |hitomi| replace:: ``hitomi``
.. _hitomi: https://github.com/naonori/hitomi
