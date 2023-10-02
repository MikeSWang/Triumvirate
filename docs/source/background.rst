**********
Background
**********

.. seealso::

    :cite:t:`Wang:2023` for the companion journal article.

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

- :mod:`triumvirate.threept` measures three-point clustering statistics,
  namely multipoles of the bispectrum in Fourier space and of the three-point
  correlation function (3PCF) in configuration space, as in eqs. (41) & (48)
  of :cite:t:`Sugiyama:2019`;

- :mod:`triumvirate.twopt` measures two-point clustering statistics, namely
  multipoles of the power spectrum and two-point correlation function (2PCF),
  as in eqs. (40) & (49) of :cite:t:`Sugiyama:2018`;

- both modules mentioned above include the local plane-parallel estimator for
  a pair of survey-like data and random catalogues, and the global
  plane-parallel estimator for a simulation-like catalogue in a cuboid box,
  as in e.g. eq. (52) of :cite:t:`Sugiyama:2019`;

- both modules mentioned above can measure the configuration-space window
  function multipoles, which are used to convolve theoretical models derived
  in Fourier space through the Hankel transform, as in e.g. eqs. (58) & (63)
  of :cite:t:`Sugiyama:2019` and eq. (56) of :cite:t:`Sugiyama:2018`.

For the global plane-parallel estimators, the simulation box is placed at the
spatial infinity (or equivalently the observer is), so that the line of sight
to each particle can be treated as the same and taken to be along
the :math:`z`-axis.

For the local plane-parallel estimators, the observer is placed at the origin
in the survey coordinates, and the line of sight is chosen to point towards
one of the particles in a doublet/triplet for two-/three-point clustering
measurements.

The geometry of the survey leaves an imprint on the clustering statistics,
where in Fourier space the effect is a convolution with the survey window
function. This convolution mixes different multipoles of the underlying
clustering statistics and the survey window, and the precise convolution
formula (i.e. the number of multipoles to include in modelling) needed to
achieve a given level of convergence depends on the precise survey geometry
including any sample weights applied.

.. bibliography::
    :style: plain


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
