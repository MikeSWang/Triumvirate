---
title: 'Triumvirate: A Python/C++ package for three-point clustering measurements'
tags:
  - Python
  - C++
  - Cython
  - cosmology
  - astronomy
authors:
  - name: Mike Shengbo Wang
    orcid: 0000-0002-2652-4043
    affiliation: 1       # (add quotes to multiple affiliations)
    corresponding: true
  - name: Florian Beutler
    orcid: 0000-0003-0467-5438
    affiliation: 1
  - name: Naonori S. Sugiyama
    affiliation: 2
affiliations:
  - name: Institute for Astronomy, University of Edinburgh,
          Royal Observatory Edinburgh, Blackford Hill,
          Edinburgh EH9 3HJ, United Kingdom
    index: 1
  - name: National Astronomical Observatory of Japan,
          Mitaka, Tokyo 181-8588, Japan
    index: 2
date: 1 February 2023
bibliography: paper.bib


---

# Summary

`Triumvirate` is a Python/C++ package for measuring the three-point
clustering statistics in large-scale structure (LSS) cosmological analyses.
Given a catalogue of discrete particles (such as galaxies) with their
spatial coordinates, it computes estimators of the multipoles of the
three-point correlation function, also known as the bispectrum in Fourier
space, in the tri-polar spherical harmonic (TripoSH) decomposition
proposed by @Sugiyama:2019. The objective of `Triumvirate` is to provide
efficient end-to-end measurement of clustering statistics which can be fed
into downstream galaxy survey analyses to constrain and test cosmological
models. To this end, it builds upon the original algorithms in the `hitomi`
code[^1] developed by Sugiyama et al. [-@Sugiyama:2018;-@Sugiyama:2019],
and supplies a user-friendly interface with flexible input/output (I/O) of
catalogue data and measurement results, with the built program
configurable through external parameter files and tracked through enhanced
logging and warning/exception handling. For completeness and
complementarity, methods for measuring two-point clustering statistics are
also included in the package.


[^1]: [github.com/naonori/hitomi](https://github.com/naonori/hitomi)


# Statement of need

The analysis of higher-order clustering statistics is a key pursuit
of the current and forthcoming generations of large galaxy surveys such as
the Dark Energy Spectroscopic Instrument (DESI)[^2] [@DESI:2016] and
_Euclid_[^3] [@Euclid:2011]. Although the matter density fluctuations in
the Universe have been observed to be almost Gaussian on large scales,
primordial non-Gaussianity in the initial conditions of structure
formation and late-time non-linear gravitational dynamics can both leave
potentially detectable signals in the higher-order moments of the galaxy
distribution on very large and small scales [@Bernardeau:2002]. Therefore
any measurement of the three-point clustering statistics, the leading
non-Gaussian moment, offers a promising probe of both early-Universe and
gravitational physics. In addition, it complements two-point clustering
statistics in constraining cosmological models by breaking down certain
parameter degeneracies [@Sefusatti:2006].

In contrast to the two-point statistic analysis which has become standard
in recent galaxy surveys [e.g. @BOSS:2017; @eBOSS:2021], three-point
clustering statistics have more degrees of freedom and thus can be
compressed in a greater number of ways with different choices of the
coordinate system. For the TripoSH decomposition mentioned above, there
is a need for a computational program that is easy to use, versatile and
suited for the large data sets expected from modern galaxy surveys, and
`Triumvirate` is designed to meet that demand. More specifically, it
can compute:

  * three-point clustering statistics, namely multipoles of the bispectrum
    in Fourier space and of the three-point correlation function (3PCF)
    in configuration space [@Sugiyama:2019];

  * two-point clustering statistics, namely multipoles of the power
    spectrum and two-point correlation function (2PCF) [@Sugiyama:2018];

  * for both three- and two-point statistics, the local plane-parallel
    estimator for a pair of survey-like data and random catalogues, and
    the global plane-parallel estimator for a cubic-box simulation
    [@Feldman:1994;@Yamamoto:2006;@Sugiyama:2019];

  * for both three- and two-point statistics, the configuration-space
    window function multipoles, which are used to convolve theoretical
    models derived in Fourier space through the Hankel transform
    [@Wilson:2016;@Sugiyama:2019].

For the global plane-parallel estimators, the simulation box is placed at
the spatial infinity (or equivalently the observer is), so that the
line of sight to each particle can be treated as the same and taken to be
along the $z$-axis. For the local plane-parallel estimators, the observer
is placed at the origin in the survey coordinates, and the line of sight
is chosen to point towards one of the particles in a triplet or pair for
three- or two-point clustering measurements respectively.

The geometry of the survey leaves an imprint on the clustering statistics,
where in Fourier space the effect is a convolution with the survey window
function. This convolution mixes different multipoles of the underlying
clustering statistics and the survey window, and the precise convolution
formula (i.e. the number of multipoles to include in modelling) needed to
achieve a given level of convergence depends on the precise survey geometry
including any sample weights applied. Therefore the functionality to
measure the window function is an integral part of this program.

These functionalities are essential to cosmological inference pipelines,
and can help validate any analytical covariance matrix predictions against
sample estimates. Since precise covariance matrix estimates usually
require clustering measurements repeated over a large number of simulated
mock catalogues, computational efficiency is an important objective.

Finally, `Triumvirate` also enables comparison studies between
alternative compressed statistics of three-point clustering, which may have
different constraining power on different cosmological parameters. There
are existing software packages for some of these alternative approaches:

  * `pylians`[^4] [@VillaescusaNavarro:2018] computes the bispectrum with
    the Scoccimarro estimator [@Scoccimarro:2015] for triangle
    configurations parametrised by two wavenumbers and the angle between
    the corresponding wavevectors;

  * `nbodykit`[^5] [@Hand:2018] computes the isotropised 3PCF with a
    pair-counting algorithm [@Slepian:2015], although in principle
    this can be generalised to anisotropic 3PCF [@Slepian:2018], or be
    implemented using FFTs [@Slepian:2016];

  * @Philcox:2021[^6] advocates a windowless cubic estimator for the
    bispectrum in the Soccimarro decomposition, which can be evaluated
    using FFTs. However, this approach requires the inversion of a Fisher
    matrix obtained from a suite of Monte Carlo realisations.

As these programs use a different decomposition of three-point clustering
statistics and focus on either configuration- or Fourier-space statistics
only, `Triumvirate` fulfills complementary needs in current galaxy
clustering analyses.


[^2]: [desi.lbl.gov](https://www.desi.lbl.gov)
[^3]: [sci.esa.int/euclid](https://sci.esa.int/web/euclid/),
      [euclid-ec.org](https://www.euclid-ec.org)
[^4]: [pylians3.readthedocs.io](https://pylians3.readthedocs.io)
[^5]: [nbodykit.readthedocs.io](https://nbodykit.readthedocs.io)
[^6]: [github.com/oliverphilcox/Spectra-Without-Windows](https://github.com/oliverphilcox/Spectra-Without-Windows)


# Implementation

Direct calculation of Fourier modes of density field fluctuations
as a sum over a large number of particles is computationally infeasible,
but the TripoSH estimators can be cast in a form amenable to
fast Fourier transforms (FFTs), which can be utilised to speed up
evaluations. In our numerical scheme, the particles are assigned to
regular mesh grids with appropriate weighting, which is a combination of
spherical harmonics and weights from the input catalogues.
Fourier-space fields are obtained by FFTs over the cubic mesh grids, and
clustering statistics are formed by multiplying the discretely sampled and
transformed fields grid-by-grid, before binning in spherical shells.
For the three-point statistics, the shot noise components are effectively
two-point statistics and can be calculated the same way.

`Triumvirate` supports mesh assignment schemes up to order 4, namely the
piecewise cubic spline (PCS) scheme, where order $p$ refers to the number
of grid points a particle is assigned to. Since mesh assignment convolves
the underlying field with a sampling window, the transformed Fourier-space
field should be compensated by dividing out the sampling window
[@Hockney:1988]. For power spectrum measurements only where no inverse FFT
is involved the interlacing technique can be used to reduce the amount of
aliasing [@Sefusatti:2016], an artifact of discrete Fourier transform where
the sampled Fourier mode at each wavenumber receives contributions from
other modes, with the effect increasingly prominent as the Nyquist
wavenumber is approached; without interlacing, the corrections from
@Jing:2005 (see eq. 20 therein) are adopted instead. For all clustering
statistics, increasing the mesh assignment order and/or the number
of grid cells can help reduce aliasing (at the expense of speed and
memory).


# Features

Besides code refactoring, `Triumvirate` has many value-added features
in comparison with the predecessor `himoti`:

* The frontend is written in Python for interactivity and convenience,
  with Cython binding the C++ backend. Although the user will typically
  use the Python interface, the C++ code can also be compiled and executed
  independently.

* Measurement pipelines can be configured through external parameter files
  (in the YAML format for the Python program), cleanly separating user
  inputs from the program itself. Alternatively, measurement parameters can
  be set for individual Python methods without the use of a parameter file.

* The reading of catalogue data is implemented via `astropy.io`
  [@Astropy:2022] and `nbodykit` [@Hand:2018], with flexible support for
  different file formats such as text and `fits` files.

* Numerical algorithms are parallelised with OpenMP, with for-loops over
  catalogue particles and mesh grid cells distributed amongst multiple
  CPU threads.

* Mesh assignment schemes from order $p = 1$ to $4$ are supported:
  nearest grid point (NGP), cloud-in-cell (CIC), triangular-shape cloud
  (TSC) and piecewise cubic spline (PCS).

* Interlacing is supported for power spectrum measurements.

* Two normalisation choices are implemented for all clustering statistics,
  one as a sum of catalogue particles [@Feldman:1994, eq. 2.4.1] and
  another as a sum over the mesh grid [@Sugiyama:2019, eq. 37].

* A customised logger is provided for runtime tracking, with enhanced
  handling of warnings and exceptions for parameter and data I/O.


# Performance

When a large number of grid cells, $N_\mathrm{mesh}$, are used to
sample the density fields from a catalogue of $N_\mathrm{part}$ particles
on a mesh with $N_\mathrm{mesh} \gg N_\mathrm{part}$, the dominant
operations are FFTs with complexity
$\mathcal{O}\left({N_\mathrm{mesh} \ln N_\mathrm{mesh}}\right)$.
Therefore the complexity for three-point clustering measurements is
$\mathcal{O}\left({N_\mathrm{bin}^2 N_\mathrm{mesh} \ln N_\mathrm{mesh}}\right)$,
where $N_\mathrm{bin}$ is the number of coordinate bins.

It is worth noting that in `Triumvirate`, the spherical harmonic weights
are applied to individual particles rather than the mesh grids. This should
result in more accurate results at the expense of memory usage, as multiple
meshes need to be stored for spherical harmonics of different degrees
and orders. We estimate the minimum memory usage for bispectrum
measurements to be $11 M$ and $9 M$ respectively for local and global
plane-parallel estimators, where $M = 16 N_\mathrm{mesh}$ bytes (roughly
$1.5\times10^{-8} N_\mathrm{mesh}$ gibibytes[^7]); for local and global
plane-parallel 3PCF estimators, the figures are $10 M$ and $9 M$
respectively.

In the table below, we show the wall time and peak memory usage for
bispectrum and 3PCF measurements of a few select multipoles and
grid numbers with $N_\mathrm{bin} = 20$, using a single core on one
AMD EPYC 7H12 processor with base frequency 2.60 GHz. With multithreading
enabled, the run time is reduced (see the last column in the table). Here
'lpp' and 'gpp' denote local and global plane-parallel approximations
respectively. For the global plane-parallel estimates, the catalogue used
is a cubic box containing $N_\mathrm{part} = 8 \times 10^6$ particles;
for the local plane-parallel estimates, the data and random catalogues
contain $N_\mathrm{part} = 6.6 \times 10^5$ and $1.3 \times 10^7$ particles
respectively. Since both the bispectrum and 3PCF are computed with FFTs,
the computation time and memory usage for them are roughly the same, with
minor differences due to the slightly different number of mesh grids needed
for evaluation.

--------------------------------------------------------------------------------------------------------------------------------------------------
Multipole/$N_\mathrm{mesh}$                        $128^3$                     $256^3$                     $512^3$            $512^3$ (32 threads)
------------------------------ --------------------------- --------------------------- --------------------------- -------------------------------
$B_{000}^{\mathrm{(lpp)}}$                   96 s, 1.8 GiB              215 s, 4.4 GiB              1247 s, 25 GiB                    85 s, 25 GiB

$B_{000}^{\mathrm{(gpp)}}$                   42 s, 0.7 GiB              172 s, 2.9 GiB              1185 s, 21 GiB                    59 s, 21 GiB

$B_{202}^{\mathrm{(lpp)}}$                  320 s, 1.8 GiB             1030 s, 4.4 GiB              6449 s, 25 GiB                   267 s, 25 GiB

$B_{202}^{\mathrm{(gpp)}}$                   42 s, 0.7 GiB              176 s, 2.9 GiB              1187 s, 21 GiB                    60 s, 21 GiB

$\zeta_{000}^{\mathrm{(lpp)}}$               90 s, 1.8 GiB              211 s, 4.4 GiB              1403 s, 21 GiB                    83 s, 21 GiB

$\zeta_{000}^{\mathrm{(gpp)}}$               43 s, 0.7 GiB              178 s, 2.9 GiB              1226 s, 19 GiB                    55 s, 19 GiB

$\zeta_{202}^{\mathrm{(lpp)}}$              267 s, 1.8 GiB              964 s, 4.4 GiB              6377 s, 21 GiB                   266 s, 21 GiB

$\zeta_{202}^{\mathrm{(gpp)}}$               43 s, 0.7 GiB              177 s, 2.9 GiB              1241 s, 19 GiB                    57 s, 19 GiB
--------------------------------------------------------------------------------------------------------------------------------------------------


[^7]: Note that 1 gibibyte (GiB) is $2^{30}$ bytes, as opposed to
      1 gigabyte (GB) which is $10^9$ bytes. GiB is the preferred unit
      by job schedulers such as Slurm for computer clusters.


# Future work

`Triumvirate` will be routinely maintained and updated depending on
user feedback. One extension of interest is the inclusion of other
three-point clustering estimators with different coordinate systems and
compression choices, and the functionality to transform between them.
The ability to measure clustering statistics from a density field already
sampled on a mesh grid may also be useful. In addition, porting the code
to graphic processing units (GPUs) can bring further parallelisation that
can enhance the performance of the code.


# Acknowledgements

This project has received funding from the European Research Council (ERC)
under the European Union’s Horizon 2020 research and innovation programme
(grant agreement 853291). FB is a Royal Society University Research
Fellow.

We thank the reviewers for their valuable feedback and suggestions, which
have improved the functionality and documentation of this code.

# References
