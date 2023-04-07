***
FAQ
***

Have you got a question? See if it has been asked before and answered below.

- *Can I specify multiple clustering statistic multipoles in one parameter
  file or one function call?*

No, unfortunately.

Unlike some other programs which assign a catalogue of particles to a
single mesh grid and subsequently apply the spherical harmonic weighting
(as a function of the line of sight) on the grid, |Triumvirate| apply
the spherical harmonic weighting to each particle before mesh assignment,
which should be more accurate (as it uses the actual line of sight rather
than the position vector of the grid cells).

Since different multipoles require different spherical harmonic weighting,
storing a large number of differently sampled mesh grids for each spherical
harmonic factor is memory-intensive. Instead, the mesh grids are recomputed
for each multipole, and there is no time saving from computing multiple
multipoles in one run.

Therefore, for the sake of syntax simplicity, only one multipole can be
specified in the parameter file or in one function call.


.. |Triumvirate| raw:: html

    <span style="font-variant: small-caps">Triumvirate</span>
