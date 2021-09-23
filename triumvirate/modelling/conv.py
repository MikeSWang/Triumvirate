"""
Window Convolution (:mod:`~triumvirate.modelling.conv`)
==========================================================================

Perform window convolution of (wide-angle corrected) bispectrum models.

"""
import numpy as np


def _interpolate_2d_grid(r_samples, f_samples, r_grid):
  """Interpolate samples of a bivariate function over a
  2-dimensional regular grid.

  Parameters
  ----------
  r_samples : array of float
    Sample points.
  f_samples : array of float
    Function samples.
  r_grid : array of float
    Grid points over which the function is interpolated.

  Returns
  -------
  f_grid : array of float
    Function values interpolated over the grid.

  """


def wide_angle_convolve():
