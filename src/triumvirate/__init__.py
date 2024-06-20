"""
##########################################################################
``Triumvirate`` | Three-Point Clustering Measurements in LSS
##########################################################################

``Triumvirate`` is a Python/C++ package for measuring three-point (and
two-point) clustering statistics in large-scale structure (LSS)
cosmological analyses.

.. topic:: GNU General Public License, version 3

    Copyright (C) 2023, Mike S Wang & Naonori S Sugiyama

    This program is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License along
    with this program.  If not, see `<https://www.gnu.org/licenses/>`_.

"""
from importlib.metadata import PackageNotFoundError, version

# Installation validation
from ._valid_install import validate_installation  # noqa: F401


__copyright__ = 'Copyright 2023, Mike S Wang & Naonori S Sugiyama'
__date__ = '2023-10-04'
__license__ = 'GPL-3.0'

try:
    __version__ = version('triumvirate')
except PackageNotFoundError:
    __version__ = '0.4.0'  # UPDATE: fallback version number
