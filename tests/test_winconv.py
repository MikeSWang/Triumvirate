"""Test :mod:`~triumvirate.winconv`.

"""
import numpy as np
import pytest
from numpy.testing import assert_allclose

from triumvirate.winconv import (
    WinConvTerm,
    WinConvFormulae,
    ThreePointWindow,
    calc_threept_ic,
    TwoPointWinConvBase,
    TwoPCFWinConv,
    PowspecWinConv,
    ThreePointWinConvBase,
    ThreePCFWinConv,
    BispecWinConv,
)
