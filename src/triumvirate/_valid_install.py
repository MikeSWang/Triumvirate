"""
Installation Validation (:mod:`~triumvirate._valid_install`)
==========================================================================

Check whether the Cython extensions have been compiled correctly by
performing test computations.

It is not meant to replace the full `pytest` suite, and does not check the
correctness of the computed results.

"""
import warnings

import numpy as np

from .catalogue import ParticleCatalogue
from .dataobjs import Binning
from .parameters import ParameterSet, fetch_paramset_template
from .twopt import compute_powspec_in_gpp_box


def _create_minimal_catalogue():
    """Create a minimal catalogue of a single particle at a unit length
    away from all positive axes.

    Returns
    -------
    :class:`~triumvirate.catalogue.ParticleCatalogue`
        Minimal catalogue.

    """
    CTLG_LENGTH = 1

    x = y = z = np.ones(CTLG_LENGTH)

    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message=".*'nz' field is None.*")
        return ParticleCatalogue(x, y, z)


def _define_minimal_binning():
    """Define a minimal binning scheme.

    Returns
    -------
    :class:`~triumvirate.dataobjs.Binning`
        Minimal binning.

    """
    SPACE, SCHEME = 'fourier', 'lin'
    BIN_MIN, BIN_MAX, NUM_BINS = 0., .1, 2

    return Binning(
        SPACE, SCHEME, bin_min=BIN_MIN, bin_max=BIN_MAX, num_bins=NUM_BINS
    )


def _create_mock_paramset():
    """Create a mock parameter set.

    Returns
    -------
    :class:`~triumvirate.parameters.ParameterSet`
        Mock parameter set.

    """
    tmpl_paramdict = fetch_paramset_template('dict')

    for ax_name in ['x', 'y', 'z']:
        tmpl_paramdict['boxsize'][ax_name] = 2.
        tmpl_paramdict['ngrid'][ax_name] = 2.

    tmpl_paramdict.update({
        'degrees': {'ell1': 0, 'ell2': 0, 'ELL': 0},
        'range': [0., 1.],
        'num_bins': 2,
        'verbose': 60,
    })

    return ParameterSet(param_dict=tmpl_paramdict)


def validate_installation():
    """Validate installation by computing a minimal test case.

    Returns
    -------
    True
        If the installation is valid.

    """
    catalogue = _create_minimal_catalogue()
    binning = _define_minimal_binning()
    paramset = _create_mock_paramset()

    compute_powspec_in_gpp_box(catalogue, binning=binning, paramset=paramset)

    print("Installation of Triumvirate has been validated.")

    return True
