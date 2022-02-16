"""
Convolve bispectrum models with window functions.

"""
import glob
import inspect
import os
import sys
from pathlib import Path

import numpy as np

root_dir = Path(
    os.path.dirname(os.path.abspath(
        inspect.getframeinfo(inspect.currentframe()).filename
    ))
).parent

try:
    from triumvirate.modelling.conv import winconv_3pcf
except ImportError:
    sys.path.insert(0, str(root_dir))  # NOTE: Path-to-str conversion necessary
    from triumvirate.modelling.conv import winconv_3pcf


def load_model_bispectra(filedir, filepath_pattern):
    """Load model bispectra.

    Parameters
    ----------
    filedir : str or :class:`pathlib.Path`
        [description]
    filepath_pattern : str
        [description]

    Returns
    -------
    dict

    """
    bispectra = {}
    for filepath in glob.glob(str(Path(filedir)/filepath_pattern)):
        k1, k2, bk = np.loadtxt(filepath, unpack=True)

        ndim = int(np.sqrt(len(bk)))

        bispectra.update({os.path.basename(filepath): (
            k1.reshape((ndim, ndim)),
            k2.reshape((ndim, ndim)),
            bk.reshape((ndim, ndim))
        )})

    return bispectra


if __name__ == '__main__':

    bispectra_ = load_model_bispectra(
        root_dir/"storage"/"input"/"models"/"DiDio",
        "BISP_num_z0_Newt_WA0_*_short.dat"
    )

    bispectra = {}
    for key, val in bispectra_.items():
        ell1, ell2, ELL = list(map(int, np.asarray(key.split('_'))[[-4, -3, -2]]))
        bispectra[(ell1, ell2, ELL)] = val
