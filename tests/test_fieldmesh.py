"""Test :mod:`~triumvirate.fieldmesh`.

"""
import numpy as np
import pytest

from triumvirate.fieldmesh import record_binned_vectors


@pytest.mark.parametrize(
    "case, statistic_type, count, last_entry",
    [
        (
            'default', 'modes', 19548,
            (3, 0.08, 0.105, 0.10053096, 0.02513274, 0.01256637),
        ),
        (
            'override', 'modes', 2468,
            (3, 0.08, 0.105, 0.10053096, 0.02513274, 0.01256637),
        ),
        (
            'empty', 'pairs', 228396,
            (3, 125., 150., 148.4375, 19.53125, 7.8125),
        ),
    ]
)
def test_record_binned_vectors(case, statistic_type, count, last_entry,
                               request):

    if case == 'default':
        test_paramset = request.getfixturevalue('test_paramset')
        test_paramset.update(statistic_type=statistic_type)
        boxsize = None
        ngrid = None
    elif case == 'override':
        test_paramset = request.getfixturevalue('test_paramset')
        test_paramset.update(statistic_type=statistic_type)
        boxsize = 500.
        ngrid = 128
    elif case == 'empty':
        test_paramset = None
        boxsize = 500.
        ngrid = 128
    else:
        raise ValueError("Invalid test case for `record_binned_vectors`.")

    if statistic_type == 'modes':
        test_binning = request.getfixturevalue('test_binning_fourier')
    elif statistic_type == 'pairs':
        test_binning = request.getfixturevalue('test_binning_config')
    else:
        raise ValueError("Invalid statistic type for `record_binned_vectors`.")

    binned_vectors = record_binned_vectors(
        test_binning, paramset=test_paramset, boxsize=boxsize, ngrid=ngrid
    )
    sorting = np.lexsort((
        binned_vectors['vecz'],
        binned_vectors['vecy'],
        binned_vectors['vecx'],
        binned_vectors['index'],
    ))

    binned_vectors_sorted = binned_vectors[sorting]

    assert len(binned_vectors_sorted) == count, (
        "Binned vectors have incorrect count."
    )
    for ifield, name in enumerate(binned_vectors_sorted.dtype.names):
        assert np.allclose(
            binned_vectors_sorted[-1][ifield], last_entry[ifield]
        ), (
            "The last entry of binned vectors has "
            f"incorrect '{name}' field value."
        )
