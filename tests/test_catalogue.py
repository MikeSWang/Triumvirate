"""Test :mod:`~triumvirate.catalogue`.

"""
import warnings

import numpy as np
import pytest

from triumvirate.catalogue import (
    DefaultValueWarning,
    ParticleCatalogue,
)

CONST_LEN = 1
CONST_COORDS = np.ones(CONST_LEN)
CONST_NZ = 1.e-3

SAMPLE_LEN = 1
SAMPLE_COORDS = {
    axis: [coord,] * SAMPLE_LEN  # noqa: E231
    for axis, coord in zip(['x', 'y', 'z'], [-1., 0., 1.])
}
SAMPLE_NZ = 1.e-3
SAMPLE_WSYS = 1.1
SAMPLE_WFKP = 0.9


@pytest.fixture
def minimal_catalogue():
    x = y = z = CONST_COORDS
    nz = None
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message=".*'nz' field is None.*")
        return ParticleCatalogue(x, y, z, nz=nz)


@pytest.fixture
def minimal_catalogue_paired():
    x = y = CONST_COORDS
    z = CONST_COORDS - 0.5
    nz = CONST_NZ
    return ParticleCatalogue(x, y, z, nz=nz)


@pytest.mark.parametrize(
    "x, y, z, nz, ws, wc",
    [
        (np.arange(5), -np.arange(5), np.arange(5) / 2., 1.e-3, 1.1, 0.9),
        (np.arange(5), -np.arange(5), np.arange(5) / 2., None, None, None),
    ]
)
def test_ParticleCatalogue___init__(x, y, z, nz, ws, wc):

    weight_kwargs = {}
    if ws is not None:
        weight_kwargs.update(ws=ws)
    if wc is not None:
        weight_kwargs.update(wc=wc)

    if nz is None:
        # Check warning is raised for missing the 'nz' field.
        with pytest.warns(DefaultValueWarning):
            ParticleCatalogue(x, y, z, nz=nz, **weight_kwargs)
    else:
        ParticleCatalogue(x, y, z, nz=nz, **weight_kwargs)


@pytest.mark.parametrize(
    "filename, reader, format, names, name_mapping",
    [
        (
            "sample_catalogue.txt",
            'astropy',
            'ascii.no_header',
            ['x', 'y', 'z', 'nz', 'ws', 'wc'],
            None,
        ),
        (
            "sample_catalogue.txt",
            'nbodykit',
            'text',
            ['x', 'y', 'z', 'nz', 'ws', 'wc'],
            None,
        ),
        (
            "sample_catalogue.fits",
            'astropy',
            'fits',
            None,
            {
                'x': 'X', 'y': 'Y', 'z': 'Z',
                'nz': 'NZ', 'ws': 'WEIGHT_SYS', 'wc': 'WEIGHT_FKP'
            },
        ),
        (
            "sample_catalogue.dat",
            'nbodykit',
            'binary',
            [
                ('X', np.float64), ('Y', np.float64), ('Z', np.float64),
                ('NZ', np.float64),
                ('WEIGHT_SYS', np.float64), ('WEIGHT_FKP', np.float64),
            ],
            {
                'x': 'X', 'y': 'Y', 'z': 'Z',
                'nz': 'NZ', 'ws': 'WEIGHT_SYS', 'wc': 'WEIGHT_FKP'
            },
        ),
        (
            "sample_catalogue.h5",
            'nbodykit',
            'hdf5',
            None,
            {
                'x': 'X', 'y': 'Y', 'z': 'Z',
                'nz': 'NZ', 'ws': 'WEIGHT_SYS', 'wc': 'WEIGHT_FKP'
            },
        ),
    ]
)
def test_ParticleCatalogue_read_from_file(filename, reader, format,
                                          names, name_mapping, test_ctlg_dir):

    if reader == 'nbodykit':
        import importlib
        try:
            importlib.import_module(reader)
        except ImportError:
            return

    catalogue = ParticleCatalogue.read_from_file(
        test_ctlg_dir/filename,
        reader=reader, names=names, format=format, name_mapping=name_mapping
    )

    for axis in ['x', 'y', 'z']:
        assert catalogue[axis] == SAMPLE_COORDS[axis], \
            "Loaded catalogue particle coordinates do not match."
    assert catalogue['nz'] == np.ones(CONST_LEN) * SAMPLE_NZ, \
        "Loaded catalogue particle number density column does not match."
    assert catalogue['ws'] == np.ones(CONST_LEN) * SAMPLE_WSYS, \
        "Loaded catalogue particle sample weight column does not match."
    assert catalogue['wc'] == np.ones(CONST_LEN) * SAMPLE_WFKP, \
        "Loaded catalogue particle clustering weight column does not match."


def test_ParticleCatalogue___str__(minimal_catalogue):
    assert 'extdata' in str(minimal_catalogue), \
        "Catalogue string representation has incorrect source."


def test_ParticleCatalogue___len__(minimal_catalogue):
    assert len(minimal_catalogue) == CONST_LEN, \
        "Catalogue length is incorrect."


def test_ParticleCatalogue___getitem__(minimal_catalogue):
    assert np.allclose(minimal_catalogue['x'], CONST_COORDS), \
        "Catalogue column not returned correctly."


def test_ParticleCatalogue_compute_los(minimal_catalogue):
    _los = np.tile(CONST_COORDS, (CONST_LEN, 3)) / np.linalg.norm(
        np.tile(CONST_COORDS, (CONST_LEN, 3)), axis=-1
    )
    assert np.allclose(minimal_catalogue.compute_los(), _los), \
        "Catalogue lines of sight not computed correctly."


@pytest.mark.parametrize(
    "volume, boxsize, nz",
    [
        (1.e3, None, CONST_NZ),
        (None, 10., CONST_NZ),
        (None, [10.,]*3, CONST_NZ),  # noqa: E231
    ]
)
def test_ParticleCatalogue_compute_mean_density(
    volume, boxsize, nz, minimal_catalogue
):
    minimal_catalogue.compute_mean_density(volume=volume, boxsize=boxsize)
    assert all(np.isclose(minimal_catalogue['nz'], nz)), \
        "Catalogue in-box mean density not computed correctly."


@pytest.mark.parametrize(
    "boxsize, ref_catalogue, extents",
    [
        (
            1000.,
            None,
            dict(zip(['x', 'y', 'z'], [[500.,]*2]*3)),  # noqa: E231
        ),
        (
            1000.,
            'minimal_catalogue_paired',
            dict(
                zip(
                    ['x', 'y', 'z'], [[500.,]*2]*2 + [[500.5,]*2]  # noqa: E231
                )
            ),
        ),
    ]
)
def test_ParticleCatalogue_centre(boxsize, ref_catalogue, extents,
                                  minimal_catalogue, request):

    if ref_catalogue is not None:
        ref_catalogue = request.getfixturevalue(ref_catalogue)

    minimal_catalogue.centre(boxsize, catalogue_ref=ref_catalogue)
    for axis in ['x', 'y', 'z']:
        assert np.allclose(minimal_catalogue.bounds[axis], extents[axis]), \
            "Catalogue particles are not centred correctly in the box."


@pytest.mark.parametrize(
    "boxsize, ngrid, boxsize_pad, ngrid_pad, ref_catalogue, extents",
    [
        (
            1000.,
            None,
            0.05,
            None,
            None,
            dict(zip(['x', 'y', 'z'], [[50.,]*2]*3)),  # noqa: E231
        ),
        (
            1000.,
            None,
            0.05,
            None,
            'minimal_catalogue_paired',
            dict(
                zip(
                    ['x', 'y', 'z'], [[50.,]*2]*2 + [[50.5,]*2]  # noqa: E231
                )
            ),
        ),
        (
            1000.,
            10,
            None,
            0.5,
            None,
            dict(zip(['x', 'y', 'z'], [[50.,]*2]*3)),  # noqa: E231
        ),
    ]
)
def test_ParticleCatalogue_pad(boxsize, ngrid, boxsize_pad, ngrid_pad,
                               ref_catalogue, extents,
                               minimal_catalogue, request):

    if ref_catalogue is not None:
        ref_catalogue = request.getfixturevalue(ref_catalogue)

    minimal_catalogue.pad(
        boxsize, ngrid=ngrid, boxsize_pad=boxsize_pad, ngrid_pad=ngrid_pad,
        catalogue_ref=ref_catalogue
    )

    for axis in ['x', 'y', 'z']:
        assert np.allclose(minimal_catalogue.bounds[axis], extents[axis]), \
            "Catalogue particles are not padded correctly in the box."


@pytest.mark.parametrize(
    "boxsize, extents",
    [
        (
            0.3,
            dict(
                zip(
                    ['x', 'y', 'z'],
                    [[0.3/2.,]*2]*3  # noqa: E231
                )
            ),
        ),
    ]
)
def test_ParticleCatalogue_periodise(boxsize, extents, minimal_catalogue):
    minimal_catalogue.periodise(boxsize)
    for axis in ['x', 'y', 'z']:
        assert np.allclose(minimal_catalogue.bounds[axis], extents[axis]), \
            "Catalogue particle coordinates do not satisfy periodicity."


@pytest.mark.parametrize(
    "origin, extents",
    [
        (
            [1., 1., 1.],
            dict(zip(['x', 'y', 'z'], [[0.,]*2]*3)),  # noqa: E231
        ),
    ]
)
def test_ParticleCatalogue_offset_coords(origin, extents, minimal_catalogue):
    minimal_catalogue.offset_coords(origin)
    for axis in ['x', 'y', 'z']:
        assert np.allclose(minimal_catalogue.bounds[axis], extents[axis]), \
            "Catalogue particle coordinates are not offset correctly."


@pytest.mark.parametrize(
    "ref_catalogue",
    [
        'minimal_catalogue_paired',
        None,
    ]
)
def test_ParticleCatalogue_write_attrs_as_header(ref_catalogue,
                                                 minimal_catalogue, request):

    if ref_catalogue is not None:
        ref_catalogue = request.getfixturevalue(ref_catalogue)

    header = minimal_catalogue.write_attrs_as_header(
        catalogue_ref=ref_catalogue
    )

    if ref_catalogue is None:
        assert (
            "Catalogue source: {}".format(minimal_catalogue._source)
        ) in header
        assert (
            "Catalogue size: ntotal = {:d}, wtotal = {:.3f}, wstotal = {:.3f}"
            .format(
                minimal_catalogue.ntotal,
                minimal_catalogue.wtotal,
                minimal_catalogue.wstotal,
            )
        ) in header
        assert (
            "Catalogue particle extents: "
            "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
            .format(
                *minimal_catalogue.bounds['x'],
                *minimal_catalogue.bounds['y'],
                *minimal_catalogue.bounds['z'],
            )
        ) in header
    else:
        for _source_type, _catalogue in zip(
            ['Data', 'Random'], [minimal_catalogue, ref_catalogue]
        ):
            assert (
                "{} catalogue source: {}"
                .format(_source_type, _catalogue._source)
            ) in header
            assert (
                "{} catalogue size: "
                "ntotal = {:d}, wtotal = {:.3f}, wstotal = {:.3f}"
                .format(
                    _source_type,
                    _catalogue.ntotal, _catalogue.wtotal, _catalogue.wstotal,
                )
            ) in header
            assert (
                "{}-source particle extents: "
                "([{:.3f}, {:.3f}], [{:.3f}, {:.3f}], [{:.3f}, {:.3f}])"
                .format(
                    _source_type,
                    *_catalogue.bounds['x'],
                    *_catalogue.bounds['y'],
                    *_catalogue.bounds['z'],
                )
            ) in header
