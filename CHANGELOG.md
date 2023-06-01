# Change Log

## v0.2.0 (2023-06-01)

### Bug fixes

Fix updating of derived parameters ``npoint`` and ``space`` in ``ParameterSet``
in [``parameters.pyx``](src/triumvirate/parameters.pyx).

### Features

- Add ``pypower``-like normalisation for two-point clustering statistics with
  the new value option 'mesh-mixed' for the ``norm_convention`` parameter.

- Separate ``wtotal`` and ``wstotal`` attributes for ``ParticleCatalogue``.
  This is also reflected in logging and output file headers.

- Add particle position spans as a new member ``pos_span``
  in ``trv::ParticleCatalogue`` as this is used for ``pypower``-like
  normalisation calculations.

### Improvements

Enhance logging and build and packaging processes.

### Maintenance

Update syntax and fix typos.

### Documentation

Rerun tutorial notebooks.


## v0.1.2 (2023-04-12)

### Bug fixes

Fix parsing of build environmental variables in [``setup.py``](setup.py).

### Improvements

Upgrade build and packaging processes.

### Maintenance

Update syntax and fix typos.

### Documentation

Update installation guide based on the enhanced build process.


## v0.1.1 (2023-04-07)

Initial full public release ([major version 0](https://semver.org/#spec-item-4)):
publish to the PyPI index and Anaconda repository.

### Bug fixes

Sort loaded measurements files in [``application/tools/comb_data_vectors.py``](
application/tools/comb_data_vectors.py).

### Improvements

Upgrade build and packaging processes.

### Documentation

Update installation guide.


## v0.1.0 (2023-03-30)

Initial public release ([major version 0](https://semver.org/#spec-item-4)):
publish to the PyPI Index (Anaconda repository pending).
