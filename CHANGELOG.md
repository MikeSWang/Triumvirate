# Change Log

## v0.5 (expected in 2024)

### Bug fixes

### Features

- Add public API for window convolution.

### Improvements

### Maintenance

### Documentation

- Add documentation for window convolution.

### Miscellaneous


## v0.4.6 (2024-06-24)

### Bug fixes

- Fix asymmetric 2-d three-point clustering statistics.

- Fix various legacy issues in the non-public methods related to
  Hankel-like transforms.

- Fix FFTW-related memory leaks.

### Features

- Add public API for Hankel-like transforms using the FFTLog algorithm
  with enhanced extrapolation options.

- Expose installation validation as a public function.

- Add non-public API for window convolution (as part of a future release).

### Improvements

- Add support for FFTW planner flags and wisdom.

- Refactor assignment and aliasing compensation.

- Refactor mesh field (re-)initialisation.

- Refactor gamma function computations.

- Refactor array operation checks.

- Refactor FFTW plans.

- Add logs to ``trv::MeshField`` and ``trv::FieldStats`` operations.

- Add tracking of (I)FFTs and 3-d grids.

- Promote `trv::ParameterSet::nmesh` from type `int` to `long long` and
  refactor related methods and variables.

### Maintenance

- Require C++17 standard.
- Upgrade build against NumPy 2 and require Python 3.10+
  ([gh-52](https://github.com/MikeSWang/Triumvirate/issues/52)).
- Refactor logger message emission.
- Enhance build recipes.
- Update syntax and fix typos.

### Documentation

- Add version context to documentation.
- Add and update details to/in various documentation components.

### Miscellaneous

Relicense under/clarify GPL-3.0-or-later in lieu of GPL-3.0 where applicable.


## v0.3.0 (2023-10-04)

### Bug fixes

Fix reused bin statistics when duplicate recording is designed to be avoided
in [``threept.cpp``](src/triumvirate/src/threept.cpp).

### Features

- Add the functionality to record binning details for wavevector modes
  and separation pairs from a mesh grid with the option to save the results
  to a file in C++.

- Add more forms of three-point statistics including the off-diagonal and
  full form (the original 'full' form is renamed to 'row')
  ([gh-22](https://github.com/MikeSWang/Triumvirate/issues/22)).

### Improvements

- Implement further parallelisation and refactoring resulting in
  significant speed-up.

- Enhance parameter validation consistency.

### Maintenance

- Update syntax and fix typos.
- Add and update test data after API changes.

### Documentation

- Add and update details to/in various documentation components.
- Update and rerun tutorial notebooks after API changes.


## v0.2.2 (2023-07-04)

### Improvements

Enhance build process and update syntax.

### Documentation

Update installation guide after build enhancement.


## v0.2.1 (2023-06-20)

### Bug fixes

Fix parity factor in three-point correlation functions in
[``threept.cpp``](src/triumvirate/src/threept.cpp).

### Maintenance

- Remove reality-condition division in mode/pair counts for generality.
- Update test data.
- Update syntax.

### Documentation

Rerun tutorial notebooks.


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
