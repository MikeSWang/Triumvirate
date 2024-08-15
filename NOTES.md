# Developer's Notes

## Pre-release checklist

- [ ] Generate version number and change log for new release
      (by running [`prerelease_ops.sh`](deploy/pkg/prerelease_ops.sh))
- [ ] Bump fallback version number
      (by running [`autoinc_version.sh`](deploy/pkg/autoinc_version.sh)) in:
    - [`monitor.hpp`](src/triumvirate/include/monitor.hpp)
    - [`__init__.py`](src/triumvirate/__init__.py)
    - [Conda recipe](deploy/pkg/conda_recipe/meta.yaml)
    - [Conda recipe cross-platform](deploy/pkg/conda_recipe_xp/meta.yaml)
    - [RTD version switcher](docs/versions.json)
    - [RTD announcement banner](docs/source/conf.py)
- [ ] Check agreement amongst read-me pages:
    - [repo read-me](README.md)
    - [PyPI read-me](README.rst)
    - [ReadtheDocs homepage](docs/source/index.rst)
- [ ] Check agreement between installation instructions:
    - [repo read-me](README.md#Installation)
    - [PyPI read-me](README.rst#Installation)
    - [ReadtheDocs](docs/source/installation.rst)
- [ ] Check agreement between release notes:
    - [change-log](CHANGELOG.md)
    - [release history](docs/source/releases.rst)


## Code tags

- The following code tags are in use:
    - 'NOTE', 'STYLE';
    - 'CAVEAT';
    - 'RFE', 'FUTURE';
    - 'TODO'.

- The following code tags are available for use:
    - 'WARNING';
    - 'HACK';
    - 'FIXME', 'QUEST'.
