#!/usr/bin/env bash
#
# @file autorel_pypi.sh
# @author Mike S Wang
# @brief Build and upload release distributions to the PyPI index.
#

# Set directories and paths.
DIST_DIR=dist/

# Clean distribution directory.
rm -rf ${DIST_DIR}

# Install distribution tools.
python -m pip install --upgrade build auditwheel twine

# Build both source and built distributions.
python -m build --sdist --wheel --outdir ${DIST_DIR} .

# Repair built wheels.
auditwheel -v repair ${DIST_DIR}/Triumvirate-*.whl

# Upload to PyPI.
python -m twine upload --verbose ${DIST_DIR}/*
