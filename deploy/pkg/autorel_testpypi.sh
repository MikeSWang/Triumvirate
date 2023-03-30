#!/usr/bin/env bash
#
# @file autorel_pypi.sh
#
# Require docker to run `cibuildwheel`.
#
# @author Mike S Wang
# @brief Build and upload release distributions to the PyPI index.
#

# Set directories and paths.
DIST_DIR=dist/

# Clean distribution directory.
rm -rf ${DIST_DIR}

# Install/upgrade distribution tools.
python -m pip install --upgrade build cibuildwheel twine

# Build both source and built distributions.
python -m build --sdist --outdir ${DIST_DIR} .
python -m cibuildwheel --platform linux --output-dir ${DIST_DIR} .

# Check and then upload to [Test]PyPI.
python -m twine check ${DIST_DIR}/* && \
python -m twine upload --repository testpypi --verbose ${DIST_DIR}/*
