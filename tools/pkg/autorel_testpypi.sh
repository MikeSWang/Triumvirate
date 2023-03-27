#!/usr/bin/env bash
#
# @file autorel_testpypi.sh
# @author Mike S Wang
# @brief Upload release to test index repositories.
#

# Clean distribution directory.
rm -rf dist/

# Install distribution tools.
python -m pip install --upgrade build auditwheel twine

# Build both 'sdist' and 'bdist'.
python -m build --sdist --wheel --outdir dist/ .

# Audit built wheels.
auditwheel -v repair dist/Triumvirate-*.whl

# # Upload 'sdist' only for now.
# python -m twine upload --verbose --repository testpypi dist/*.tar.gz
