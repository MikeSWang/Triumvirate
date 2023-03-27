#!/usr/bin/env bash
#
# @file autorel_testpypi.sh
# @author Mike S Wang
# @brief Upload release to test index repositories.
#

rm -rf dist/
python -m pip install --upgrade wheel build auditwheel twine
python -m build --sdist --wheel --outdir dist/ .
auditwheel -v repair dist/Triumvirate-*.whl
python -m twine upload --verbose --repository testpypi dist/*
