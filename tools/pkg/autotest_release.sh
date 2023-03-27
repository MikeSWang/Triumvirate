#!/usr/bin/env bash
#
# @file autotest_vers.sh
# @author Mike S Wang
# @brief Upload release to test index repositories.
#

rm -rf dist/
python -m pip install --upgrade wheel build auditwheel twine
python -m build
auditwheel -v repair dist/Triumvirate-*.whl
python -m twine upload --verbose --repository testpypi dist/*
