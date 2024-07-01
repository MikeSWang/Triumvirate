#!/usr/bin/env bash
#
# @file autorel_testpypi.sh
# @author Mike S Wang
# @brief Build and upload release distributions to the TestPyPI index.
#
# Require docker to run `cibuildwheel`.
#

# Parse options.
build_wheel=false
upload_flag=false
while getopts :wu opt; do
    case $opt in
        w) build_wheel=true;;
        u) upload_flag=true;;
        :) echo "Missing argument for option -$OPTARG"; exit 1;;
       \?) echo "Unknown option -$OPTARG"; exit 1;;
    esac
done
shift "$((OPTIND-1))"

# Set directories and paths.
DIST_DIR=dist/

# Clean distribution directory.
rm -rf "${DIST_DIR}"

# Install/upgrade distribution tools.
# python -m pip install --upgrade build cibuildwheel twine

# Build source-distribution.
python -m build --sdist --outdir "${DIST_DIR}" .

# Optionally build built-distribution.
if ${build_wheel}; then
    python -m cibuildwheel --platform linux --output-dir "${DIST_DIR}" .
fi

# Check
python -m twine check --strict "${DIST_DIR}"/*.tar.gz "${DIST_DIR}"/*.whl

# Optionally upload to [Test]PyPI.
if ${upload_flag}; then
    python -m twine upload --repository testpypi --verbose "${DIST_DIR}"/*
fi
