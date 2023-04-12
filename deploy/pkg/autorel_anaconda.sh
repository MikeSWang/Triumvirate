#!/usr/bin/env bash
#
# @file autorel_anaconda.sh
#
# @author Mike S Wang
# @brief Build and verify release distributions for the Anaconda repository.
# @arg Non-native architecture case, {"", "arm64"}.
#

# Parse command-line arguments.
if [[ -z "${1}" ]]; then
  arch_suffix=
else
  arch_suffix=_${1}  # non-native architecture case: {"arm64"}
fi

# Set directories and paths.
RECIPE_DIR=deploy/pkg/conda_recipe${arch_suffix}
DIST_DIR=dist/

# Clean distribution directory.
rm -rf ${DIST_DIR}
if [[ ! -d ${DIST_DIR} ]]; then mkdir -p ${DIST_DIR}; fi

# Install/upgrade distribution tools.
# conda update -y conda-build conda-verify conda-package-handling

# Build and verify built-distribution.
conda build purge
conda build --strict-verify --no-anaconda-upload ${RECIPE_DIR} \
  --output-folder ${DIST_DIR} \
  --variants "{'python': ['3.8', '3.9', '3.10']}"

# Transmute compression formats.
find ${DIST_DIR} -name "*.tar.bz2" \
  -exec cph transmute {} .conda --out-folder ${DIST_DIR} \;
