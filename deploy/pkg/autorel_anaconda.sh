#!/usr/bin/env bash
#
# @file autorel_anaconda.sh
#
# @author Mike S Wang
# @brief Build and verify release distributions for the Anaconda repository.
#

# Set directories and paths.
RECIPE_DIR=deploy/pkg/conda_recipe
DIST_DIR=dist/

# Clean distribution directory.
rm -rf ${DIST_DIR}
if [[ ! -d ${DIST_DIR} ]]; then mkdir -p ${DIST_DIR}; fi

# Install/upgrade distribution tools.
# conda update conda-build

# Build built-distribution.
conda build purge
conda build ${RECIPE_DIR} \
  --output-folder ${DIST_DIR} \
  --variants "{'python': ['3.8', '3.9', '3.10']}"

# Check.
conda verify ${DIST_DIR}
