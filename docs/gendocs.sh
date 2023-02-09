#!/usr/bin/env bash

# Change to the docs directory.
DOCS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

cd ${DOCS_DIR}

# Set API doc directories and paths.
APIDOC_DIR=./source/apidoc
TMPL_DIR=./source/_templates/
EXCL_DIRS="../triumvirate ../**/tests/*"

RM_FILES="${APIDOC_DIR}/triumvirate.rst"

# Build docs.
make clean

sphinx-apidoc -efEMT -d 1 -t ${TMPL_DIR} -o ${APIDOC_DIR} ${EXCL_DIRS}
rm ${RM_FILES}

make html

# Return to original directory.
cd -
