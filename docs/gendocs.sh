#!/usr/bin/env bash

DOCS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

cd ${DOCS_DIR}
if [ ! -d api ]; then mkdir -p ./source/api; fi

make clean
sphinx-apidoc -efEMT -d 1 -t ./source/_templates/ -o ./source/apidoc ../triumvirate ../**/tests/*
rm ./source/apidoc/triumvirate.rst
make html
cd ..
