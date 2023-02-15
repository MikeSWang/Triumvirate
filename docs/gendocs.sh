#!/usr/bin/env bash

# -- Directories & Paths -------------------------------------------------

# Change to the docs directory.
DOCS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

cd ${DOCS_DIR}

# Set API doc directories and paths.
DOXY_CONF_FILE=./source/Doxyfile.conf

APIDOC_PY_DIR=./source/apidoc_py
TMPL_DIR=./source/_templates/
EXCL_DIRS="../triumvirate ../**/tests/*"

RM_FILES="${APIDOC_PY_DIR}/triumvirate.rst"


# -- Build Docs ----------------------------------------------------------

# Clean up.
make clean

# Build Doxygen docs.
doxygen ${DOXY_CONF_FILE}

# Build Sphinx docs with Doxygen+Breathe+Exhale.
cp ${DOXY_CONF_FILE} ./source/Doxyfile

SPHINX_DOXYFILE_OUTPUT_DIR=".\/source\/apidoc_cpp\/"

sed -i\
    "s/OUTPUT_DIRECTORY .*=.*/OUTPUT_DIRECTORY       = ${SPHINX_DOXYFILE_OUTPUT_DIR}/g"\
    ./source/Doxyfile
sed -i\
    "s/GENERATE_HTML .*=.*/GENERATE_HTML          = NO/g"\
    ./source/Doxyfile
sed -i\
    "s/GENERATE_XML .*=.*/GENERATE_XML          = YES/g"\
    ./source/Doxyfile
sed -i\
    "s/HTML_HEADER .*=.*/HTML_HEADER            =/g"\
    ./source/Doxyfile
sed -i\
    "s/USE_MDFILE_AS_MAINPAGE .*=.*/USE_MDFILE_AS_MAINPAGE =/g"\
    ./source/Doxyfile
sed -i\
    "s/..\/README.md//g"\
    ./source/Doxyfile
sed -i\
    "s/.\/source/..\/source/g"\
    ./source/Doxyfile
sed -i\
    "s/..\/triumvirate/..\/..\/triumvirate/g"\
    ./source/Doxyfile
sphinx-apidoc -efEMT -d 1\
    -t ${TMPL_DIR} -o ${APIDOC_PY_DIR} ${EXCL_DIRS}

rm ${RM_FILES}
make html
rm ./source/Doxyfile

# Return to original directory.
cd -
