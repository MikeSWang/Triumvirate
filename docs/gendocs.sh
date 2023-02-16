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

SPHINX_SOURCE_STATIC_DIR=./source/_static/
SPHINX_BUILD_HTML_DIR=./build/sphinx/html/
DOXYGEN_BUILD_HTML_DIR=./build/doxygen/html/

APIDOC_DOXY_DIR=${SPHINX_SOURCE_STATIC_DIR}/apiref_doxy/


# -- Build Docs ----------------------------------------------------------

# Clean up.
make clean

# Build Doxygen docs.
doxygen ${DOXY_CONF_FILE}

# Bridge Doxygen and Sphinx docs.
mkdir -p ${APIDOC_DOXY_DIR}
cp -r ${DOXYGEN_BUILD_HTML_DIR}/** ${APIDOC_DOXY_DIR}

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


# -- Organise Docs -------------------------------------------------------

HTML_PUBLIC_DIR=./build/public_html/
IMG_MIRROR_DIR_STEM=_static/_static/apiref_doxy/docs/source/_static/

if [[ "${READTHEDOCS}" != "True" ]]; then
    rm ./source/Doxyfile
    mkdir -p ${HTML_PUBLIC_DIR}
    cp -r ${SPHINX_BUILD_HTML_DIR}/** ${HTML_PUBLIC_DIR}
    mkdir -p ${HTML_PUBLIC_DIR}/${IMG_MIRROR_DIR_STEM}
    find ${SPHINX_SOURCE_STATIC_DIR} -maxdepth 1 -type f \
        -exec cp {} ${HTML_PUBLIC_DIR}/${IMG_MIRROR_DIR_STEM} \;
fi

# Return to original directory.
cd -
