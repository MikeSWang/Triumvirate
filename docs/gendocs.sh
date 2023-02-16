#!/usr/bin/env bash

# -- Directories & Paths -------------------------------------------------

# Change to the docs directory.
DOCS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

cd ${DOCS_DIR}

# Set Doxygen configuration file.
DOXY_CONF_FILE=./source/Doxyfile.conf
DOXYFILE_FOR_EXHALE=./source/Doxyfile

# Sey API directories.
APIREF_DOXY_DIRNAME=apiref_doxy
APIREF_DOXY_DIR_FOR_SPHINX=./source/${APIREF_DOXY_DIRNAME}/
APIDOC_PY_DIR=./source/apidoc_py/
APIDOC_CPP_DIR_FOR_EXHALE=".\/source\/apidoc_cpp\/"

# Set source subdirectories.
STAT_DIR=./source/_static/
TMPL_DIR=./source/_templates/

# Set built HTML directories.
SPHINX_BUILD_HTML_DIR=./build/sphinx/html/
DOXYGEN_BUILD_HTML_DIR=./build/doxygen/html/

# Set exclusion path patterns.
EXCL_DIRS="../triumvirate ../**/tests/*"
RM_FILES="${APIDOC_PY_DIR}/triumvirate.rst"


# -- Build Docs ----------------------------------------------------------

# Pre-configure.
recycle_doxyfile () {
    str_confline=$1
    str_confvar=$(printf $1 | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s/${str_confvar} .*=.*/${str_confline}/g" ${DOXYFILE_FOR_EXHALE}
}

if [[ "${READTHEDOCS}" == "True" ]]; then
    sed -i "s/\$darkmode//g" ./docs/source/_themes/doxygen-header.html
fi

# Clean up.
make clean

# Build Doxygen docs.
doxygen ${DOXY_CONF_FILE}

# Bridge Doxygen and Sphinx docs.
mkdir -p ${APIREF_DOXY_DIR_FOR_SPHINX}
cp -r ${DOXYGEN_BUILD_HTML_DIR}/** ${APIREF_DOXY_DIR_FOR_SPHINX}

# Build Sphinx docs with Doxygen+Breathe+Exhale.
cp ${DOXY_CONF_FILE} ${DOXYFILE_FOR_EXHALE}

recycle_doxyfile "OUTPUT_DIRECTORY       = ${APIDOC_CPP_DIR_FOR_EXHALE}"
recycle_doxyfile "GENERATE_HTML          = NO"
recycle_doxyfile "GENERATE_XML           = YES"
recycle_doxyfile "HTML_HEADER            ="
recycle_doxyfile "USE_MDFILE_AS_MAINPAGE ="
sed -i "s/..\/README.md//g" ${DOXYFILE_FOR_EXHALE}
sed -i "s/.\/source/..\/source/g" ${DOXYFILE_FOR_EXHALE}
sed -i "s/..\/triumvirate/..\/..\/triumvirate/g" ${DOXYFILE_FOR_EXHALE}

sphinx-apidoc -efEMT -d 1\
    -t ${TMPL_DIR} -o ${APIDOC_PY_DIR} ${EXCL_DIRS}

rm ${RM_FILES}
make html

# Clean up.
if [[ "${READTHEDOCS}" != "True" ]]; then rm ${DOXYFILE_FOR_EXHALE}; fi


# -- Organise Docs -------------------------------------------------------

HTML_PUBLIC_DIR=./build/public_html/
IMG_MIRROR_SUBDIR=${APIREF_DOXY_DIRNAME}/docs/source/_static/

if [[ "${READTHEDOCS}" != "True" ]]; then
    # Move Sphinx-build HTML to public HTML directory.
    mkdir -p ${HTML_PUBLIC_DIR}
    cp -r ${SPHINX_BUILD_HTML_DIR}/** ${HTML_PUBLIC_DIR}
    # Promote Doxygen-build HTML to top-level subdirectory
    # under the public HTML directory.
    mv ${HTML_PUBLIC_DIR}/_static/${APIREF_DOXY_DIRNAME}/ ${HTML_PUBLIC_DIR}/
    # Move static assets to a mirrored subdirectory
    # in the public HTML directory for Doxygen mainpage.
    mkdir -p ${HTML_PUBLIC_DIR}/${IMG_MIRROR_SUBDIR}
    find ${STAT_DIR} -maxdepth 1 -type f \
        -exec cp {} ${HTML_PUBLIC_DIR}/${IMG_MIRROR_SUBDIR} \;
    echo -e "RewriteEngine On\nRewriteRule (.*)\$ index.html" \
        > ${HTML_PUBLIC_DIR}.htaccess
fi

# Return to original directory.
cd -
