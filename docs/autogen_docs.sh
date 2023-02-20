#!/usr/bin/env bash
#
# @file autogen_docs.sh
# @author Mike S Wang
# @brief Automate documentation generation locally and on `RTD <rtfd.io>`_.
# @args Exclusion case, {"excl_doxy", "excl_sphinx"}.
#

# Parse command-line arguments.
excl=${1#"excl_"}  # {"doxy", "sphinx"}


# -- Directories & Paths -------------------------------------------------

# Change to the ``docs`` directory.
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

# @func replace_in_file
#
# Replace strings in file.
#
# @args File, string to be replaced, string replacement.
#
replace_in_file () {
    sed -i "s/${2}/${3}/g" $1
}

# @func recycle_doxyfile
#
# Replace key--value pairs in Doxyfile[.conf].
#
# @args Replacement string.
# @globals DOXYFILE_FOR_EXHALE
# @locals str_confline, str_confvar
#
recycle_doxyfile () {
    str_confline=$1
    str_confvar=$(printf $1 | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s/${str_confvar} .*=.*/${str_confline}/g" ${DOXYFILE_FOR_EXHALE}
}

# HACK: Ensure RTD Doxygen backward compatibility.
if [[ "${READTHEDOCS}" == "True" ]]; then
    replace_in_file ./source/_themes/doxygen-header.html "\$darkmode"
    replace_in_file ./source/_themes/doxygen-header.html \
        "top: 0; right: 0;" "bottom: 0; left: 0; transform: scale(-1, -1);"
    replace_in_file ${DOXY_CONF_FILE} "doxygen-awesome-sidebar-only.css"
    replace_in_file ${DOXY_CONF_FILE} \
        "= https://cdn.jsdelivr.net/npm/mathjax@3" \
        "= https://cdn.jsdelivr.net/npm/mathjax@2"
fi

# Clean up.
make clean excl=$excl

# Build Doxygen docs.
if [[ "$excl" != 'doxy' ]]; then
    doxygen ${DOXY_CONF_FILE}
fi

# Build Sphinx docs.
if [[ "$excl" != 'sphinx' ]]; then
    # Bridge Doxygen and Sphinx docs.
    mkdir -p ${APIREF_DOXY_DIR_FOR_SPHINX}
    cp -r ${DOXYGEN_BUILD_HTML_DIR}/** ${APIREF_DOXY_DIR_FOR_SPHINX}

    # Configure Doxygen for Breathe+Exhale.
    cp ${DOXY_CONF_FILE} ${DOXYFILE_FOR_EXHALE}

    recycle_doxyfile "OUTPUT_DIRECTORY       = ${APIDOC_CPP_DIR_FOR_EXHALE}"
    recycle_doxyfile "GENERATE_HTML          = NO"
    recycle_doxyfile "GENERATE_XML           = YES"
    recycle_doxyfile "HAVE_DOT               = NO"
    recycle_doxyfile "HTML_HEADER            ="
    recycle_doxyfile "USE_MDFILE_AS_MAINPAGE ="
    sed -i "s/..\/README.md//g" ${DOXYFILE_FOR_EXHALE}
    sed -i "s/.\/source/..\/source/g" ${DOXYFILE_FOR_EXHALE}
    sed -i "s/..\/triumvirate/..\/..\/triumvirate/g" ${DOXYFILE_FOR_EXHALE}

    # Build Sphinx docs with Doxygen+Breathe+Exhale.
    sphinx-apidoc -efEMT -d 1 -t ${TMPL_DIR} -o ${APIDOC_PY_DIR} ${EXCL_DIRS}

    rm ${RM_FILES}

    make html

    # Clean up.
    if [[ "${READTHEDOCS}" != "True" ]]; then rm ${DOXYFILE_FOR_EXHALE}; fi
fi


# -- Organise Docs -------------------------------------------------------

if [[ "${READTHEDOCS}" != "True" ]]; then
    # Set directories.
    HTML_PUBLIC_DIR=./build/public_html/
    STAT_APIREF_DOXY_SUBDIR=_static/${APIREF_DOXY_DIRNAME}/
    IMG_MIRROR_SUBDIR=${APIREF_DOXY_DIRNAME}/docs/source/_static/

    # Move Sphinx-build HTML to public HTML directory.
    mkdir -p ${HTML_PUBLIC_DIR}
    if [[ -d ${SPHINX_BUILD_HTML_DIR} ]]; then
        cp -r ${SPHINX_BUILD_HTML_DIR}/** ${HTML_PUBLIC_DIR}
    fi

    # Promote Doxygen-build HTML to top-level subdirectory
    # under the public HTML directory.
    if [[ -d ${HTML_PUBLIC_DIR}/${STAT_APIREF_DOXY_SUBDIR} ]]; then
        mv ${HTML_PUBLIC_DIR}/${STAT_APIREF_DOXY_SUBDIR} ${HTML_PUBLIC_DIR}/
    fi

    # Move static assets to a mirrored subdirectory
    # in the public HTML directory for Doxygen 'mainpage'.
    mkdir -p ${HTML_PUBLIC_DIR}/${IMG_MIRROR_SUBDIR}
    find ${STAT_DIR} -maxdepth 1 -type f \
        -exec cp {} ${HTML_PUBLIC_DIR}/${IMG_MIRROR_SUBDIR} \;

    # Redirect listings to index pages.
    echo -e "RewriteEngine On\nRewriteRule (.*)\$ index.html" \
        > ${HTML_PUBLIC_DIR}.htaccess
fi

# Return to original directory.
cd -
