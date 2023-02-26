#!/usr/bin/env bash
#
# @file autogen_docs.sh
# @author Mike S Wang
# @brief Automate documentation generation locally and on `RTD <rtfd.io>`_.
# @arg Exclusion case, {"excl_doxy", "excl_sphinx"}.
#

# Parse command-line arguments.
excl=${1#"excl_"}  # {"doxy", "sphinx"}


# -- Directories & Paths -------------------------------------------------

PROJ_NAME=triumvirate

# Set the ``docs`` directory (the working directory).
DOCS_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# Set Doxygen configuration file.
DOXYFILE=./source/Doxyfile.conf
DOXYFILE_SPHINX=./source/Doxyfile

# Set Sphinx subdirectories.
STAT_DIR=./source/_static/
TMPL_DIR=./source/_templates/

APIREF_DOXY_DIRNAME=apiref_doxy
APIREF_DOXY_DIR=./source/${APIREF_DOXY_DIRNAME}/
APIDOC_PY_DIR=./source/apidoc_py/
APIDOC_CPP_DIR_STR_ESC=".\/source\/apidoc_cpp\/"

# Set HTML build directories.
BUILD_HTML_DOXY_DIR=./build/doxygen/html/
BUILD_HTML_SPHINX_DIR=./build/sphinx/html/

# Set resource origin directory.
RESOURCE_DIR=../${PROJ_NAME}/resources/

# Set exclusion path patterns.
EXCL_DIRS="../${PROJ_NAME} ../**/tests/*"
RM_FILES="${APIDOC_PY_DIR}/${PROJ_NAME}.rst"


# -- Build Docs ----------------------------------------------------------

# @brief Get the project release number.
#
get_version_release () {
    echo $(python -c "import ${PROJ_NAME}; print(${PROJ_NAME}.__version__)")
}

# @brief Replace strings in file.
#
# @arg File.
# @arg String to be replaced.
# @arg String replacement.
#
replace_in_file () {
    sed -i "s/${2}/${3}/g" ${1}
}

# @brief Replace key--value pairs in Doxyfile.conf.
#
# @arg Replacement string.
# @globals DOXYFILE
# @locals str_confline, str_confvar
#
replace_in_doxyfile () {
    str_confline=$1
    str_confvar=$(printf $1 | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s/${str_confvar} .*=.*/${str_confline}/g" ${DOXYFILE}
}

# @brief Replace key--value pairs in Doxyfile[.conf].
#
# @arg Replacement string.
# @globals DOXYFILE_SPHINX
# @locals str_confline, str_confvar
#
recycle_doxyfile () {
    str_confline=$1
    str_confvar=$(printf $1 | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s/${str_confvar} .*=.*/${str_confline}/g" ${DOXYFILE_SPHINX}
}

# Change to working directory.
cd ${DOCS_DIR}

# Clean up.
make clean excl=$excl

# Pre-configure Doxygen-related docs.
replace_in_doxyfile "PROJECT_NUMBER         = $(get_version_release)"

# HACK: Ensure RTD Doxygen backward compatibility.
if [[ "${READTHEDOCS}" == "True" ]]; then
    MATHJAX_RELPATH="https:\/\/cdn.jsdelivr.net\/npm\/mathjax@2"
    replace_in_file ./source/_themes/doxygen-header.html "\$darkmode"
    replace_in_doxyfile "MATHJAX_RELPATH        = ${MATHJAX_RELPATH}"
fi

# Prepare static assets.
cp -r ${RESOURCE_DIR}/params_template.* ${STAT_DIR}

# Build Doxygen docs.
if [[ "$excl" != 'doxy' ]]; then
    doxygen ${DOXYFILE}
fi

# Build Sphinx docs.
if [[ "$excl" != 'sphinx' ]]; then
    # Bridge Doxygen and Sphinx docs.
    mkdir -p ${APIREF_DOXY_DIR}
    cp -r ${BUILD_HTML_DOXY_DIR}/** ${APIREF_DOXY_DIR}

    # Configure Doxygen for Breathe+Exhale.
    cp ${DOXYFILE} ${DOXYFILE_SPHINX}

    recycle_doxyfile "OUTPUT_DIRECTORY       = ${APIDOC_CPP_DIR_STR_ESC}"
    recycle_doxyfile "GENERATE_HTML          = NO"
    recycle_doxyfile "GENERATE_XML           = YES"
    recycle_doxyfile "HAVE_DOT               = NO"
    recycle_doxyfile "HTML_HEADER            ="
    recycle_doxyfile "USE_MDFILE_AS_MAINPAGE ="
    sed -i "s/..\/README.md//g" ${DOXYFILE_SPHINX}
    sed -i "s/.\/source/..\/source/g" ${DOXYFILE_SPHINX}
    sed -i "s/..\/${PROJ_NAME}/..\/..\/${PROJ_NAME}/g" ${DOXYFILE_SPHINX}

    # Build docs with Breathe+Exhale.
    sphinx-apidoc -efEMT -d 1 -t ${TMPL_DIR} -o ${APIDOC_PY_DIR} ${EXCL_DIRS}

    rm ${RM_FILES}

    make html

    # Clean up.
    if [[ "${READTHEDOCS}" != "True" ]]; then rm ${DOXYFILE_SPHINX}; fi
fi


# -- Organise Docs -------------------------------------------------------

if [[ "${READTHEDOCS}" != "True" ]]; then
    # Set directories.
    HTML_PUBLIC_DIR=./build/public_html/
    STAT_APIREF_DOXY_SUBDIR=_static/${APIREF_DOXY_DIRNAME}/
    STAT_MIRROR_SUBDIR=${APIREF_DOXY_DIRNAME}/docs/source/_static/

    # Move Sphinx-build HTML to public HTML directory.
    mkdir -p ${HTML_PUBLIC_DIR}
    if [[ -d ${BUILD_HTML_SPHINX_DIR} ]]; then
        cp -r ${BUILD_HTML_SPHINX_DIR}/** ${HTML_PUBLIC_DIR}
    fi

    # Promote Doxygen-build HTML to top-level subdirectory
    # under the public HTML directory.
    if [[ -d ${HTML_PUBLIC_DIR}/${STAT_APIREF_DOXY_SUBDIR} ]]; then
        mv ${HTML_PUBLIC_DIR}/${STAT_APIREF_DOXY_SUBDIR} ${HTML_PUBLIC_DIR}/
    fi

    # Move static assets to a mirrored subdirectory
    # in the public HTML directory for Doxygen 'mainpage'.
    mkdir -p ${HTML_PUBLIC_DIR}/${STAT_MIRROR_SUBDIR}
    find ${STAT_DIR} -maxdepth 1 -type f \
        -exec cp {} ${HTML_PUBLIC_DIR}/${STAT_MIRROR_SUBDIR} \;

    # Redirect HTML listings to the index page.
    echo -e "RewriteEngine On\nRewriteRule (.*)\$ index.html" \
        > ${HTML_PUBLIC_DIR}.htaccess
fi

# Return to original directory.
cd -
