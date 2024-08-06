#!/usr/bin/env bash
#
# @file autogen_docs.sh
# @author Mike S Wang
# @brief Automate documentation generation locally and on `RTD <rtfd.io>`_.
# @arg Exclusion case, {'excl_doxy', 'excl_sphinx'}.
#

# Parse command-line arguments.
excl="${1#'excl_'}"  # exclusion case: {'doxy', 'sphinx'}


# -- Directories & Paths -------------------------------------------------

PKG_NAME=triumvirate

# Set the ``docs`` directory (the working directory).
THIS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
DOCS_DIR="${THIS_DIR}/../../docs"

# Set the project repository directory.
ROOT_DIR="${DOCS_DIR}/.."

# Set Doxygen configuration file.
DOXYFILE=./source/Doxyfile.conf
DOXYFILE_SPHINX=./source/Doxyfile

# Set Sphinx subdirectories.
STAT_DIR=./source/_static/
TMPL_DIR=./source/_templates/
THMS_DIR=./source/_themes/

# Set API subdirectories.
APIREF_DOXY_DIRNAME=apiref_doxy
APIREF_DOXY_DIR="./source/${APIREF_DOXY_DIRNAME}/"
APIDOC_PY_DIR=./source/apidoc_py/
APIDOC_CPP_DIR_STR=".\/source\/apidoc_cpp\/"

# Set HTML build subdirectories.
BUILD_HTML_DOXY_DIR=./build/doxygen/html/
BUILD_HTML_SPHINX_DIR=./build/sphinx/html/
BUILD_HTML_PUBLIC_DIR=./build/public_html/

# Set source directory and subdirectories.
SRC_DIR=../src/
PKG_RESOURCES_DIR="${SRC_DIR}/${PKG_NAME}/resources/"

# Set exclusion path patterns.
EXCL_APIDOC_FILES="${APIDOC_PY_DIR}/${PKG_NAME}.rst"


# -- Build Docs ----------------------------------------------------------

# @brief Get the project release number.
#
get_version_release () {
    echo $(python "${ROOT_DIR}/deploy/pkg/describe_release.py")
}

# @brief Replace strings in file.
#
# @arg File.
# @arg String to be replaced.
# @arg String replacement.
#
replace_in_file () {
    sed -i "s|${2}|${3}|g" "$1"
}

# @brief Replace key--value pairs in Doxyfile.conf.
#
# @arg Replacement string.
# @globals DOXYFILE
# @locals str_confline, str_confvar
#
replace_in_doxyfile () {
    str_confline="$1"
    str_confvar=$(printf "${str_confline}" | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s|${str_confvar} .*=.*|${str_confline}|g" "${DOXYFILE}"
}

# @brief Replace key--value pairs in Doxyfile[.conf].
#
# @arg Replacement string.
# @globals DOXYFILE_SPHINX
# @locals str_confline, str_confvar
#
recycle_doxyfile () {
    str_confline="$1"
    str_confvar=$(printf "${str_confline}" | tr -s ' ' | cut -d ' ' -f 1)
    sed -i "s|${str_confvar} .*=.*|${str_confline}|g" "${DOXYFILE_SPHINX}"
}

# @brief Modify README.md file temporarily.
#
# @arg Replacement strings.
# @globals README_FILE
# @locals str_line, str_line_new
#
recycle_readme () {
    str_line="$1"
    str_line_new="$2"
    sed -i "s|.*${str_line}.*|${str_line_new}|g" "${README_FILE}"
}

# Change to working directory.
cd "${DOCS_DIR}"

# Clean up.
make clean excl="${excl}"

# Extract version number.
vers="$(get_version_release)"
if [[ "${READTHEDOCS}" == 'True' ]]; then
    echo "${vers}" > RTD_VERSION.tmp
fi

# Backup Doxyfile.
cp "${DOXYFILE}" "${DOXYFILE}.bak"

# Pre-configure Doxygen-related docs.
replace_in_doxyfile "PROJECT_NUMBER         = ${vers}"

# HACK: Ensure RTD Doxygen backward compatibility.
if [[ "${READTHEDOCS}" == 'True' ]]; then
    MATHJAX_RELPATH="https:\/\/cdn.jsdelivr.net\/npm\/mathjax@2"
    replace_in_file ${THMS_DIR}/doxygen-header.html "\$darkmode"
    replace_in_doxyfile "MATHJAX_RELPATH        = ${MATHJAX_RELPATH}"
fi

# Prepare static assets.
cp -r "${PKG_RESOURCES_DIR}"/* "${STAT_DIR}"
cp "${ROOT_DIR}/Makefile" "${STAT_DIR}"

# Build Doxygen docs.
if [[ ${excl} != 'doxy' ]]; then
    # Temporarily alter README.md for Doxygen mainpage.
    README_FILE="${ROOT_DIR}/README.md"
    cp "${README_FILE}" "${README_FILE}.bak"
    # recycle_readme "\`\`\`bash session" "\`\`\`console"
    recycle_readme "ERC-Logo-Flag.png#gh-light-mode-only"
    # recycle_readme "\[Release\]"
    recycle_readme "\[CI\]"
    recycle_readme "\[Docs\]"
    recycle_readme "\[pre-commit.ci-Status\]"
    recycle_readme "\[Codacy-Badge\]"
    recycle_readme "\[Release-Date\]"
    recycle_readme "\[Commits-Since\]"
    recycle_readme "\[Build-Issues\]"
    recycle_readme "\[Bug-Issues\]"
    recycle_readme "\[Feature-Issues\]"
    recycle_readme "\[Pull-Requests\]"
    recycle_readme "\[pre-commit\]"

    # Build docs.
    doxygen "${DOXYFILE}"

    # Restore README.md.
    mv "${README_FILE}.bak" "${README_FILE}"
fi

# Build Sphinx docs.
if [[ ${excl} != 'sphinx' ]]; then
    # Bridge Doxygen and Sphinx docs.
    mkdir -p "${APIREF_DOXY_DIR}"
    cp -r "${BUILD_HTML_DOXY_DIR}"/** "${APIREF_DOXY_DIR}"

    # Configure Doxygen for Breathe+Exhale.
    cp "${DOXYFILE}" "${DOXYFILE_SPHINX}"

    recycle_doxyfile "OUTPUT_DIRECTORY       = ${APIDOC_CPP_DIR_STR}"
    recycle_doxyfile "GENERATE_HTML          = NO"
    recycle_doxyfile "GENERATE_XML           = YES"
    recycle_doxyfile "HAVE_DOT               = NO"
    recycle_doxyfile "HTML_HEADER            ="
    recycle_doxyfile "USE_MDFILE_AS_MAINPAGE ="
    recycle_doxyfile "EXCLUDE_PATTERNS       = triumvirate.cpp"
    sed -i "s|\.\/source|..\/source|g" "${DOXYFILE_SPHINX}"
    sed -i "s|\.\.\/src|..\/..\/src|g" "${DOXYFILE_SPHINX}"
    sed -i "s|\.\.\/README.md||g" "${DOXYFILE_SPHINX}"

    # Build docs with Breathe+Exhale.
    sphinx-apidoc -efEMT -d 1 -t "${TMPL_DIR}" -o "${APIDOC_PY_DIR}" "${SRC_DIR}"

    rm ${EXCL_APIDOC_FILES}

    if [[ "${READTHEDOCS}" != 'True' ]]; then
        make html SPHINXOPTS='-j auto'
    fi

    # Clean up.
    if [[ "${READTHEDOCS}" != 'True' ]]; then
        rm "${DOXYFILE_SPHINX}"
    fi
fi


# -- Organise Docs -------------------------------------------------------

if [[ "${READTHEDOCS}" != 'True' ]]; then
    # Set directories.
    _APIREF_DOXY="${APIREF_DOXY_DIRNAME}"
    _APIREF_DOXY_IN_HTML_SPHINX="_static/${_APIREF_DOXY}/"
    _STAT_MIRROR_IN_HTML_PUBLIC="${_APIREF_DOXY}/docs/${STAT_DIR}"

    # Move Sphinx-build HTML to public HTML directory.
    mkdir -p "${BUILD_HTML_PUBLIC_DIR}"
    if [[ -d "${BUILD_HTML_SPHINX_DIR}" ]]; then
        cp -r "${BUILD_HTML_SPHINX_DIR}"/** "${BUILD_HTML_PUBLIC_DIR}"
    fi

    # Promote Doxygen-build HTML to top-level subdirectory
    # under the public HTML directory.
    if [[ -d "${BUILD_HTML_PUBLIC_DIR}/${_APIREF_DOXY_IN_HTML_SPHINX}" ]]; then
        mv \
            "${BUILD_HTML_PUBLIC_DIR}/${_APIREF_DOXY_IN_HTML_SPHINX}" \
            "${BUILD_HTML_PUBLIC_DIR}/"
    fi

    # Move static assets to a mirrored subdirectory
    # in the public HTML directory for Doxygen 'mainpage'.
    mkdir -p "${BUILD_HTML_PUBLIC_DIR}/${_STAT_MIRROR_IN_HTML_PUBLIC}"
    find "${STAT_DIR}" -maxdepth 1 -type d -exec  \
        cp -r {} "${BUILD_HTML_PUBLIC_DIR}/${_STAT_MIRROR_IN_HTML_PUBLIC}" \;

    # Redirect HTML listings to the index page.
    echo -e "RewriteEngine On\nRewriteRule (.*)\$ index.html" \
        > "${BUILD_HTML_PUBLIC_DIR}.htaccess"
fi

# Restore backed-up files.
mv "${DOXYFILE}.bak" "${DOXYFILE}"

# Return to original directory.
cd -
