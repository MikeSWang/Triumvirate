#!/usr/bin/env bash
#
# @file autoinc_vers.sh
# @author Mike S Wang
# @brief Increment fallback version number based on git tags.
#

# @brief Get version number of the latest release from git tags.
#
get_latest_release () {
    releases=$(git tag -l --sort=-v:refname | grep -E "^v[[:digit:]]+")
    release_latest=$(printf ${releases} | cut -d ' ' -f 1)
    printf $(printf ${release_latest} | sed -e "s/^v//")
}

# @brief Set fallback version number in relevant files.
#
# @arg Version number.
#
set_fallback_version () {
    vers=$1
    versfile=src/triumvirate/__init__.py
    versline="__version__ = '.*'.*"
    versfill="__version__ = '${vers}'  # fallback version number"
    sed -i "s/${versline}/${versfill}/g" ${versfile}
}

# Increment fallback version based on the latest-release git tag.
set_fallback_version $(get_latest_release)
