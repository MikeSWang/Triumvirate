#!/usr/bin/env bash
#
# @file autoinc_version.sh
# @author Mike S Wang
# @brief Increment fallback version number based on git tags.
#

# @brief Get version number of the latest release from git tags.
#
get_latest_release () {
    # release_tags=$(git tag -l --sort=-v:refname | grep -E "^v[[:digit:]]+")
    release_local=$(git describe --tag)
    # release_local=$(python deploy/pkg/describe_release.py)

    # release_latest=$(printf ${release_tags} | cut -d ' ' -f 1)
    release_latest=$(printf ${release_local} | cut -d '-' -f 1)

    printf $(printf ${release_latest} | sed -e "s/^v//")
}

# @brief Set fallback version number in ``__init__.py`` file.
#
# @arg Version number.
#
set_fallback_version_initpy () {
    vers=$1
    versfile=src/triumvirate/__init__.py
    versline="__version__ = '.*'.*"
    versfill="__version__ = '${vers}'  # fallback version number"
    sed -i "s/${versline}/${versfill}/g" ${versfile}
}

# @brief Set fallback version number in ``meta.yaml`` files.
#
# @arg Version number.
#
set_fallback_version_metayaml () {
    vers=$1
    versline="{% set version = environ.get('GIT_DESCRIBE_TAG', .*) %}"
    versfill="{% set version = environ.get('GIT_DESCRIBE_TAG', 'v${vers}') %}"
    for versfile in $(find deploy/pkg -type f -name "meta.yaml"); do
        sed -i "s/${versline}/${versfill}/g" ${versfile}
    done
}

# Increment fallback version based on the latest-release git tag.
set_fallback_version_initpy $(get_latest_release)
set_fallback_version_metayaml $(get_latest_release)
