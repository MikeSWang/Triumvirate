#!/usr/bin/env bash
#
# @file autoinc_version.sh
# @author Mike S Wang
# @brief Increment fallback version number based on git tags.
#

# @brief Get version number of the latest release from git tags.
#
# @locals release_local, release_latest
#
get_latest_release () {
    # release_tags=$(git tag -l --sort=-v:refname | grep -E "^v[[:digit:]]+")
    release_local="$(git describe --tag)"
    # release_local=$(python deploy/pkg/describe_release.py)

    # release_latest=$(echo ${release_tags} | head -n 1 | cut -d ' ' -f 1)
    release_latest=$(echo "${release_local}" | cut -d '-' -f 1)

    echo "${release_latest#v}"
}

# @brief Set fallback version number in C++ header file.
#
# @arg Version number.
# @locals vers, versfile, versline, versfill
#
set_fallback_version_trv () {
    vers="$1"
    versfile=src/triumvirate/include/monitor.hpp
    versline="#define __TRV_VERSION__ .*"
    versfill="#define __TRV_VERSION__ \"${vers}\"  \/\/ (fallback) version number"
    sed -i "s/${versline}/${versfill}/g" "${versfile}"
}

# @brief Set fallback version number in ``__init__.py`` file.
#
# @arg Version number.
# @locals vers, versfile, versline, versfill
#
set_fallback_version_initpy () {
    vers="$1"
    versfile=src/triumvirate/__init__.py
    versline="__version__ = '.*'.*"
    versfill="__version__ = '${vers}'  # fallback version number"
    sed -i "s/${versline}/${versfill}/g" "${versfile}"
}

# @brief Set fallback version number in ``meta.yaml`` files.
#
# @arg Version number.
# @locals vers, versline, versfill
#
set_fallback_version_metayaml () {
    vers="$1"

    versline="{% set version = environ.get('GIT_DESCRIBE_TAG', .*) %}"
    versfill="{% set version = environ.get('GIT_DESCRIBE_TAG', 'v${vers}') %}"
    for versfile in $(find deploy/pkg -type f -name 'meta.yaml'); do
        sed -i "s/${versline}/${versfill}/g" "${versfile}"
    done

    versline="# git_rev: .*"
    versfill="# git_rev: v${vers}"
    for versfile in $(find deploy/pkg -type f -name 'meta.yaml'); do
        sed -i "s/${versline}/${versfill}/g" "${versfile}"
    done
}

# @brief Set fallback version number in ``versions.json`` files.
#
# @arg Version number.
# @locals vers, version_switcher_file, versline, versfill
#
set_fallback_version_switcher () {
    vers="$1"
    version_switcher_file=docs/versions.json

    versline="\"name\": \".* (stable)\","
    versfill="\"name\": \"${vers} (stable)\","
    sed -i "s/${versline}/${versfill}/g" "${version_switcher_file}"

    versline="\"version\": \"[[:digit:]].*\","
    versfill="\"version\": \"${vers}\","
    sed -i "s/${versline}/${versfill}/g" "${version_switcher_file}"
}

# Increment fallback version based on the latest-release git tag.
set_fallback_version_trv "$(get_latest_release)"
set_fallback_version_initpy "$(get_latest_release)"
set_fallback_version_metayaml "$(get_latest_release)"
set_fallback_version_switcher "$(get_latest_release)"
