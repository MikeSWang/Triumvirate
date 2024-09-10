#!/usr/bin/env bash
#
# @file enh_codee.sh
# @author Mike S Wang
# @brief Improve code quality with Codee.
# @note Requires ``bear`` (e.g. from conda-forge) and ``codee``.
#

THIS_DIR="$(dirname -- "${BASH_SOURCE[0]}")"

make cppclean

bear --output "${THIS_DIR}/compile_commands.json" -- \
    make -j cppappbuild useomp=true

codee checks --config "${THIS_DIR}/compile_commands.json"
