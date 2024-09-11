#!/usr/bin/env bash
#
# @file enh_codee.sh
# @author Mike S Wang
# @brief Improve code quality with Codee.
# @note Requires ``bear`` (e.g. from conda-forge) and ``codee``.
#

THIS_DIR="$(dirname -- "${BASH_SOURCE[0]}")"
AUX_DIR="${THIS_DIR}/auxiliary"

mkdir -p "${AUX_DIR}"

make cppclean

bear --output "${AUX_DIR}/compile_commands.json" -- \
    make -j cppappbuild useomp=true

codee checks --config "${AUX_DIR}/compile_commands.json" 2>&1 | tee "${AUX_DIR}/codee.log"
