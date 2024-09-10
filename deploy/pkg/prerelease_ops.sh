#!/usr/bin/env bash
#
# @file prerelease_ops.sh
# @author Mike S Wang
# @brief Pre-release operations in dry-run mode.
#

THIS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
CONFIG_FILE="${THIS_DIR}/../../.pyproject.toml"
export PYTHONPATH="${PYTHONPATH}:${THIS_DIR}"

semantic-release -v --config ${CONFIG_FILE} changelog 2> /dev/null
semantic-release -v --config ${CONFIG_FILE} --noop version # 2> /dev/null
