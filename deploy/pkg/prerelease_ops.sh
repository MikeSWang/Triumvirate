#!/usr/bin/env bash
#
# @file prerelease_ops.sh
# @author Mike S Wang
# @brief Pre-release operations in dry-run mode.
#

THIS_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)"
export PYTHONPATH="${PYTHONPATH}:${THIS_DIR}"

semantic-release -v changelog 2> /dev/null
semantic-release -v --noop version --no-push # 2> /dev/null
