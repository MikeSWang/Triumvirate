#!/usr/bin/env bash
#
# @file prerelease_ops.sh
# @author Mike S Wang
# @brief Pre-release operations.
#

semantic-release --noop version --print 2> /dev/null
semantic-release changelog 2> /dev/null
