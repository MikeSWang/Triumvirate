#!/usr/bin/env sh
#
# @file onCreate.sh
# @author Mike S Wang
#

# Replace the notice message upon container creation.
sudo cp \
    .devcontainer/welcome_notice.txt \
    /usr/local/etc/vscode-dev-containers/first-run-notice.txt
