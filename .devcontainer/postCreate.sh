#!/usr/bin/env sh
#
# @file postCreate.sh
# @author Mike S Wang
#

# Initialise Conda for each shell.
conda init bash
conda init zsh

# Imitate `postStartCommand` using configuration files for each shell.
# cat .devcontainer/postStart.sh >> "$HOME"/.bashrc
# cat .devcontainer/postStart.sh >> "$HOME"/.zshrc
NOTICE_FILE="/usr/local/etc/vscode-dev-containers/first-run-notice.txt"
echo "alias welcome-notice='cat ${NOTICE_FILE}'" >> "$HOME"/.bashrc
echo "alias welcome-notice='cat ${NOTICE_FILE}'" >> "$HOME"/.zshrc

# Fetch Git tags.
git fetch --tags
