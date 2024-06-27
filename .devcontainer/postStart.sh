#!/usr/bin/env sh
#
# @file postStart.sh
# @author Mike S Wang
#

# Check if in an interactive terminal in Codespaces VS Code and
# first-run notice has already been shown.
if [ -t 1 ] \
    && [[ "${TERM_PROGRAM}" = "vscode" || "${TERM_PROGRAM}" = "codespaces" ]] \
    && [ -f "$HOME/.config/vscode-dev-containers/first-run-notice-already-displayed" ]; then
    # If so, show the notice again.
    NOTICE_FILE="/usr/local/etc/vscode-dev-containers/first-run-notice.txt"
    if [ -f "${NOTICE_FILE}" ]; then
        cat "${NOTICE_FILE}"
    fi
fi
