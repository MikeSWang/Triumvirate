
# Check if in an interactive terminal in Codespaces VS Code and
# first-run notice has already been shown.
if [ -t 1 ] \
    && [[ "${TERM_PROGRAM}" = "vscode" || "${TERM_PROGRAM}" = "codespaces" ]] \
    && [ -f "$HOME/.config/vscode-dev-containers/first-run-notice-already-displayed" ]; then
    # If so, show the notice again.
    if [ -f "/usr/local/etc/vscode-dev-containers/first-run-notice.txt" ]; then
        cat "/usr/local/etc/vscode-dev-containers/first-run-notice.txt"
    fi
fi
