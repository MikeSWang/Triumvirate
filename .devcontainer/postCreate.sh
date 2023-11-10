# Initialise Conda for each shell.
conda init bash
conda init zsh

# Imitate `postStartCommand` using configuration files for each shell.
# cat .devcontainer/postStart.sh >> $HOME/.bashrc
# cat .devcontainer/postStart.sh >> $HOME/.zshrc
echo "alias welcome-notice='cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt'" >> $HOME/.bashrc
echo "alias welcome-notice='cat /usr/local/etc/vscode-dev-containers/first-run-notice.txt'" >> $HOME/.zshrc

# Fetch Git tags.
git fetch --tags
