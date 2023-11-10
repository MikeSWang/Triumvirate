# Initialise Conda for each shell.
conda init bash
conda init zsh

# # Imitate `postStartCommand` using configuration files for each shell.
# cat .devcontainer/postStart.sh >> $HOME/.bashrc
# cat .devcontainer/postStart.sh >> $HOME/.zshrc

# Fetch Git tags.
git fetch --tags
