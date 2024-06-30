# Contribution Guidelines

User feedback and contributions are highly valued and very welcome. If
you have found an issue, please refer to '[Reporting issues](#issue)';
if you would like to contribute, please refer to
'[Contributing](#contribution)'.


## <a name="issue"></a> Reporting issues

Before reporting a new issue, please check it has not already been reported
in [Issues](https://github.com/MikeSWang/Triumvirate/issues).

If so,
[open a new issue](https://github.com/MikeSWang/Triumvirate/issues/new/choose)
with one of the templates including any relevant information. If applicable,
show a minimum (non-)working example as well as expected behaviour.

If the issue relates to build/installation, please also include the
information about your build environment.

> [!NOTE]
> If you are reporting a security issue, please refer to the
> [security policy](SECURITY.md) instead.


## <a name="contribution"></a> Contributing

To make a contribution such as a bug fix or a feature, please fork
``Triumvirate`` into your GitHub repositories and clone it to your local
machine. Always create a new branch for your edits and once complete,
commit and push it to your forked repository. Finally, [open a pull
request](https://github.com/MikeSWang/Triumvirate/pulls/compare).

If ``Triumvirate`` has been updated while you are editing your forked copy,
please merge the updated
[``main``](https://github.com/MikeSWang/Triumvirate/tree/main) branch first
before completing your edits.

Please ensure your edits conform to [PEP 8](
https://www.python.org/dev/peps/pep-0008/) standards as closely as possible,
and always include doc-strings (ideally in the [``numpydoc``](
https://numpydoc.readthedocs.io/en/latest/format.html) format) and unit tests
if you are adding new features.

A number of quality control tools are used for this project, including
GitHub Actions, Codacy and ``pre-commit``. At the minimum level, please install
``pre-commit`` in your developmental environment and run it (or install as
a Git hook) against the [``pre-commit`` configuration file](.pre-commit-config.yaml)
(which should be discovered and used by default). For more details,
please see [``pre-commit``](https://pre-commit.com) and
[``pre-commit.ci``](https://pre-commit.ci) documentation.


## Development Container

[![Codespaces](https://github.com/codespaces/badge.svg)](https://codespaces.new/MikeSWang/Triumvirate?hide_repo_select=true&ref=main)

A development container (dev container) is made available through
GitHub Codespaces. It provides a web-based cloud development environment
in the form of Visual Studio Code, and is pre-configured for building
and developing this project.
