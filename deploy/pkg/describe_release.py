"""Describe the version of the current repository as a release.

"""
import argparse

from setuptools_scm import get_version


def parse_options():
    """Parse command-line options.

    Returns
    -------
    dict
        Parsed command-line options as a dictionary.

    """
    parser = argparse.ArgumentParser(
        prog='describe-repository-release',
        argument_default=argparse.SUPPRESS
    )

    parser.add_argument(
        '--cd', type=str,
        help="Git root relative path to the script"
    )
    parser.add_argument(
        '--scheme', type=str,
        help="version scheme (see ``setuptools_scm``)"
    )
    parser.add_argument(
        '--local-scheme', type=str,
        help="local version scheme (see ``setuptools_scm``)"
    )

    return vars(parser.parse_args())


def describe_version(cd="../../",
                     version_scheme='post-release',
                     local_scheme='node-and-date'):
    """Describe the version of the current repository as a release.

    Parameters
    ----------
    cd : str
        Repository directory root relative to this file.
    version_scheme, local_scheme : str
        Version scheme and local version scheme. See
        :func:`setuptools_scm.get_version` for details.

    Returns
    -------
    str
        Version of the current repository as release.

    """
    version = get_version(
        root=cd,
        relative_to=__file__,
        version_scheme=version_scheme,
        local_scheme=local_scheme,
    )

    print(version)

    return version


if __name__ == '__main__':

    opts = parse_options()
    vers = describe_version(**opts)
