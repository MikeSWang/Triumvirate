"""List package artifacts for a release tag.

"""
import json
import os
import re
import subprocess as subp
from argparse import ArgumentParser
from functools import reduce
from operator import getitem
from pprint import pprint


def configure():
    """Configure command-line arguments.

    Returns
    -------
    args : argparse.Namespace
        Parsed command-line arguments.

    """
    parser = ArgumentParser(
        description='List package artifact IDs for a release tag.'
    )

    parser.add_argument('tag', help='release tag')
    parser.add_argument(
        '--github-owner',
        default=os.getenv('GITHUB_OWNER', 'MikeSWang'),
        help="GitHub owner"
    )
    parser.add_argument(
        '--github-repo',
        default=os.getenv('GITHUB_REPO', 'Triumvirate'),
        help="GitHub repository"
    )
    parser.add_argument(
        '--filter', action='append', nargs=2,
        metavar=('KEYCHAIN', 'VALUE_PATTERN'),
        default=[],
        help="additive filter keychain and value pattern, e.g. "
             "``--filter 'name' 'pypi_dist_*' "
             "--filter 'workflow_run.head_branch' 'main'`` "
             "filters artifacts with name starting with 'pypi_dist_' "
             "and workflow run head branch 'main'"
    )
    parser.add_argument(
        '--id-only', action='store_true',
        help="list only artifact IDs"
    )

    args = parser.parse_args()

    return args


def get_artifact_data(config):
    """Get package artifact data for a release tag.

    Parameters
    ----------
    config : argparse.Namespace
        Configured parameters.

    Returns
    -------
    artifact_data : dict
        Package artifact data.

    """
    github_header_accept = os.getenv(
        'GITHUB_HEADER_ACCEPT',
        "Accept: application/vnd.github+json"
    )

    github_header_auth = os.getenv('GITHUB_HEADER_AUTH')
    if github_header_auth is None:
        github_header_auth_prefix = "Authorization: token "
        github_token_artifacts = \
            os.getenv('GITHUB_TOKEN_ARTIFACTS') or os.getenv('GITHUB_TOKEN')
        if github_token_artifacts is None:
            raise RuntimeError(
                "Neither of the following environment variable is set: "
                "GITHUB_TOKEN, GITHUB_TOKEN_ARTIFACTS"
            )
        github_header_auth = github_header_auth_prefix + github_token_artifacts

    github_url_actions = os.getenv(
        'GITHUB_URL_ACTIONS',
        "https://api.github.com/repos/"
        f"{config.github_owner}/{config.github_repo}/"
        "actions/artifacts"
    )

    cmd_fields = [
        'curl', '-s', '-L',
        '-H', github_header_accept,
        '-H', github_header_auth,
        github_url_actions
    ]

    std_result = subp.run(cmd_fields, stdout=subp.PIPE)
    artifact_data = json.loads(std_result.stdout.decode('utf-8'))

    return artifact_data


def filter_pkg_artifacts(artifact_data, config):
    """Filter package artifacts.

    Parameters
    ----------
    artifact_data : dict
        Package artifact data.
    config : argparse.Namespace
        Configured parameters.

    Returns
    -------
    artifact_list : list of dict
        Filtered package artifact data.

    """
    filters = config.filter or []

    artifact_list = []
    for artifact in artifact_data['artifacts']:
        for (keychain, pattern) in filters:
            excluded = not bool(
                re.match(
                    pattern,
                    reduce(getitem, keychain.split('.'), artifact),
                )
            )
            if excluded:
                break
        else:
            artifact_list.append(artifact)

    return artifact_list


def get_latest_artifacts(artifacts_list):
    """Get only the latest artifacts when multiple IDs have the same name.

    Parameters
    ----------
    artifacts_list : list of dict
        Package artifact list.

    Returns
    -------
    artifacts_latest_list : list of dict
        Latest package artifact list.

    """
    unique_names = {artifact['name'] for artifact in artifacts_list}

    latest_ids = {}
    for name in unique_names:
        latest_ids[name] = max(
            artifact['id']
            for artifact in artifacts_list
            if artifact['name'] == name
        )

    artifacts_latest_list = [
        artifact for artifact in artifacts_list
        if artifact['id'] in latest_ids.values()
    ]

    return artifacts_latest_list


def list_pkg_artifacts(artifacts_list, id_only=False):
    """List package artifacts.

    Parameters
    ----------
    artifacts_list : list of dict
        Package artifact list.
    id_only : bool, optional
        If True, only print artifact IDs.

    """
    if id_only:
        for artifact in artifacts_list:
            print(artifact['id'])
    else:
        for artifact in artifacts_list:
            pprint(artifact)


if __name__ == '__main__':

    cfg = configure()
    artifact_data = get_artifact_data(cfg)
    artifact_list = get_latest_artifacts(
        filter_pkg_artifacts(artifact_data, cfg)
    )
    list_pkg_artifacts(artifact_list, id_only=cfg.id_only)
