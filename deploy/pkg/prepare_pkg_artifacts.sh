#!/usr/bin/env sh
#
# @file prepare_pkg_artifacts.sh
# @author Mike S Wang
# @brief Prepare package artifacts from GitHub Actions.
# @arg release_tag, the release tag to download artifacts for.
#

release_tag=$1

tmp_dir=$(mktemp -d)

# GitHub attributes
GITHUB_OWNER=MikeSWang
GITHUB_REPO=Triumvirate
GITHUB_TOKEN_ARTIFACTS="$(sed '2q;d' ~/.github_auth)"
GITHUB_URL_ACTIONS="https://api.github.com/repos/${GITHUB_OWNER}/${GITHUB_REPO}/actions/artifacts"

GITHUB_HEADER_ACCEPT="Accept: application/vnd.github+json"
GITHUB_HEADER_AUTH="Authorization: Bearer ${GITHUB_TOKEN_ARTIFACTS}"

# Query artifact IDs.
# e.g.
# artifact_names="pypi_dist_macos-13_v0.4.7.post1 conda_dist_macos-13_v0.4.7.post1"
artifact_ids=()
artifact_names="pypi_dist.*${release_tag} conda_dist.*${release_tag}"
for artifact_name in ${artifact_names};
do
    id=$(python deploy/pkg/list_pkg_artifacts.py ${release_tag} --filter name $artifact_name --id-only)
    artifact_ids+=(${id})
done

# Download and unzip artifacts.
dist_dir=${tmp_dir}/dist
mkdir -p ${dist_dir}

echo "Downloading artifacts to temporary directory for extraction: ${tmp_dir}"
for artifact_id in "${artifact_ids[@]}";
do
    archive_file=${dist_dir}/${artifact_id}.zip
    curl -L -H "${GITHUB_HEADER_ACCEPT}" -H "${GITHUB_HEADER_AUTH}" \
        -o ${archive_file} "${GITHUB_URL_ACTIONS}/${artifact_id}/zip"
    unzip ${archive_file} -d ${dist_dir} && rm ${archive_file}
done

# Organise artifacts.
release_dir=${tmp_dir}/releases

mkdir -p ${release_dir}/github-release
find ${dist_dir} -type f -name '*.conda' -exec mv {} ${release_dir}/github-release \;

mkdir -p ${release_dir}/pypi-release
find ${dist_dir} -type f \( -name '*.tar.gz' -or -name '*.whl' \) -exec mv {} ${release_dir}/pypi-release \;

mkdir -p ${release_dir}/conda-release
find ${dist_dir} -type f -name '*.tar.bz2' -exec mv {} ${release_dir}/conda-release \;

find ${release_dir} -mindepth 1 -type d ! -name '*-release' -exec rm -rf {} +

mv ${release_dir} ~/Downloads/

rm -rf ${tmp_dir}
