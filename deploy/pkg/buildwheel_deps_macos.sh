#!/usr/bin/env bash
#
# @file buildwheel_deps_macos.sh
# @author Mike S Wang
# @brief Build dependencies of wheels from source natively on macOS.
#

# ---- Configuration -----------------------------------------------------

set -x

# Detect macOS architecture and set deployment target.
ARCH=$(uname -m)
if [[ ${ARCH} == "arm64" ]]; then
    export MACOSX_DEPLOYMENT_TARGET=11.0
elif [[ ${ARCH} == "x86_64" ]]; then
    export MACOSX_DEPLOYMENT_TARGET=10.13
fi

# Set up paths and directories.
ROOT="$(pwd)/tmp"
SOURCE_DIR="${ROOT}/source"
BUILD_DIR="${ROOT}/build"
INSTALL_DIR="${ROOT}/deps"

mkdir -p "${SOURCE_DIR}" "${BUILD_DIR}" "${INSTALL_DIR}/include" "${INSTALL_DIR}/lib"
cd "${ROOT}"


# ---- OMP ---------------------------------------------------------------

# Download and extract OpenMP libraries.
if [[ ${ARCH} == "arm64" ]]; then
    SOURCE_ARCH="osx-arm64"
    SOURCE_FILENAME="llvm-openmp-18.1.5-hde57baf_0.conda"
elif [[ ${ARCH} == "x86_64" ]]; then
    SOURCE_ARCH="osx-64"
    SOURCE_FILENAME="llvm-openmp-18.1.5-h39e0ece_0.conda"
fi

SOURCE_FILE="${SOURCE_DIR}/${SOURCE_FILENAME}"

mkdir -p "${SOURCE_DIR}/llvm-openmp" "${BUILD_DIR}/llvm-openmp"
curl -L https://anaconda.org/conda-forge/llvm-openmp/18.1.5/download/${SOURCE_ARCH}/${SOURCE_FILENAME} -o "${SOURCE_FILE}"
tar -zxvf "${SOURCE_FILE}" -C "${SOURCE_DIR}/llvm-openmp" && rm "${SOURCE_FILE}"
tar -zxvf "${SOURCE_DIR}"/llvm-openmp/pkg-llvm-openmp-* -C "${BUILD_DIR}/llvm-openmp"

# Copy OpenMP libraries to installation directory.
cp -r "${BUILD_DIR}/llvm-openmp/include/." "${INSTALL_DIR}/include"
cp -r "${BUILD_DIR}/llvm-openmp/lib/." "${INSTALL_DIR}/lib"


# ---- FFTW --------------------------------------------------------------

# ------------------------------------------------------------------------
# # Download and extract FFTW source.
# SOURCE_FILENAME="fftw-3.3.10.tar.gz"
# SOURCE_FILE="${ROOT}/${SOURCE_FILENAME}"

# curl -L ftp://ftp.fftw.org/pub/fftw/${SOURCE_FILENAME} -o "${SOURCE_FILE}"
# tar -zxvf "${SOURCE_FILE}" -C "${SOURCE_DIR}" && rm "${SOURCE_FILE}"

# # Configure and build GSL.
# cd "${SOURCE_DIR}"/fftw-*
# mkdir -p "${BUILD_DIR}/fftw"
# CC=${CXX} ./configure --enable-shared --enable-openmp --prefix="${BUILD_DIR}/fftw"
# make
# make install
# cd -
# ------------------------------------------------------------------------

# Download and extract FFTW libraries.
if [[ ${ARCH} == "arm64" ]]; then
    SOURCE_ARCH="osx-arm64"
    SOURCE_FILENAME="fftw-3.3.10-nompi_h3046061_108.conda"
elif [[ ${ARCH} == "x86_64" ]]; then
    SOURCE_ARCH="osx-64"
    SOURCE_FILENAME="fftw-3.3.10-nompi_h4fa670e_108.conda"
fi

SOURCE_FILE="${ROOT}/${SOURCE_FILENAME}"

mkdir -p "${SOURCE_DIR}/fftw-3.3.10" "${BUILD_DIR}/fftw-3.3.10"
curl -L https://anaconda.org/conda-forge/fftw/3.3.10/download/${SOURCE_ARCH}/${SOURCE_FILENAME} -o "${SOURCE_FILE}"
tar -zxvf "${SOURCE_FILE}" -C "${SOURCE_DIR}/fftw-3.3.10" && rm "${SOURCE_FILE}"
tar -zxvf "${SOURCE_DIR}"/fftw-3.3.10/pkg-fftw-3.3.10-* -C "${BUILD_DIR}"/fftw-3.3.10

# Move FFTW libraries to installation directory.
cp -r "${BUILD_DIR}/fftw-3.3.10/include/." "${INSTALL_DIR}/include"
cp -r "${BUILD_DIR}/fftw-3.3.10/lib/." "${INSTALL_DIR}/lib"


# ---- GSL ---------------------------------------------------------------

# ------------------------------------------------------------------------
# # Download and extract GSL source.
# SOURCE_FILENAME="gsl-latest.tar.gz"
# SOURCE_FILE="${ROOT}/${SOURCE_FILENAME}"

# curl -L https://mirror.ibcp.fr/pub/gnu/gsl/${SOURCE_FILENAME} -o "${SOURCE_FILE}"
# tar -zxvf "${SOURCE_FILE}" -C "${SOURCE_DIR}" && rm "${SOURCE_FILE}"

# # Configure and build GSL.
# cd "${SOURCE_DIR}"/gsl-*
# mkdir -p "${BUILD_DIR}/gsl"
# CC=${CXX} ./configure --prefix="${BUILD_DIR}/gsl"
# make
# make install
# cd -
# ------------------------------------------------------------------------

# Download and extract GSL libraries.
if [[ ${ARCH} == "arm64" ]]; then
    SOURCE_ARCH="osx-arm64"
    SOURCE_FILENAME="gsl-2.7-h6e638da_0.tar.bz2"
elif [[ ${ARCH} == "x86_64" ]]; then
    SOURCE_ARCH="osx-64"
    SOURCE_FILENAME="gsl-2.7-h93259b0_0.tar.bz2"
fi

SOURCE_FILE="${SOURCE_DIR}/${SOURCE_FILENAME}"

mkdir -p "${BUILD_DIR}/gsl-2.7"
curl -L https://anaconda.org/conda-forge/gsl/2.7/download/${SOURCE_ARCH}/${SOURCE_FILENAME} -o "${SOURCE_FILE}"
tar -zxvf "${SOURCE_FILE}" -C "${BUILD_DIR}/gsl-2.7"

# Copy GSL libraries to installation directory.
cp -r "${BUILD_DIR}/gsl-2.7/include/." "${INSTALL_DIR}/include"
find "${BUILD_DIR}"/gsl-2.7/lib/* -type d -exec cp -r {} "${INSTALL_DIR}/lib" \;
find "${BUILD_DIR}"/gsl-2.7/lib/* -type f \( -name "*.a" -or -name "libgsl.*.dylib" \) -exec cp -r {} "${INSTALL_DIR}/lib" \;
