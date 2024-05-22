#!/usr/bin/env bash
#
# @file buildwheel_deps_macos.sh
#
# @author Mike S Wang
# @brief Build dependencies of wheels from source natively on macOS.
#

# ---- Configuration -----------------------------------------------------

set -e -x

# Detect macOS architecture and set deployment target.
ARCH=$(uname -m)
if [[ $ARCH == "arm64" ]]; then
    export MACOSX_DEPLOYMENT_TARGET=11.0
elif [[ $ARCH == "x86_64" ]]; then
    export MACOSX_DEPLOYMENT_TARGET=10.9
fi

# Set up paths and directories.
ROOT=$(pwd)/tmp
SOURCE_DIR=${ROOT}/source
INSTALL_DIR=${ROOT}/build

mkdir -p ${SOURCE_DIR} ${INSTALL_DIR}
cd ${ROOT}


# ---- FFTW --------------------------------------------------------------

# Download and extract FFTW source.
SOURCE_FILENAME="fftw-3.3.10.tar.gz"
SOURCE_FILE=${ROOT}/${SOURCE_FILENAME}

curl -L ftp://ftp.fftw.org/pub/fftw/${SOURCE_FILENAME} -o ${SOURCE_FILE}
tar -zxvf ${SOURCE_FILE} -C ${SOURCE_DIR} && rm ${SOURCE_FILE}

# Configure and build GSL.
cd ${SOURCE_DIR}/fftw-*
mkdir -p ${INSTALL_DIR}/fftw
CC=${CXX} ./configure --enable-shared --enable-openmp --prefix=${INSTALL_DIR}/fftw
make
make install
cd -

# ---- GSL ---------------------------------------------------------------

# Download and extract GSL source.
SOURCE_FILENAME="gsl-latest.tar.gz"
SOURCE_FILE=${ROOT}/${SOURCE_FILENAME}

curl -L https://mirror.ibcp.fr/pub/gnu/gsl/${SOURCE_FILENAME} -o ${SOURCE_FILE}
tar -zxvf ${SOURCE_FILE} -C ${SOURCE_DIR} && rm ${SOURCE_FILE}

# Configure and build GSL.
cd ${SOURCE_DIR}/gsl-*
mkdir -p ${INSTALL_DIR}/gsl
CC=${CXX} ./configure --prefix=${INSTALL_DIR}/gsl
make
make install
cd -
