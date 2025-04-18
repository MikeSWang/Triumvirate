# @file Makefile
# @brief `make` instructions for building and testing Triumvirate.
# @authors Mike S Wang (https://github.com/MikeSWang)
#

# ========================================================================
# Configuration
# ========================================================================

REPONAME := Triumvirate
PKGNAME := Triumvirate
MODNAME := triumvirate
PROGNAME := triumvirate
LIBNAME := trv


# ------------------------------------------------------------------------
# Preamble
# ------------------------------------------------------------------------

# Versioning
SCM_VER_SCHEME ?= no-guess-dev
SCM_LOC_SCHEME ?= node-and-date

PKG_VER := $(shell \
	python deploy/pkg/describe_release.py 2>/dev/null || \
	git describe --tag | sed -e 's/^v//'\
)

# Escape characters
COMMA := ,


# ------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------

# Repository root
DIR_ROOT := $(shell pwd)

# Package, dependency submodule, build, test and distribution directories
DIR_PKG := ${DIR_ROOT}/src/${MODNAME}
DIR_DEPDSMOD := ${DIR_ROOT}/depds
DIR_BUILD := ${DIR_ROOT}/build
DIR_TESTS := ${DIR_ROOT}/tests
DIR_DIST := ${DIR_ROOT}/dist

# Package subdirectories
DIR_PKG_INCLUDE := ${DIR_PKG}/include
DIR_PKG_SRC := ${DIR_PKG}/src
DIR_PKG_SRCPROG := ${DIR_PKG}/main

# Build subdirectories
DIR_BUILDOBJ := ${DIR_BUILD}/obj
DIR_BUILDLIB := ${DIR_BUILD}/lib
DIR_BUILDBIN := ${DIR_BUILD}/bin

# Test subdirectories
DIR_TESTBUILD := ${DIR_TESTS}/test_build
DIR_TESTOUT := ${DIR_TESTS}/test_output


# ------------------------------------------------------------------------
# Options
# ------------------------------------------------------------------------

# -- Make options --------------------------------------------------------

# Extract the '-j' or '--jobs' option (possibly empty).
# Note the space added after ``${MAKEFLAGS}``.
PATTERN_JOBS = "\-j[[:digit:][:space:]]*[^a-z[:punct:]]"
MAKEFLAGS_JOBS = $(shell echo "${MAKEFLAGS} " | grep -Eo ${PATTERN_JOBS})


# -- Custom options ------------------------------------------------------

# NOTE: Use `undefine` for make>=3.82.

# CUDA: enabled with ``usecuda=(true|1)``; disabled otherwise
ifdef usecuda
ifeq ($(strip ${usecuda}), $(filter $(strip ${usecuda}), true 1))
usecuda := true
else   # usecuda != (true|1)
unexport usecuda
endif  # usecuda == (true|1)
endif  # usecuda

# HIP: enabled with ``usehip=(true|1)``; disabled otherwise
ifdef usehip
ifeq ($(strip ${usehip}), $(filter $(strip ${usehip}), true 1))
usehip := true
else   # usehip != (true|1)
unexport usehip
endif  # usehip == (true|1)
endif  # usehip

ifdef usehip
ifdef usecuda
$(error ERROR: HIP-ported CUDA is not supported yet.)
endif  # usecuda
endif  # usehip

# OpenMP: enabled with ``useomp=(true|1)``; disabled otherwise
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))
useomp := true
else   # useomp != (true|1)
unexport useomp
endif  # useomp == (true|1)
endif  # useomp

# HDF5: enabled with ``usehdf5=(true|1)``; disabled otherwise
ifdef usehdf5
ifeq ($(strip ${usehdf5}), $(filter $(strip ${usehdf5}), true 1))
usehdf5 := true
$(warning CAUTION: HDF5 support is enabled. \
	If your HDF5 library is MPI-parallelised, check whether it is threadsafe, \
	and whether your compiler is MPI-compatible. \
	If in doubt, use a non-parallel HDF5 library.)
else   # usehdf5 != (true|1)
unexport usehdf5
endif  # usehdf5 == (true|1)
endif  # usehdf5

# Visual display: enabled with ``usedisp=(true|1)``; disabled otherwise
ifdef usedisp
ifeq ($(strip ${usedisp}), $(filter $(strip ${usedisp}), true 1))
usedisp := true
else   # usedisp != (true|1)
unexport usedisp
endif  # usedisp == (true|1)
endif  # usedisp

# Profiling: enabled with ``useprof=(true|1)``; disabled otherwise
ifdef useprof
ifeq ($(strip ${useprof}), $(filter $(strip ${useprof}), true 1))
useprof := true
else   # useprof != (true|1)
unexport useprof
endif  # useprof == (true|1)
endif  # useprof

# Parameter debugging: enabled with ``dbgpars=(true|1)``; disabled otherwise
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))
dbgpars := true
else   # dbgpars != (true|1)
unexport dbgpars
endif  # dbgpars == (true|1)
endif  # dbgpars


# ------------------------------------------------------------------------
# Compilation
# ------------------------------------------------------------------------

# -- OS ------------------------------------------------------------------

OS := $(shell uname -s)


# -- Compiler ------------------------------------------------------------

# Assume explicitly GCC compiler by default. [adapt]
ifeq (${OS}, Linux)

## If using CUDA/HIP, use CUDA/HIP compiler.
	ifdef usehip
	CXX ?= hipcc
	else   # !usehip
	ifdef usecuda
	CXX ?= nvcc
	else   # !usehip && !usecuda
	CXX ?= g++
	endif  # !usehip && usecuda
	endif  # usehip

else  # OS != Linux
ifeq (${OS}, Darwin)

	ifdef usecuda
	$(error ERROR: CUDA is not supported on macOS.)
	endif  # usecuda
	ifdef usehip
	$(error ERROR: HIP is not supported on macOS.)
	endif  # usehip

## Use GCC compiler from Homebrew (brew formula 'gcc').
## The compiler binary may have suffix '-<version>';
## check the version number with ``brew info gcc``.
	CXX ?= $(shell find $(brew --prefix gcc)/bin -type f -name 'g++*')

## Use LLVM compiler from Homebrew (brew formula 'llvm').
    # CXX ?= $(shell brew --prefix llvm)/bin/clang++

else  # OS != Linux && OS != Darwin

## If using CUDA/HIP, use CUDA/HIP compiler.
	ifdef usehip
	CXX ?= hipcc
	else   # !usehip
	ifdef usecuda
	CXX ?= nvcc
	else   # !usehip && !usecuda
	CXX ?= g++
	endif  # !usehip && usecuda
	endif  # usehip

endif  # OS != Linux && OS == Darwin
endif  # OS == Linux

# Assume default archiver. [adapt]
AR ?= ar
ARFLAGS ?= -rcsv

# Assume default remover. [adapt]
RM ?= rm -f


# -- Dependencies --------------------------------------------------------

DEPDS := gsl fftw3

# If using HDF5, add its dependencies.
ifdef usehdf5
DEPDS += hdf5
endif  # usehdf5

# Dependencies are searched for by `pkg-config`.  Ensure the set-up of
# `pkg-config` matches that of the dependencies (e.g. both are installed
# by Conda in the same Conda environment).
DEPD_INCLUDES := $(foreach dep,${DEPDS},$(shell pkg-config --silence-errors --cflags-only-I $(dep)))
DEPD_CXXFLAGS := $(foreach dep,${DEPDS},$(shell pkg-config --silence-errors --cflags-only-other $(dep)))
DEPD_LDFLAGS := $(foreach dep,${DEPDS},$(shell pkg-config --silence-errors --libs-only-other --libs-only-L $(dep)))
DEPD_LDLIBS := $(foreach dep,${DEPDS},$(shell pkg-config --silence-errors --libs-only-l $(dep)))

# If using cuFFT/hipFFT, add its dependencies.
ifdef usehip
DEPD_LDLIBS += -lhipfft
else   # !usehip
ifdef usecuda
DEPD_LDLIBS += -lcufft # -lcufftw
endif  # !usehip && usecuda
endif  # usehip


# -- Dependencies (test) -------------------------------------------------

DEPDS_TEST := gtest

DEPD_TEST_INCLUDES := $(foreach deptest,${DEPDS_TEST},$(shell pkg-config --silence-errors --cflags-only-I $(deptest)))
DEPD_TEST_CXXFLAGS := $(foreach deptest,${DEPDS_TEST},$(shell pkg-config --silence-errors --cflags-only-other $(deptest)))
DEPD_TEST_LDFLAGS := $(foreach deptest,${DEPDS_TEST},$(shell pkg-config --silence-errors --libs-only-other --libs-only-L $(deptest)))
DEPD_TEST_LDLIBS := $(foreach deptest,${DEPDS_TEST},$(shell pkg-config --silence-errors --libs-only-l $(deptest)))


# -- Options -------------------------------------------------------------

INCLUDES += -I${DIR_PKG_INCLUDE} ${DEPD_INCLUDES}

ifdef usehip
ifdef usecuda
CPPFLAGS += -D__HIP_PLATFORM_NVIDIA__
else   # usehip && !usecuda
CPPFLAGS += -D__HIP_PLATFORM_AMD__
endif  # usehip && usecuda
endif  # usehip

CPPFLAGS += -MMD -MP \
    -D__TRV_VERSION__=\"${PKG_VER}\" \
	-D__TZOFFSET__=\"$(shell date +%z)\"

ifdef usehip
ifdef usecuda
CXXFLAGS += -std=c++17 -Xcompiler -Wall,-O3 ${DEPD_CXXFLAGS}
else   # usehip && !usecuda
CXXFLAGS += -std=c++17 -Wall -Wno-vla-cxx-extension -O3 ${DEPD_CXXFLAGS}
endif  # usehip && usecuda
else   # !usehip
ifdef usecuda
CXXFLAGS += -std=c++17 -Xcompiler -Wall,-O3 ${DEPD_CXXFLAGS}
else   # !usehip && !usecuda
CXXFLAGS += -std=c++17 -Wall -O3 ${DEPD_CXXFLAGS}
endif  # !usehip && usecuda
endif  # usehip

ifeq ('$(findstring clang++,$(CXX))', '')
else
CXXFLAGS += -Wno-vla-cxx-extension
endif

ifdef usehip
LDFLAGS += \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEPD_LDFLAGS})) \
	${DEPD_LDFLAGS}
else   # !usehip
ifdef usecuda
LDFLAGS += \
	$(addprefix -Xlinker -rpath${COMMA},$(patsubst -L%,%,${DEPD_LDFLAGS})) \
	${DEPD_LDFLAGS}
else   # !usehip && !usecuda
LDFLAGS += \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEPD_LDFLAGS})) \
	${DEPD_LDFLAGS}
endif  # !usehip && usecuda
endif  # usehip

LDLIBS_CORE := -lgsl -lgslcblas -lfftw3 -lm
ifdef useomp
LDLIBS_CORE += -lfftw3_omp
endif  # useomp
ifdef usehdf5
LDLIBS_CORE += -lhdf5
endif  # usehdf5

LDLIBS += $(if ${DEPD_LDLIBS},${DEPD_LDLIBS},${LDLIBS_CORE})

PIPOPTS ?= --user


# -- Options (test) ------------------------------------------------------

INCLUDES_TEST = ${INCLUDES} ${DEPD_TEST_INCLUDES}
CXXFLAGS_TEST = ${CXXFLAGS} ${DEPD_TEST_CXXFLAGS}

ifdef usehip
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEPD_TEST_LDFLAGS})) \
	${DEPD_TEST_LDFLAGS}
else   # !usehip
ifdef usecuda
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Xlinker -rpath${COMMA},$(patsubst -L%,%,${DEPD_TEST_LDFLAGS})) \
	${DEPD_TEST_LDFLAGS}
else   # !usehip && !usecuda
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEPD_TEST_LDFLAGS})) \
	${DEPD_TEST_LDFLAGS}
endif  # !usehip && usecuda
endif  # usehip

LDLIBS_TEST = -l${LIBNAME} ${LDLIBS} \
	$(if ${DEPD_TEST_LDLIBS},${DEPD_TEST_LDLIBS},-lgtest -lpthread)


# -- Environment ---------------------------------------------------------

# NERSC computer cluster: an example environment [adapt]
ifdef NERSC_HOST

	ifdef usehip
	CXX := hipcc
	else   # !usehip
	ifdef usecuda
	CXX := nvcc
	endif  # usecuda
	endif  # usehip

## GSL library [deprecated]
	# ifdef GSL_ROOT

	# 	INCLUDES += -I${GSL_ROOT}/include

	# 	ifdef usehip
	# 	LDFLAGS += -Wl,-rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	else   # !usehip
	# 	ifdef usecuda
	# 	LDFLAGS += -Xlinker -rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	else   # !usehip && !usecuda
	# 	LDFLAGS += -Wl,-rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	endif  # !usehip && usecuda
	# 	endif  # usehip

	# endif  # GSL_ROOT

## FFTW library [deprecated]
	# ifdef FFTW_ROOT
	# ifndef usehip
	# ifndef usecuda
	# INCLUDES += -I${FFTW_INC}
	# LDFLAGS += -Wl,-rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# endif  # !usecuda
	# endif  # usehip
	# endif  # FFTW_ROOT

## cuFFT/hipFFT library
	ifdef usehip
	INCLUDES += -I${HIP_PATH}/include
	LDFLAGS += -Wl,-rpath,${HIP_PATH}/lib -L${HIP_PATH}/lib
	else   # !usehip
	ifdef usecuda
	INCLUDES += -I${NVIDIA_PATH}/math_libs/include
	LDFLAGS += -Xlinker -rpath,${NVIDIA_PATH}/math_libs/lib64 -L${NVIDIA_PATH}/math_libs/lib64
	endif  # !usehip && usecuda
	endif  # usehip

## GTEST library
	ifdef GTEST_ROOT
	INCLUDES_TEST += -I${GTEST_ROOT}/include
	ifdef usehip
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else   # !usehip
	ifdef usecuda
	LDFLAGS_TEST += -Xlinker -rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else   # !usehip && !usecuda
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	endif  # !usehip && usecuda
	endif  # usehip
	endif  # GTEST_ROOT

endif  # NERSC_HOST

# DiRAC computer cluster [adapt]
ifdef DIRAC_HOST

	ifdef usehip
	CXX := hipcc
	else   # !usehip
	ifdef usecuda
	CXX := nvcc
	endif  # usecuda
	endif  # usehip

	ifdef usehip
	HIPCLANG_PATH := /opt/rocm-6.3.2/lib/llvm
	INCLUDES += -I${HIPCLANG_PATH}/include
	LDFLAGS += -Wl,-rpath,${HIPCLANG_PATH}/lib -L${HIPCLANG_PATH}/lib
	endif  # usehip

## cuFFT/hipFFT library
	# ifdef usehip
	# else   # !usehip
	# ifdef usecuda
	# endif  # !usehip && usecuda
	# endif  # usehip

## GTEST library
	ifdef GTEST_ROOT
	INCLUDES_TEST += -I${GTEST_ROOT}/include
	ifdef usehip
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else   # !usehip
	ifdef usecuda
	LDFLAGS_TEST += -Xlinker -rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else   # !usehip && !usecuda
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	endif  # !usehip && usecuda
	endif  # usehip
	endif  # GTEST_ROOT

endif  # DIRAC_HOST


# -- Customisation -------------------------------------------------------

# OpenMP
ifdef useomp

## Assume GCC implementation by default. [adapt]
	ifeq (${OS}, Linux)

### If using CUDA, add preprocessing flags.
	    ifdef usehip
			ifdef usecuda

			CXXFLAGS_OMP ?= -Xcompiler -fopenmp
			LDFLAGS_OMP ?= -Xcompiler -fopenmp
			LDLIBS_OMP ?= -lgomp

			else   # !usecuda

			CXXFLAGS_OMP ?= -fopenmp
			LDFLAGS_OMP ?= -fopenmp
			# LDLIBS_OMP ?= -lgomp

			endif  # usecuda
		else  # !usehip
			ifdef usecuda

			CXXFLAGS_OMP ?= -Xcompiler -fopenmp
			LDFLAGS_OMP ?= -Xcompiler -fopenmp
			LDLIBS_OMP ?= -lgomp

			else   # !usecuda

#### Use GCC implementation.
			CXXFLAGS_OMP ?= -fopenmp
			LDFLAGS_OMP ?= -fopenmp
			# LDLIBS_OMP ?= -lgomp

#### Use Intel implementation.
			# CXXFLAGS_OMP ?= -qopenmp
			# LDFLAGS_OMP ?= -qopenmp
			# # LDLIBS_OMP ?= -liomp5

			endif  # usecuda
		endif  # usehip

	else  # OS != Linux
	ifeq (${OS}, Darwin)

#### Use GCC implementation.
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		# LDLIBS_OMP ?= -lgomp

#### Use LLVM implementation from Homebrew (brew formula 'libomp').
		# CXXFLAGS_OMP ?= -I$(shell brew --prefix libomp)/include -Xpreprocessor -fopenmp
		# LDFLAGS_OMP ?= -Wl,-rpath,$(shell brew --prefix libomp)/lib -L$(shell brew --prefix libomp)/lib
		# LDLIBS_OMP ?= -lomp

	else   # OS != Linux && OS != Darwin

### If using CUDA, add preprocessing flags.
#### Use GCC implementation.
		ifdef usehip
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		else   # !usehip
		ifdef usecuda  # !usehip && usecuda
		CXXFLAGS_OMP ?= -Xcompiler -fopenmp
		LDFLAGS_OMP ?= -Xcompiler -fopenmp
		else   # !usehip && !usecuda
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		endif  # !usehip && usecuda
		endif  # usehip

	endif  # OS != Linux && OS == Darwin
	endif  # OS == Linux

	ifdef usehip
	CPPFLAGS += -DTRV_USE_OMP -DTRV_USE_FFTWOMP
	else   # !usehip
	ifdef usecuda
	CPPFLAGS += -DTRV_USE_OMP -DTRV_USE_FFTWOMP
	else   # !usehip && !usecuda
	CPPFLAGS += -DTRV_USE_OMP -DTRV_USE_FFTWOMP
	endif  # !usehip && usecuda
	endif  # usehip

	CXXFLAGS += ${CXXFLAGS_OMP}
	LDFLAGS += ${LDFLAGS_OMP}

	ifdef usehip
	LDLIBS += -lfftw3_omp ${LDLIBS_OMP}
	else   # !usehip
	ifdef usecuda
	LDLIBS += -lfftw3_omp ${LDLIBS_OMP}
	else   # !usehip && !usecuda
	LDLIBS += -lfftw3_omp ${LDLIBS_OMP}
	endif  # !usehip && usecuda
	endif  # usehip

	CPPFLAGS_TEST +=
	CXXFLAGS_TEST += ${CXXFLAGS_OMP}
	LDFLAGS_TEST += ${LDFLAGS_OMP}
	LDLIBS_TEST += ${LDLIBS_OMP}

	WOMP := with

else  # !useomp

WOMP := without

endif  # useomp

# CUDA/HIP
ifdef usehip
CPPFLAGS += -DTRV_USE_HIP
else   # !usehip
ifdef usecuda
CPPFLAGS += -DTRV_USE_CUDA
endif  # !usehip && usecuda
endif  # usehip

# HDF5
ifdef usehdf5
CPPFLAGS += -DTRV_USE_H5
INCLUDES += -I${DIR_DEPDSMOD}/highfive/include
endif  # usehdf5

# Visual display
ifdef usedisp
CPPFLAGS += -DTRV_USE_DISP
endif  # usedisp

# Profiling
ifdef useprof
## Linaro MAP profiler
ifdef usehip
CXXFLAGS += -g1 -O3 -fno-inline -fno-optimize-sibling-calls
else   # !usehip
ifdef usecuda
CXXFLAGS += -g -O3 -lineinfo
else   # !usehip && !usecuda
CXXFLAGS += -g1 -O3 -fno-inline -fno-optimize-sibling-calls
endif  # !usehip && usecuda
endif  # usehip
endif  # useprof

# Parameter debugging
ifdef dbgpars
CPPFLAGS += -DDBG_MODE -DDBG_PARS
endif  # dbgpars

# Add options: e.g. {-DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_FLAG_NOAC, ...}.


# -- Parsing -------------------------------------------------------------

# Python: export build options as environmental variables.
export PY_CXX=${CXX}
export PY_INCLUDES=${INCLUDES}
export PY_CXXFLAGS=${CXXFLAGS}
export PY_LDFLAGS=${LDFLAGS}

ifndef useomp
export PY_NO_OMP
else   # useomp
export PY_CXXFLAGS_OMP=${CXXFLAGS_OMP}
export PY_LDFLAGS_OMP=${LDFLAGS_OMP} ${LDLIBS_OMP}
ifdef usehip
export PY_OMP=true
endif  # usehip
ifdef usecuda
export PY_OMP=true
endif  # usecuda
endif  # !useomp

ifdef usehip
export PY_HIP=true
endif  # usehip
ifdef usecuda
export PY_CUDA=true
endif  # usecuda

export PY_BUILD_PARALLEL=${MAKEFLAGS_JOBS}

export PY_SCM_VER_SCHEME=${SCM_VER_SCHEME}
export PY_SCM_LOC_SCHEME=${SCM_LOC_SCHEME}

# C++: strip whitespace.
CPPFLAGS := $(strip ${CPPFLAGS}) $(strip ${INCLUDES})
CXXFLAGS := $(strip ${CXXFLAGS})
LDFLAGS := $(strip ${LDFLAGS})
LDLIBS := $(strip ${LDLIBS})

CPPFLAGS_TEST := $(strip ${CPPFLAGS_TEST}) $(strip ${INCLUDES_TEST})
CXXFLAGS_TEST := $(strip ${CXXFLAGS_TEST})
LDFLAGS_TEST := $(strip ${LDFLAGS_TEST})
LDLIBS_TEST := $(strip ${LDLIBS_TEST})


# ========================================================================
# Recipes
# ========================================================================

# ------------------------------------------------------------------------
# Targets
# ------------------------------------------------------------------------

ifdef usehip
PKGSUFFIX := -HIP
PROGSUFFIX := _hip
LIBSUFFIX := _hip
ifdef usecuda
PKGSUFFIX := -HIPCUDA
PROGSUFFIX := _hipcuda
LIBSUFFIX := _hipcuda
endif  # usecuda
else   # !usehip
ifdef usecuda
PKGSUFFIX := -CUDA
PROGSUFFIX := _cuda
LIBSUFFIX := _cuda
endif  # usecuda
endif  # usehip


# ------------------------------------------------------------------------
# Building
# ------------------------------------------------------------------------

SRCS := $(wildcard ${DIR_PKG_SRC}/*.cpp)

ifdef usehip
ifdef usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_hipcuda.o)
else   # usehip && !usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_hip.o)
endif  # usehip && usecuda
else   # !usehip
ifdef usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_cuda.o)
else   # !usehip && !usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%.o)
endif  # !usehip && usecuda
endif  # usehip

DEPS := $(OBJS:.o=.d)

PROGSRC := ${DIR_PKG_SRCPROG}/${PROGNAME}.cpp
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}${PROGSUFFIX}.o

PROGEXE := ${DIR_BUILDBIN}/${PROGNAME}${PROGSUFFIX}
PROGLIB := ${DIR_BUILDLIB}/lib${LIBNAME}${LIBSUFFIX}.a


# -- Installation --------------------------------------------------------

.PHONY: install cppinstall pyinstall cpplibinstall cppappbuild \
        uninstall cppuninstall pyuninstall

install: cppinstall pyinstall

cppinstall: cppinstall_ cpplibinstall cppappbuild

cppinstall_:
	@echo "Installing ${PKGNAME}${PKGSUFFIX} C++ library/program..."

cpplibinstall: library

cppappbuild: executable

ifdef usehip
ifdef usecuda
pyinstall:
	@echo "Installing ${PKGNAME}${PKGSUFFIX} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp deploy/pkg/pyproject/.pyproject_hipcuda.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
else   # usehip && !usecuda
pyinstall:
	@echo "Installing ${PKGNAME}${PKGSUFFIX} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp deploy/pkg/pyproject/.pyproject_hip.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
endif  # usehip && usecuda
else
ifdef usecuda
pyinstall:
	@echo "Installing ${PKGNAME}${PKGSUFFIX} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp deploy/pkg/pyproject/.pyproject_cuda.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
else  # !usehip && usecuda
pyinstall:
	@echo "Installing ${PKGNAME}${PKGSUFFIX} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp deploy/pkg/pyproject/.pyproject.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
endif  # !usehip && !usecuda
endif

uninstall: cppuninstall pyuninstall

cppuninstall:
	@echo "Uninstalling Triumvirate(-CUDA/HIP/HIPCUDA) C++ library/program..."
	@echo "  removing builds..."
	@find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".git*" -exec rm -r {} +

pyuninstall:
	@echo "Uninstalling ${PKGNAME}${PKGSUFFIX} Python package (in pip mode)..."
	python -m pip uninstall -y ${PKGNAME}${PKGSUFFIX}


# -- Components ----------------------------------------------------------

.PHONY: executable library objects_

executable: ${PROGEXE}

${PROGEXE}: $(OBJS) ${PROGOBJ}
	@echo "Compiling ${PKGNAME}${PKGSUFFIX} C++ program ${WOMP} OpenMP..."
	@if [ ! -d ${DIR_BUILDBIN} ]; then \
	    echo "  making bin subdirectory in build directory..."; \
	    mkdir -p ${DIR_BUILDBIN}; \
	fi
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

library: ${PROGLIB}

${PROGLIB}: $(OBJS)
	@echo "Creating ${PKGNAME}${PKGSUFFIX} C++ library ${WOMP} OpenMP..."
	@if [ ! -d ${DIR_BUILDLIB} ]; then \
	    echo "  making lib subdirectory in build directory..."; \
	    mkdir -p ${DIR_BUILDLIB}; \
	fi
	$(AR) $(ARFLAGS) $@ $^

objects_:
	@echo "Creating ${PKGNAME}${PKGSUFFIX} C++ object files..."
	@if [ ! -d ${DIR_BUILDOBJ} ]; then \
	    echo "  making obj subdirectory in build directory..."; \
	    mkdir -p ${DIR_BUILDOBJ}; \
	fi

${PROGOBJ}: ${PROGSRC}
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

ifdef usehip
ifdef usecuda
$(OBJS): ${DIR_BUILDOBJ}/%_hipcuda.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
else   # usehip && !usecuda
$(OBJS): ${DIR_BUILDOBJ}/%_hip.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
endif  # usehip && usecuda
else   # !usehip
ifdef usecuda
$(OBJS): ${DIR_BUILDOBJ}/%_cuda.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
else   # !usehip && !usecuda
$(OBJS): ${DIR_BUILDOBJ}/%.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
endif  # !usehip && usecuda
endif  # usehip

-include $(DEPS)


# -- Configuration -------------------------------------------------------

.PHONY: checkopts

checknames:
	@echo "Checking names parsed by Makefile..."
	@echo "  REPONAME=${REPONAME}"
	@echo "  PKGNAME=${PKGNAME}${PKGSUFFIX}"
	@echo "  MODNAME=${MODNAME}"
	@echo "  PROGNAME=${PROGNAME}${PROGSUFFIX}"
	@echo "  LIBNAME=${LIBNAME}${LIBSUFFIX}"

checkopts:
	@echo "Checking options parsed by Makefile..."
	@echo "  MAKEFLAGS=${MAKEFLAGS}"
	@echo "  MAKEFLAGS_JOBS=${MAKEFLAGS_JOBS}"
	@echo "  SCM_VER_SCHEME=${SCM_VER_SCHEME}"
	@echo "  SCM_LOC_SCHEME=${SCM_LOC_SCHEME}"
	@echo "  CXX=${CXX}"
	@echo "  INCLUDES=${INCLUDES}"
	@echo "  CPPFLAGS=${CPPFLAGS}"
	@echo "  CXXFLAGS=${CXXFLAGS}"
	@echo "  LDFLAGS=${LDFLAGS}"
	@echo "  LDLIBS=${LDLIBS}"
	@echo "  CXXFLAGS_OMP=${CXXFLAGS_OMP}"
	@echo "  LDFLAGS_OMP=${LDFLAGS_OMP}"
	@echo "  LDLIBS_OMP=${LDLIBS_OMP}"
	@echo "  AR=${AR}"
	@echo "  ARFLAGS=${ARFLAGS}"
	@echo "  RM=${RM}"
	@echo "  PIPOPTS=${PIPOPTS}"


# ------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------

TEST_SRCS := $(wildcard ${DIR_TESTS}/*.cpp)
TEST_EXES := $(TEST_SRCS:${DIR_TESTS}/%.cpp=${DIR_TESTBUILD}/%)

.PHONY: test cpptest cpptest_ pytest

test: cpptest pytest

cpptest: cpptest_ library ${TEST_EXES}
	@echo "  running tests..."
	@sh -c ${TEST_EXES}

cpptest_:
	@echo "Performing ${PKGNAME} C++ tests..."
	@if [ ! -d ${DIR_TESTBUILD} ]; then \
	    echo "  making build subdirectory in test directory..."; \
	    mkdir -p ${DIR_TESTBUILD}; \
	fi
	@echo "  compiling tests..."

${TEST_EXES}: ${TEST_SRCS}
	$(CXX) $(CPPFLAGS_TEST) $(CXXFLAGS_TEST) $< -o $@ $(LDFLAGS_TEST) $(LDLIBS_TEST)

pytest:
	@echo "Performing ${PKGNAME} Python tests..."
	@if [ ! -d ${DIR_TESTOUT} ]; then \
	    echo "  making output subdirectory in test directory..."; \
	    mkdir -p ${DIR_TESTOUT}; \
	fi
	@echo "  running tests..."
	pytest -vvv


# ------------------------------------------------------------------------
# Cleaning
# ------------------------------------------------------------------------

.PHONY: clean buildclean testclean distclean runclean cppclean pyclean

clean: buildclean testclean distclean runclean

buildclean: cppclean pyclean

cppclean:
	@echo "Cleaning up Triumvirate(-CUDA/HIP/HIPCUDA) C++ build..."
	@echo "  removing builds..."
	@find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".git*" -exec rm -r {} +

pyclean:
	@echo "Cleaning up Triumvirate(-CUDA/HIP/HIPCUDA) Python build..."
	@echo "  removing Cythonised C/C++ scripts..."
	@find ${DIR_PKG} -maxdepth 1 -name "*.cpp" -exec rm {} +
	@echo "  removing Cythonised extensions..."
	@find ${DIR_PKG} -maxdepth 1 -name "*.so" -exec rm {} +
	@echo "  removing eggs..."
	@find . -type d -name "*.eggs" -exec rm -r {} +
	@find . -type d -name "*.egg-info" -exec rm -r {} +
	@echo "  removing compiled bytecode..."
	@find . -type d -name "__pycache__" -exec rm -r {} +
	@echo "  removing Jupyter notebook checkpoints..."
	@find . -type d -name ".ipynb_checkpoints" -exec rm -r {} +

testclean:
	@echo "Cleaning up Triumvirate(-CUDA/HIP/HIPCUDA) tests..."
	@echo "  removing test builds and outputs..."
	@$(RM) -r ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/*
	@echo "  removing pytest cache..."
	@find . -type d -name ".pytest_cache" -exec rm -r {} +
	@echo "  removing compiled bytecode..."
	@find . -type d -name "__pycache__" -exec rm -r {} +
	@echo "  removing core dumps..."
	@$(RM) -r core

distclean:
	@echo "Cleaning up Triumvirate(-CUDA/HIP/HIPCUDA) distributions..."
	@echo "  removing distribution outputs..."
	@$(RM) -r ${DIR_DIST}/
	@echo "  removing wheels..."
	@find . -type d -name "wheelhouse" -exec rm -r {} +
	@echo "  removing eggs..."
	@find . -name "*.eggs" -exec rm -r {} +
	@find . -name "*.egg-info" -exec rm -r {} +

runclean:
	@echo "Cleaning up Triumvirate(-CUDA/HIP/HIPCUDA) runs..."
	@echo "  removing compiled bytecode..."
	@find . -type d -name "__pycache__" -exec rm -r {} +
	@echo "  removing core dumps..."
	@$(RM) -r core
