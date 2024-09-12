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

# Package, build, test and distribution directories
DIR_PKG := ${DIR_ROOT}/src/${MODNAME}
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

# CUDA: enabled with ``useomp=(true|1)``; disabled otherwise
ifdef usecuda
ifeq ($(strip ${usecuda}), $(filter $(strip ${usecuda}), true 1))
usecuda := true
	PKGNAME := ${PKGNAME}-CUDA
	PROGNAME := ${PROGNAME}_cuda
	LIBNAME := ${LIBNAME}_cuda
else   # usecuda != (true|1)
unexport usecuda
endif  # usecuda == (true|1)
endif  # usecuda

# HIP: enabled with ``usehip=(true|1)``; disabled otherwise
ifdef usehip
ifeq ($(strip ${usehip}), $(filter $(strip ${usehip}), true 1))
usehip := true
	ifndef usecuda
	PKGNAME := ${PKGNAME}-HIP
	PROGNAME := ${PROGNAME}_hip
	LIBNAME := ${LIBNAME}_hip
	else   # usecuda
	PKGNAME := ${PKGNAME}-HIPCUDA
	PROGNAME := ${PROGNAME}_hipcuda
	LIBNAME := ${LIBNAME}_hipcuda
	endif  # !usecuda
else   # usehip != (true|1)
unexport usehip
endif  # usehip == (true|1)
endif  # usehip

# OpenMP: enabled with ``useomp=(true|1)``; disabled otherwise
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))
useomp := true
else   # useomp != (true|1)
unexport useomp
endif  # useomp == (true|1)
endif  # useomp

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
	else ifdef usecuda  # !usehip && usecuda
	CXX ?= nvcc
	else                # !usehip && !usecuda
	CXX ?= g++
	endif               # usehip

else ifeq (${OS}, Darwin)

	ifdef usecuda
	$(error "CUDA is not supported on macOS.")
	endif  # usecuda
	ifdef usehip
	$(error "HIP is not supported on macOS.")
	endif  # usehip

## Use GCC compiler from Homebrew (brew formula 'gcc').
## The compiler binary may have suffix '-<version>';
## check the version number with ``brew info gcc``.
	CXX ?= $(shell find $(brew --prefix gcc)/bin -type f -name 'g++*')

## Use LLVM compiler from Homebrew (brew formula 'llvm').
    # CXX ?= $(shell brew --prefix llvm)/bin/clang++

else  # OS

## If using CUDA/HIP, use CUDA/HIP compiler.
	ifdef usehip
	CXX ?= hipcc
	else ifdef usecuda  # !usehip && usecuda
	CXX ?= nvcc
	else                # !usehip && !usecuda
	CXX ?= g++
	endif               # usehip

endif  # OS

# Assume default archiver. [adapt]
AR ?= ar
ARFLAGS ?= -rcsv

# Assume default remover. [adapt]
RM ?= rm -f


# -- Dependencies --------------------------------------------------------

# If using cuFFT/hipFFT, remove standard FFTW dependency.
ifdef usehip
DEPS := gsl
else ifdef usecuda  # !usehip && usecuda
DEPS := gsl
else                # !usehip && !usecuda
DEPS := gsl fftw3
endif               # usehip

# Dependencies are searched for by `pkg-config`.  Ensure the set-up of
# `pkg-config` matches that of the dependencies (e.g. both are installed
# by Conda in the same Conda environment).
DEP_INCLUDES := $(shell pkg-config --silence-errors --cflags-only-I ${DEPS})
DEP_CXXFLAGS := $(shell pkg-config --silence-errors --cflags-only-other ${DEPS})
DEP_LDFLAGS := $(shell pkg-config --silence-errors --libs-only-other --libs-only-L ${DEPS})
DEP_LDLIBS := $(shell pkg-config --silence-errors --libs-only-l ${DEPS})

# If using cuFFT/hipFFT, add its dependencies.
ifdef usehip
DEP_LDLIBS += -lhipfft
else ifdef usecuda  # !usehip && usecuda
DEP_LDLIBS += -lcufft -lcufftw
endif               # usehip


# -- Dependencies (test) -------------------------------------------------

DEPS_TEST := gtest

DEP_TEST_INCLUDES := $(shell pkg-config --silence-errors --cflags-only-I ${DEPS_TEST})
DEP_TEST_CXXFLAGS := $(shell pkg-config --silence-errors --cflags-only-other ${DEPS_TEST})
DEP_TEST_LDFLAGS := $(shell pkg-config --silence-errors --libs-only-other --libs-only-L ${DEPS_TEST})
DEP_TEST_LDLIBS := $(shell pkg-config --silence-errors --libs-only-l ${DEPS_TEST})


# -- Options -------------------------------------------------------------

INCLUDES += -I${DIR_PKG_INCLUDE} ${DEP_INCLUDES}
CPPFLAGS += -MMD -MP -D__TRV_VERSION__=\"${PKG_VER}\"

ifdef usehip
CXXFLAGS += -std=c++17 -Wall -O3 ${DEP_CXXFLAGS}
else ifdef usecuda  # !usehip && usecuda
CXXFLAGS += -std=c++17 -Xcompiler -Wall,-O3 ${DEP_CXXFLAGS}
else                # !usehip && !usecuda
CXXFLAGS += -std=c++17 -Wall -O3 ${DEP_CXXFLAGS}
endif               # usehip

ifdef usehip
LDFLAGS += \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEP_LDFLAGS})) \
	${DEP_LDFLAGS}
else ifdef usecuda  # !usehip && usecuda
LDFLAGS += \
	$(addprefix -Xlinker -rpath${COMMA},$(patsubst -L%,%,${DEP_LDFLAGS})) \
	${DEP_LDFLAGS}
else                # !usehip && !usecuda
LDFLAGS += \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEP_LDFLAGS})) \
	${DEP_LDFLAGS}
endif               # usehip

LDLIBS += $(if ${DEP_LDLIBS},${DEP_LDLIBS},-lgsl -lgslcblas -lfftw3 -lm)

PIPOPTS ?= --user


# -- Options (test) ------------------------------------------------------

INCLUDES_TEST = ${INCLUDES} ${DEP_TEST_INCLUDES}
CXXFLAGS_TEST = ${CXXFLAGS} ${DEP_TEST_CXXFLAGS}

ifdef usehip
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEP_TEST_LDFLAGS})) \
	${DEP_TEST_LDFLAGS}
else ifdef usecuda  # !usehip && usecuda
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Xlinker -rpath${COMMA},$(patsubst -L%,%,${DEP_TEST_LDFLAGS})) \
	${DEP_TEST_LDFLAGS}
else                # !usehip && !usecuda
LDFLAGS_TEST = -L${DIR_BUILDLIB} ${LDFLAGS} \
	$(addprefix -Wl${COMMA}-rpath${COMMA},$(patsubst -L%,%,${DEP_TEST_LDFLAGS})) \
	${DEP_TEST_LDFLAGS}
endif               # usehip

LDLIBS_TEST = -l${LIBNAME} ${LDLIBS} \
	$(if ${DEP_TEST_LDLIBS},${DEP_TEST_LDLIBS},-lgtest -lpthread)


# -- Environment ---------------------------------------------------------

# NERSC computer cluster: an example environment [adapt]
ifdef NERSC_HOST

## GSL library [deprecated]
	# ifdef GSL_ROOT

	# 	INCLUDES += -I${GSL_ROOT}/include

	# 	ifdef usehip
	# 	LDFLAGS += -Wl,-rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	else ifdef usecuda  # !usehip && usecuda
	# 	LDFLAGS += -Xlinker -rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	else                # !usehip && !usecuda
	# 	LDFLAGS += -Wl,-rpath,${GSL_ROOT}/lib -L${GSL_ROOT}/lib
	# 	endif               # usehip

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
	else ifdef usecuda  # !usehip && usecuda
	INCLUDES += -I${NVIDIA_PATH}/math_libs/include
	LDFLAGS += -Xlinker -rpath,${NVIDIA_PATH}/math_libs/lib64 -L${NVIDIA_PATH}/math_libs/lib64
	endif               # usehip

## GTEST library
	ifdef GTEST_ROOT
	INCLUDES_TEST += -I${GTEST_ROOT}/include
	ifdef usehip
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else ifdef usecuda  # !usehip && usecuda
	LDFLAGS_TEST += -Xlinker -rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else                # !usehip && !usecuda
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	endif               # usehip
	endif  # GTEST_ROOT

endif  # NERSC_HOST

# DiRAC computer cluster [adapt]
ifdef DIRAC_HOST

## GTEST library
	ifdef GTEST_ROOT
	INCLUDES_TEST += -I${GTEST_ROOT}/include
	ifdef usehip
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else ifdef usecuda  # !usehip && usecuda
	LDFLAGS_TEST += -Xlinker -rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	else                # !usehip && !usecuda
	LDFLAGS_TEST += -Wl,-rpath,${GTEST_ROOT}/lib -L${GTEST_ROOT}/lib
	endif               # usehip
	endif  # GTEST_ROOT

endif  # DIRAC_HOST


# -- Customisation -------------------------------------------------------

# OpenMP
ifdef useomp

## Assume GCC implementation by default. [adapt]
	ifeq (${OS}, Linux)

### If using CUDA, add preprocessing flags.
	    ifdef usehip

#### Use GCC implementation.
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		# LDLIBS_OMP ?= -lgomp

		else ifdef usecuda  # !usehip && usecuda

#### Use GCC implementation.
		CXXFLAGS_OMP ?= -Xcompiler -fopenmp
		LDFLAGS_OMP ?= -Xcompiler -fopenmp
		LDLIBS_OMP ?= -lgomp

		else                # !usehip && !usecuda

#### Use GCC implementation.
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		# LDLIBS_OMP ?= -lgomp

#### Use Intel implementation.
		# CXXFLAGS_OMP ?= -qopenmp
		# LDFLAGS_OMP ?= -qopenmp
		# # LDLIBS_OMP ?= -liomp5

		endif               # !usehip

	else ifeq (${OS}, Darwin)

### Use GCC implementation.
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		# LDLIBS_OMP ?= -lgomp

### Use LLVM implementation from Homebrew (brew formula 'libomp').
		# CXXFLAGS_OMP ?= -I$(shell brew --prefix libomp)/include -Xpreprocessor -fopenmp
		# LDFLAGS_OMP ?= -Wl,-rpath,$(shell brew --prefix libomp)/lib -L$(shell brew --prefix libomp)/lib
		# LDLIBS_OMP ?= -lomp

	else  # OS

### Use GCC implementation.
### If using CUDA, add preprocessing flags.
		ifdef usehip
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		else ifdef usecuda  # !usehip && usecuda
		CXXFLAGS_OMP ?= -Xcompiler -fopenmp
		LDFLAGS_OMP ?= -Xcompiler -fopenmp
		else                # !usehip && !usecuda
		CXXFLAGS_OMP ?= -fopenmp
		LDFLAGS_OMP ?= -fopenmp
		endif               # usehip

	endif  # OS

## If using cuFFT/hipFFT, remove macros for FFTW.
	ifdef usehip
	CPPFLAGS += -DTRV_USE_OMP
	else ifdef usecuda  # !usehip && usecuda
	CPPFLAGS += -DTRV_USE_OMP
	else                # !usehip && !usecuda
	CPPFLAGS += -DTRV_USE_OMP -DTRV_USE_FFTWOMP
	endif               # usehip

	CXXFLAGS += ${CXXFLAGS_OMP}
	LDFLAGS += ${LDFLAGS_OMP}

## If using CUDA FFT, do not include OpenMP FFTW dependency.
	ifdef usehip
	LDLIBS += ${LDLIBS_OMP}
	else ifdef usecuda  # !usehip && usecuda
	LDLIBS += ${LDLIBS_OMP}
	else                # !usehip && !usecuda
	LDLIBS += -lfftw3_omp ${LDLIBS_OMP}
	endif               # usehip

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
else ifdef usecuda  # !usehip && usecuda
CPPFLAGS += -DTRV_USE_CUDA
endif               # usehip

# Visual display
ifdef usedisp
CPPFLAGS += -DTRV_USE_DISP
endif  # usedisp

# Profiling
ifdef useprof
## Linaro MAP profiler
ifdef usehip
CXXFLAGS += -g1 -O3 -fno-inline -fno-optimize-sibling-calls
else ifdef usecuda  # !usehip && usecuda
CXXFLAGS += -g -O3 -lineinfo
else                # !usehip && !usecuda
CXXFLAGS += -g1 -O3 -fno-inline -fno-optimize-sibling-calls
endif               # usehip
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
else ifdef usecuda  # !usehip && usecuda
export PY_CUDA=true
endif               # usehip

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
# Building
# ------------------------------------------------------------------------

SRCS := $(wildcard ${DIR_PKG_SRC}/*.cpp)
ifdef usehip
ifndef usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_hip.o)
else                # usehip && usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_hipcuda.o)
endif               # usehip && !usecuda
else ifdef usecuda  # !usehip && usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%_cuda.o)
else                # !usehip && !usecuda
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%.o)
endif               # usehip
DEPS := $(OBJS:.o=.d)

PROGSRC := ${DIR_PKG_SRCPROG}/${PROGNAME}.cpp
ifdef usehip
ifndef usecuda
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}_hip.o
else                # usehip && usecuda
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}_hipcuda.o
endif               # usehip && !usecuda
else ifdef usecuda  # !usehip && usecuda
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}_cuda.o
else                # !usehip && !usecuda
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}.o
endif               # usehip

PROGEXE := ${DIR_BUILDBIN}/${PROGNAME}
PROGLIB := ${DIR_BUILDLIB}/lib${LIBNAME}.a


# -- Installation --------------------------------------------------------

.PHONY: install cppinstall pyinstall cpplibinstall cppappbuild \
        uninstall cppuninstall pyuninstall

install: cppinstall pyinstall

cppinstall: cppinstall_ cpplibinstall cppappbuild

cppinstall_:
	@echo "Installing ${PKGNAME} C++ library/program..."

cpplibinstall: library

cppappbuild: executable

ifdef usehip
ifdef usecuda
pyinstall:
	@echo "Installing ${PKGNAME} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp .pyproject_hipcuda.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
else   # usehip && !usecuda
pyinstall:
	@echo "Installing ${PKGNAME} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp .pyproject_hip.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
endif  # usehip && usecuda
else ifdef usecuda
pyinstall:
	@echo "Installing ${PKGNAME} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp .pyproject_cuda.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
else  # !usehip && usecuda
pyinstall:
	@echo "Installing ${PKGNAME} Python package ${WOMP} OpenMP (in pip dev mode)..."
	@cp .pyproject.toml pyproject.toml
	python -m pip install ${PIPOPTS} --editable . -vvv
endif  # !usehip && !usecuda

uninstall: cppuninstall pyuninstall

cppuninstall:
	@echo "Uninstalling Triumvirate(-CUDA/HIP/HIPCUDA) C++ library/program..."
	@echo "  removing builds..."
	@find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".git*" -exec rm -r {} +

pyuninstall:
	@echo "Uninstalling ${PKGNAME} Python package (in pip mode)..."
	python -m pip uninstall -y ${PKGNAME}


# -- Components ----------------------------------------------------------

.PHONY: executable library objects_

executable: ${PROGEXE}

${PROGEXE}: $(OBJS) ${PROGOBJ}
	@echo "Compiling ${PKGNAME} C++ program ${WOMP} OpenMP..."
	@if [ ! -d ${DIR_BUILDBIN} ]; then \
	    echo "  making bin subdirectory in build directory..."; \
	    mkdir -p ${DIR_BUILDBIN}; \
	fi
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

library: ${PROGLIB}

${PROGLIB}: $(OBJS)
	@echo "Creating ${PKGNAME} C++ library ${WOMP} OpenMP..."
	@if [ ! -d ${DIR_BUILDLIB} ]; then \
	    echo "  making lib subdirectory in build directory..."; \
	    mkdir -p ${DIR_BUILDLIB}; \
	fi
	$(AR) $(ARFLAGS) $@ $^

objects_:
	@echo "Creating ${PKGNAME} C++ object files..."
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
else                # usehip && !usecuda
$(OBJS): ${DIR_BUILDOBJ}/%_hip.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
endif               # usehip && usecuda
else ifdef usecuda
$(OBJS): ${DIR_BUILDOBJ}/%_cuda.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
else                # !usehip && usecuda
$(OBJS): ${DIR_BUILDOBJ}/%.o: ${DIR_PKG_SRC}/%.cpp | objects_
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@
endif               # !usehip && !usecuda

-include $(DEPS)


# -- Configuration -------------------------------------------------------

.PHONY: checkopts

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
