# @file Makefile
# @brief `make` instructions for building and testing Triumvirate.
# @authors Mike S Wang (https://github.com/MikeSWang)
#

# ========================================================================
# Configuration
# ========================================================================

PROGNAME = triumvirate
LIBNAME = trv

# ------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------

# Repository root
DIR_ROOT := $(shell pwd)

# Package, build and test directories
DIR_PKG := ${DIR_ROOT}/src/${PROGNAME}
DIR_BUILD := ${DIR_ROOT}/build
DIR_TESTS := ${DIR_ROOT}/tests

# Package subdirectories
DIR_PKG_INCLUDE := ${DIR_PKG}/include
DIR_PKG_SRC := ${DIR_PKG}/src

# Build subdirectories
DIR_BUILDLIB := ${DIR_BUILD}/lib
DIR_BUILDOBJ := ${DIR_BUILD}/obj
DIR_BUILDBIN := ${DIR_BUILD}/bin

# Test subdirectories
DIR_TESTBUILD := ${DIR_TESTS}/test_build
DIR_TESTOUT := ${DIR_TESTS}/test_output

# Package source module subdirectory
DIR_PKG_SRCMODULES := ${DIR_PKG_SRC}/modules


# ------------------------------------------------------------------------
# Make Options
# ------------------------------------------------------------------------

# Extract the '-j' or '--jobs' option (possibly empty). Note the space
# added after ``${MAKEFLAGS}``.
MAKEFLAGS_JOBS=$(shell echo "${MAKEFLAGS} " | grep -Eo "\-j[[:digit:][:space:]]*[^a-z[:punct:]]")


# ------------------------------------------------------------------------
# OS-dependent Compilation
# ------------------------------------------------------------------------

# -- Compiler ------------------------------------------------------------

# Assume GCC compiler by default. [modify]
ifeq ($(shell uname -s), Linux)

CXX ?= g++

else ifeq ($(shell uname -s), Darwin)

# Use GCC compiler from Homebrew (brew formula 'gcc').
# The compiler binary provided may have the suffix '-<version>';
# check the version number with `brew info gcc` (here 'g++-11' is assumed).
CXX_DIR = $(shell brew --prefix gcc)/bin
CXX_BIN = $(shell ls ${CXX_DIR} | grep '^g++')
CXX ?= ${CXX_DIR}/${CXX_BIN}

# Others: GCC compiler by default. [modify]
else  # uname -s

CXX ?= g++

endif  # uname -s


# -- Options -------------------------------------------------------------

INCLUDES += -I${DIR_PKG_INCLUDE}
CXXFLAGS += -O3 -Wall $(shell pkg-config --cflags gsl fftw3)
LDFLAGS += $(if $(shell pkg-config --libs gsl fftw3),\
                $(shell pkg-config --libs gsl fftw3),\
                $(-lgsl -lgslcblas -lm -lfftw3))


# ------------------------------------------------------------------------
# Environment-specific Settings
# ------------------------------------------------------------------------

# NERSC computer cluster: an example environment. [adapt]
ifdef NERSC_HOST

FFTW_DIR = ${FFTW_ROOT}

# GSL library
ifdef GSL_DIR
INCLUDES += -I${GSL_DIR}/include
LDFLAGS = -L${GSL_DIR}/lib ${LDFLAGS}
endif  # GSL_DIR

# FFTW library
ifdef FFTW_DIR
INCLUDES += -I${FFTW_DIR}/include
LDFLAGS = -L${FFTW_DIR}/lib ${LDFLAGS}
endif  # FFTW_DIR

endif  # NERSC_HOST


# ------------------------------------------------------------------------
# Customisation
# ------------------------------------------------------------------------

# OpenMP: enabled with `useomp=(true|1)`; disabled otherwise.
ifdef useomp

ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

# Assume GCC OpenMP implementation by default. [modify]
ifeq ($(shell uname -s), Linux)

# Pass...

else ifeq ($(shell uname -s), Darwin)

# For LLVM OpenMP implementation, use 'libomp' from Homebrew
# (brew formula 'libomp'). Set also the following flags.
# CXXFLAGS += -Xpreprocessor
# LDFLAGS_OMP = -L$(shell brew --prefix libomp)/lib -lomp

# Pass...

else  # uname -s

# Pass...

endif  # uname -s

CXXFLAGS += -fopenmp -DTRV_USE_OMP -DTRV_USE_FFTWOMP
LDFLAGS += -lfftw3_omp ${LDFLAGS_OMP}

endif  # useomp=(true|1)

else  # useomp

# NOTE: Use `undefine` for make>=3.82.
unexport useomp

endif  # useomp

ifdef useomp
WOMP=with
else  # useomp
WOMP=without
endif  # useomp

# Visual enhancements: enabled with `uselogo=(true|1)`; disabled otherwise.
ifdef uselogo
ifeq ($(strip ${uselogo}), $(filter $(strip ${uselogo}), true 1))
CXXFLAGS += -DTRV_USE_LOGO
endif  # uselogo=(true|1)
endif  # uselogo

# Parameter debugging: enabled with `dbgpars=(true|1)`; disabled otherwise.
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))
CXXFLAGS += -DDBG_MODE -DDBG_PARS
endif  # dbgpars=(true|1)
endif  # dbgpars

# Other options:
# e.g. {-g, -DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_FLAG_NOAC, ...}.


# ------------------------------------------------------------------------
# Language-specific Settings
# ------------------------------------------------------------------------

# Python: export CXX compilation options as environmental variables.
export PY_CXX=${CXX}
export PY_CXXFLAGS=${CXXFLAGS}
export PY_LDFLAGS=${LDFLAGS}
export PY_INCLUDES=${INCLUDES}
# export PY_OPTS_OMP=${LDFLAGS_OMP}
ifndef useomp
export PY_NO_OMP
endif  # !useomp
export PY_BUILD_PARALLEL=${MAKEFLAGS_JOBS}


# ========================================================================
# Build
# ========================================================================

.PHONY: ${PROGNAME}

PROGEXE := ${PROGNAME}
PROGLIB := ${DIR_BUILDLIB}/lib${LIBNAME}.a

PROGSRC := ${DIR_PKG_SRC}/${PROGEXE}.cpp
MODULESRC := $(wildcard ${DIR_PKG_SRCMODULES}/*.cpp)
PROGOBJ := ${DIR_BUILDOBJ}/${PROGEXE}.o
MODULEOBJ := $(patsubst ${DIR_PKG_SRCMODULES}/%.cpp,${DIR_BUILDOBJ}/%.o,${MODULESRC})

# ------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------

install: cppinstall pyinstall

cppinstall: cpplibinstall cppappbuild

pyinstall:
	@echo "Installing Triumvirate Python package ${WOMP} OpenMP (in development mode)..."
	python -m pip install --user --editable . --verbose

cppappbuild: ${PROGEXE}

cpplibinstall: ${PROGLIB}


# ------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------

test: pytest

pytest:
	@echo "Peforming Triumvirate Python unit tests..."
	if [ ! -d ${DIR_TESTOUT} ]; then mkdir -p ${DIR_TESTOUT}; fi
	pytest


# ------------------------------------------------------------------------
# Components
# ------------------------------------------------------------------------

${PROGEXE}: ${PROGOBJ} ${MODULEOBJ}
	@echo "Compiling Triumvirate C++ program ${WOMP} OpenMP..."
	if [ ! -d build/bin ]; then mkdir -p build/bin; fi
	$(CXX) $(CXXFLAGS) -o $(addprefix $(DIR_BUILDBIN)/, $(notdir $@)) $^ $(LDFLAGS)

${PROGLIB}: ${MODULEOBJ}
	@echo "Building Triumvirate C++ library ${WOMP} OpenMP..."
	if [ ! -d build/lib ]; then mkdir -p build/lib; fi
	ar -rcsv build/lib/libtrv.a $^

${PROGOBJ}: ${PROGSRC}
	if [ ! -d build/obj ]; then mkdir -p build/obj; fi
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<

${MODULEOBJ}: ${DIR_BUILDOBJ}/%.o: ${DIR_PKG_SRCMODULES}/%.cpp
	if [ ! -d build/obj ]; then mkdir -p build/obj; fi
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $@ -c $<


# ========================================================================
# Clean
# ========================================================================

clean: cppclean pyclean

cppclean:
	@echo "Cleaning up Triumvirate C++ build..."
	rm -rf core
	find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".gitignore" -exec rm -r {} +

pyclean:
	@echo "Cleaning up Triumvirate Python/Cython build..."
	rm -rf **/*.egg-info
	find ${DIR_PKG} -maxdepth 1 \( -name "*.cpp" -or -name "*.so" \) -exec rm {} +
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

cleantest:
	@echo "Cleaning up Triumvirate tests..."
	rm -rf core
	rm -rf ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/* ${DIR_TESTS}/*_temp*
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
