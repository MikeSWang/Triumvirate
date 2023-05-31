# @file Makefile
# @brief `make` instructions for building and testing Triumvirate.
# @authors Mike S Wang (https://github.com/MikeSWang)
#

# ========================================================================
# Configuration
# ========================================================================

REPONAME := Triumvirate
PKGNAME := triumvirate
PROGNAME := triumvirate
LIBNAME := trv

SCM_VER_SCHEME ?= no-guess-dev
SCM_LOC_SCHEME ?= node-and-date


# ------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------

# Repository root
DIR_ROOT := $(shell pwd)

# Package, build, test and dist directories
DIR_PKG := ${DIR_ROOT}/src/${PKGNAME}
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

# Extract the '-j' or '--jobs' option (possibly empty).
# Note the space added after ``${MAKEFLAGS}``.
PATTERN_JOBS = "\-j[[:digit:][:space:]]*[^a-z[:punct:]]"
MAKEFLAGS_JOBS = $(shell echo "${MAKEFLAGS} " | grep -Eo ${PATTERN_JOBS})


# ------------------------------------------------------------------------
# Compilation
# ------------------------------------------------------------------------

# -- OS ------------------------------------------------------------------

OS := $(shell uname -s)


# -- Compiler ------------------------------------------------------------

# Assume explicitly GCC compiler by default. [adapt]
ifeq (${OS}, Linux)

CXX ?= g++

else ifeq (${OS}, Darwin)

# Use GCC compiler from Homebrew (brew formula 'gcc').
# The compiler binary may have suffix '-<version>';
# check the version number with ``brew info gcc``.
CXX ?= $(shell find $(brew --prefix gcc)/bin -type f -name 'g++*')

# # Use alternatively LLVM compiler from Homebrew (brew formula 'llvm').
# CXX ?= $(shell brew --prefix llvm)/bin/clang++

else  # OS

CXX ?= g++

endif  # OS

# Assume default achiver. [adapt]
AR ?= ar
ARFLAGS ?= -rcsv

# Assume default remover. [adapt]
RM ?= rm -f


# -- Dependencies --------------------------------------------------------

DEPS := gsl fftw3

DEP_INCLUDES := $(shell pkg-config --cflags-only-I ${DEPS})
DEP_CXXFLAGS := $(shell pkg-config --cflags-only-other ${DEPS})
DEP_LDFLAGS := $(shell pkg-config --libs-only-other --libs-only-L ${DEPS})
DEP_LDLIBS := $(shell pkg-config --libs-only-l ${DEPS})


# -- Options -------------------------------------------------------------

INCLUDES += -I${DIR_PKG_INCLUDE} ${DEP_INCLUDES}
CPPFLAGS += -MMD -MP
CXXFLAGS += -O3 -Wall ${DEP_CXXFLAGS}
LDFLAGS += ${DEP_LDFLAGS}
LDLIBS += $(if ${DEP_LDLIBS},${DEP_LDLIBS},$(-lgsl -lgslcblas -lm -lfftw3))

PIPOPTS ?= --user


# -- Environment ---------------------------------------------------------

# NERSC computer cluster: an example environment [adapt]
ifdef NERSC_HOST

FFTW_DIR = ${FFTW_ROOT}

# GSL library
ifdef GSL_DIR
INCLUDES += -I${GSL_DIR}/include
LDFLAGS += -L${GSL_DIR}/lib
endif  # GSL_DIR

# FFTW library
ifdef FFTW_DIR
INCLUDES += -I${FFTW_DIR}/include
LDFLAGS += -L${FFTW_DIR}/lib
endif  # FFTW_DIR

endif  # NERSC_HOST


# -- Customisation -------------------------------------------------------

# OpenMP: enabled with ``useomp=(true|1)``; disabled otherwise
ifdef useomp

ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

# Assume GCC OpenMP implementation by default. [adapt]
ifeq (${OS}, Linux)

CXXFLAGS_OMP ?= -fopenmp
LDFLAGS_OMP ?= -fopenmp
# LDLIBS_OMP ?= -lgomp

else ifeq (${OS}, Darwin)

CXXFLAGS_OMP ?= -fopenmp
LDFLAGS_OMP ?= -fopenmp
# LDLIBS_OMP ?= -lgomp

# # Use alternatively LLVM OpenMP implementation from Homebrew
# # (brew formula 'libomp').
# CXXFLAGS_OMP ?= -Xpreprocessor -fopenmp
# LDFLAGS_OMP ?= -L$(shell brew --prefix libomp)/lib
# LDLIBS_OMP ?= -lomp

else  # OS

CXXFLAGS_OMP ?= -fopenmp
LDFLAGS_OMP ?= -fopenmp

endif  # OS

CPPFLAGS += -DTRV_USE_OMP -DTRV_USE_FFTWOMP
CXXFLAGS += ${CXXFLAGS_OMP}
LDFLAGS += ${LDFLAGS_OMP}
LDLIBS += -lfftw3_omp ${LDLIBS_OMP}

WOMP := with

else  # useomp!=(true|1)

# NOTE: Use `undefine` for make>=3.82.
unexport useomp

WOMP := without

endif  # useomp==(true|1)

else  # !useomp

WOMP := without

endif  # useomp

# Visual enhancements: enabled with `uselogo=(true|1)`; disabled otherwise
ifdef uselogo
ifeq ($(strip ${uselogo}), $(filter $(strip ${uselogo}), true 1))
CPPFLAGS += -DTRV_USE_LOGO
endif  # uselogo==(true|1)
endif  # uselogo

# Parameter debugging: enabled with `dbgpars=(true|1)`; disabled otherwise
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))
CPPFLAGS += -DDBG_MODE -DDBG_PARS
endif  # dbgpars==(true|1)
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
else  # useomp
export PY_CXXFLAGS_OMP=${CXXFLAGS_OMP}
export PY_LDFLAGS_OMP=${LDFLAGS_OMP}
endif  # !useomp

export PY_BUILD_PARALLEL=${MAKEFLAGS_JOBS}

export PY_SCM_VER_SCHEME=${SCM_VER_SCHEME}
export PY_SCM_LOC_SCHEME=${SCM_LOC_SCHEME}

# C++: strip whitespace.
CPPFLAGS := $(strip ${CPPFLAGS}) $(strip ${INCLUDES})
CXXFLAGS := $(strip ${CXXFLAGS})
LDFLAGS := $(strip ${LDFLAGS})
LDLIBS := $(strip ${LDLIBS})


# ========================================================================
# Recipes
# ========================================================================

# ------------------------------------------------------------------------
# Building
# ------------------------------------------------------------------------

SRCS := $(wildcard ${DIR_PKG_SRC}/*.cpp)
OBJS := $(SRCS:${DIR_PKG_SRC}/%.cpp=${DIR_BUILDOBJ}/%.o)
DEPS := $(OBJS:.o=.d)

PROGSRC := ${DIR_PKG_SRCPROG}/${PROGNAME}.cpp
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}.o

PROGEXE := ${DIR_BUILDBIN}/${PROGNAME}
PROGLIB := ${DIR_BUILDLIB}/lib${LIBNAME}.a


# -- Installation --------------------------------------------------------

.PHONY: install cppinstall pyinstall cpplibinstall cppappbuild \
        uninstall cppuninstall pyuninstall

install: cppinstall pyinstall

cppinstall: cpplibinstall cppappbuild

cpplibinstall: ${PROGLIB}

cppappbuild: ${PROGEXE}

pyinstall:
	@echo "Installing Triumvirate Python package ${WOMP} OpenMP (in pip dev mode)..."
	python -m pip install ${PIPOPTS} --editable . -vvv

uninstall: cppuninstall pyuninstall

cppuninstall:
	@echo "Uninstalling Triumvirate C++ library/program..."
	find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".gitignore" -exec rm -r {} +

pyuninstall:
	@echo "Uninstalling Triumvirate Python package (in pip mode)..."
	python -m pip uninstall -y ${PKGNAME}


# -- Components ----------------------------------------------------------

${PROGEXE}: ${PROGOBJ} $(OBJS)
	@echo "Compiling Triumvirate C++ program ${WOMP} OpenMP..."
	if [ ! -d ${DIR_BUILDBIN} ]; then mkdir -p ${DIR_BUILDBIN}; fi
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS) $(LDLIBS)

${PROGLIB}: $(OBJS)
	@echo "Creating Triumvirate C++ library ${WOMP} OpenMP..."
	if [ ! -d ${DIR_BUILDLIB} ]; then mkdir -p ${DIR_BUILDLIB}; fi
	$(AR) $(ARFLAGS) $@ $^

${PROGOBJ}: ${PROGSRC}
	if [ ! -d ${DIR_BUILDOBJ} ]; then mkdir -p ${DIR_BUILDOBJ}; fi
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(OBJS): ${DIR_BUILDOBJ}/%.o: ${DIR_PKG_SRC}/%.cpp
	if [ ! -d ${DIR_BUILDOBJ} ]; then mkdir -p ${DIR_BUILDOBJ}; fi
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

-include $(DEPS)


# -- Configuration -------------------------------------------------------

.PHONY: checkopts

checkopts:
	@echo "Checking options parsed by Makefile..."
	@echo "MAKEFLAGS: ${MAKEFLAGS}"
	@echo "MAKEFLAGS_JOBS: ${MAKEFLAGS_JOBS}"
	@echo "SCM_VER_SCHEME: ${SCM_VER_SCHEME}"
	@echo "SCM_LOC_SCHEME: ${SCM_LOC_SCHEME}"
	@echo "CXX: ${CXX}"
	@echo "INCLUDES: ${INCLUDES}"
	@echo "CPPFLAGS: ${CPPFLAGS}"
	@echo "CXXFLAGS: ${CXXFLAGS}"
	@echo "LDFLAGS: ${LDFLAGS}"
	@echo "LDLIBS: ${LDLIBS}"
	@echo "CXXFLAGS_OMP: ${CXXFLAGS_OMP}"
	@echo "LDFLAGS_OMP: ${LDFLAGS_OMP}"
	@echo "LDLIBS_OMP: ${LDLIBS_OMP}"
	@echo "AR: ${AR}"
	@echo "ARFLAGS: ${ARFLAGS}"
	@echo "RM: ${RM}"
	@echo "PIPOPTS: ${PIPOPTS}"


# ------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------

.PHONY: test pytest

test: pytest

pytest:
	@echo "Peforming Triumvirate Python tests..."
	if [ ! -d ${DIR_TESTOUT} ]; then mkdir -p ${DIR_TESTOUT}; fi
	pytest


# ------------------------------------------------------------------------
# Cleaning
# ------------------------------------------------------------------------

.PHONY: clean cppclean pyclean testclean distclean

clean: cppclean pyclean testclean distclean

cppclean:
	@echo "Cleaning up Triumvirate C++ build..."
	$(RM) -r core
	find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".gitignore" -exec rm -r {} +

pyclean:
	@echo "Cleaning up Triumvirate Python/Cython build..."
	find ${DIR_PKG} -maxdepth 1 \( -name "*.cpp" -or -name "*.so" \) -exec rm {} +
	find . -type d -name "__pycache__" -exec rm -r {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -r {} +

testclean:
	@echo "Cleaning up Triumvirate tests..."
	$(RM) -r ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/*
	$(RM) -r core
	find . -type d -name ".pytest_cache" -exec rm -r {} +

distclean:
	@echo "Cleaning up Triumvirate distributions..."
	$(RM) -r ${DIR_DIST}/
	$(RM) -r wheelhouse/
	find . -name ".egg-info" -exec rm -r {} +
