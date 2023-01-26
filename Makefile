# ========================================================================
# Configuration
# ========================================================================

PROGNAME := triumvirate

# -- Directories ---------------------------------------------------------

DIR_ROOT := $(shell pwd)

# Repository directories.
DIR_PROG = ${DIR_ROOT}/${PROGNAME}
DIR_BUILD = ${DIR_ROOT}/build

# Package directories.
DIR_INCLUDE = ${DIR_PROG}/include
DIR_SRC = ${DIR_PROG}/src
DIR_TESTS = ${DIR_PROG}/tests

# Module directories.
DIR_SRCMODULES = ${DIR_SRC}/modules

# Test directories.
DIR_TESTBUILD = ${DIR_TESTS}/test_build
DIR_TESTOUT = ${DIR_TESTS}/test_output


# -- Common configuration ------------------------------------------------

ifeq ($(shell uname -s),Darwin)
# Here we use LLVM compiler to load the Mac OpenMP
ifndef HOMEBREW_PREFIX
HOMEBREW_PREFIX = /usr/local
endif
CC = ${HOMEBREW_PREFIX}/opt/llvm/bin/clang++
else
# default (Linux) case
CC = g++
endif
INCLUDES = -I${DIR_INCLUDE}
CFLAGS = -O3 -Wall $(shell pkg-config --cflags gsl fftw3)
LIBS = $(shell pkg-config --libs gsl fftw3)
CLIBS =


# -- System-specific configuration ---------------------------------------

# This sub-section is specific to NERSC clusters.

SYSTYPE := $(if ${NERSC_HOST}, cluster, local)

# Non-NERSC systems.
ifeq ($(strip ${SYSTYPE}), local)

# >>>>>>>>>>
# (edit/insert here)
# <<<<<<<<<<

endif

# NERSC system.
ifeq ($(strip ${SYSTYPE}), cluster)

# >>>>>>>>>>
# (edit/insert here)
# <<<<<<<<<<

FFTW_DIR = ${FFTW_ROOT}

endif


# -- Library-specific configuration --------------------------------------

# GSL library.
ifdef GSL_DIR

INCLUDES += -I${GSL_DIR}/include
LIBS += -L${GSL_DIR}/lib

endif

# FFTW library.
ifdef FFTW_DIR

INCLUDES += -I${FFTW_DIR}/include
LIBS += -L${FFTW_DIR}/lib

endif


# -- Compilation-specific configurations ---------------------------------

export PY_INCLUDES=${INCLUDES}
export PY_LIBS=${LIBS}

# Enable OpenMP by setting `useomp=true` or `useomp=1`, which adds
# `-fopenmp -DTRV_USE_OMP` and `-lfftw3_omp`.
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

CFLAGS += -fopenmp -DTRV_USE_OMP -DTRV_USE_FFTWOMP
LIBS += -lfftw3_omp

export PY_USEOMP=1

endif
endif

# Enable parameter debugging by setting `dbgpars=true` or `dbgpars=1`.
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))

CFLAGS += -DDBG_MODE -DDBG_PARS

export PY_DBGPARS=1

endif
endif

# Enable visual enhancements.
ifdef uselogo
ifeq ($(strip ${uselogo}), $(filter $(strip ${uselogo}), true 1))

CFLAGS += -DTRV_USE_LOGO

endif
endif

# Add other options (for developers only), e.g.
# {-g, -DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_NOAC, ...}.
# >>>>>>>>>>
CFLAGS +=
# <<<<<<<<<<


# ========================================================================
# Build
# ========================================================================

.PHONY: ${PROGNAME}

MODULESRC = $(wildcard ${DIR_SRCMODULES}/*.cpp)

# -- Installation build --------------------------------------------------

install: cppinstall pyinstall

cppinstall: ${PROGNAME}

pyinstall:
	@echo "Installing Triumvirate Python package."
	pip install --user --editable .


# -- Testing build -------------------------------------------------------

test: cpptest pytest

testit:
	@echo "Performing integration tests. See ${DIR_TESTOUT}/$@.log for log."
	@bash ${DIR_TESTS}/$@.sh > ${DIR_TESTOUT}/$@.log

cpptest: test_fftlog

pytest:


# -- Invididual build ----------------------------------------------------

${PROGNAME}: ${DIR_SRC}/${PROGNAME}.cpp
	@echo "Building Triumvirate C++ program."
	$(CC) $(CFLAGS) \
	-o $(addprefix $(DIR_BUILD)/, $(notdir $@)) \
	$< $(MODULESRC) $(INCLUDES) $(LIBS) $(CLIBS)

test_fftlog: ${DIR_TESTS}/test_fftlog.cpp \
						 ${DIR_SRCMODULES}/fftlog.cpp \
						 ${DIR_SRCMODULES}/monitor.cpp \
						 ${DIR_SRCMODULES}/maths.cpp \
						 ${DIR_SRCMODULES}/arrayops.cpp
	$(CC) $(CFLAGS) \
	-o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) \
	$^ $(INCLUDES) $(LIBS) $(CLIBS)


# ========================================================================
# Clean
# ========================================================================

DIR_PROG := $(or ${DIR_PROG}, '.')
DIR_BUILD := $(or ${DIR_BUILD}, '.')
DIR_TESTS := $(or ${DIR_TESTS}, '.')
DIR_TESTBUILD := $(or ${DIR_TESTBUILD}, '.')
DIR_TESTOUT := $(or ${DIR_TESTOUT}, '.')

clean:
	@echo "Cleaning up Triumvirate builds."
	rm -rf ${DIR_PROG}/*.cpp ${DIR_PROG}/*.o ${DIR_PROG}/*.so ${DIR_BUILD}/*
	rm -rf *.egg-info
	rm -rf core
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

cleantest:
	@echo "Cleaning up Triumvirate tests."
	rm -rf ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/* ${DIR_TESTS}/*_temp*
	rm -rf core
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
