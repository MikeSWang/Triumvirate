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

INCLUDES = -I${DIR_INCLUDE}
LIBS = -lgsl -lgslcblas -lfftw3
CLIBS =


# -- System-specific configuration ---------------------------------------

# This sub-section is specific to NERSC clusters.

SYSTYPE := $(if ${NERSC_HOST}, cluster, local)

# Non-NERSC systems.
ifeq ($(strip ${SYSTYPE}), local)

# >>>>>>>>>>
CC = g++
CFLAGS =
# <<<<<<<<<<

endif

# NERSC system.
ifeq ($(strip ${SYSTYPE}), cluster)

# >>>>>>>>>>
CC = g++
CFLAGS =
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

# Enable OpenMP by setting `useomp=true` or `useomp=1`, which adds
# `-fopenmp -DTRV_USE_OMP` and `-lfftw3_omp`.
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

CFLAGS += -fopenmp -DTRV_USE_OMP
LIBS += -lfftw3_omp

PY_USEOMP = 1

endif
endif

# Enable parameter debugging by setting `dbgpars=true` or `dbgpars=1`.
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))

CFLAGS += -DDBG_MODE -DDBG_PARS

PY_DBGPARS = 1

endif
endif

# Add other options (for developers only), e.g.
# {-DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_NOAC, ...}.
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
	export PY_INCLUDES="${INCLUDES}"
	export PY_USEOMP
	export PY_DBGPARS
	pip install --user -e .


# -- Testing build -------------------------------------------------------

test: cpptest pytest

cpptest: test_fftlog

pytest:


# -- Invididual build ----------------------------------------------------

${PROGNAME}: ${DIR_SRC}/${PROGNAME}.cpp
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

clean:
	rm -rf ${DIR_PROG}/*.cpp ${DIR_PROG}/*.o ${DIR_PROG}/*.so
	rm -rf core ${DIR_BUILD}/*
	rm -rf *.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

cleantest:
	rm -rf core ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/*
	find . -type d -name ".pytest_cache" -exec rm -rf {} +