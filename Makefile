# ========================================================================
# Configuration
# ========================================================================

PROGNAME := triumvirate

# -- Directories ---------------------------------------------------------

DIR_ROOT := $(shell pwd)

DIR_PROG = ${DIR_ROOT}/${PROGNAME}
DIR_BUILD = ${DIR_ROOT}/build

DIR_INCLUDE = ${DIR_PROG}/include
DIR_SRC = ${DIR_PROG}/src
DIR_TESTS = ${DIR_PROG}/tests

DIR_SRCMODULES = ${DIR_SRC}/modules

DIR_TESTBUILD = ${DIR_TESTS}/test_build
DIR_TESTOUT = ${DIR_TESTS}/test_output


# -- Common configuration ------------------------------------------------

INCLUDES = -I${DIR_INCLUDE}
LIBS = -lgsl -lgslcblas -lfftw3
CLIBS =


# -- System-specific configuration ---------------------------------------

SYSTYPE := $(if ${NERSC_HOST}, cluster, local)

# CC: compilers, {g++,}
# CFLAGS: compilation flags, {-DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_PARS, -DDBG_NOAC,}
ifeq ($(strip ${SYSTYPE}), local)

CC = g++
CFLAGS =

endif

ifeq ($(strip ${SYSTYPE}), cluster)

CC = g++
CFLAGS =

FFTW_DIR = ${FFTW_ROOT}

endif


# -- Library-specific configuration --------------------------------------

ifdef GSL_DIR

INCLUDES += -I${GSL_DIR}/include
LIBS += -L${GSL_DIR}/lib

endif

ifdef FFTW_DIR

INCLUDES += -I${FFTW_DIR}/include
LIBS += -L${FFTW_DIR}/lib

endif


# ========================================================================
# Build
# ========================================================================

MODULESRC = $(wildcard ${DIR_SRCMODULES}/*.cpp)

# -- Installation build --------------------------------------------------

install: cppinstall pyinstall

cppinstall: ${PROGNAME}

pyinstall:
	echo "${INCLUDES}" > includes.txt
	pip install --user -e .
	rm includes.txt


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
						 ${DIR_SRCMODULES}/monitor.cpp ${DIR_SRCMODULES}/maths.cpp ${DIR_SRCMODULES}/arrayops.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)


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
