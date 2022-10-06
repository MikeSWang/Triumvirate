# ========================================================================
# Configuration
# ========================================================================

# -- Common configuration ------------------------------------------------

INCLUDES = -I./triumvirate/include
LIBS = -lgsl -lgslcblas -lfftw3


# -- System-specific configuration ---------------------------------------

# CC: compilers, {g++,}
# CFLAGS: compilation flags, {-DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_PARS, -DDBG_NOAC,}

# Check for the NERSC computing cluster here.
SYSTYPE := $(if ${NERSC_HOST}, cluster, local)

# Adapt for local computing here.
ifeq ($(strip ${SYSTYPE}), local)

CC = g++
CFLAGS =

endif

# Adapt for the NERSC computing cluster here.
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
LIBS += -L$(FFTW_DIR)/lib

endif


# ========================================================================
# Build
# ========================================================================

DIR_INCLUDES=triumvirate/include
DIR_SOURCES=triumvirate/src
DIR_TESTS=triumvirate/tests

DIR_BUILD=build
DIR_TESTBUILD=$(DIR_TESTS)/test_build
DIR_TESTOUT=$(DIR_TESTS)/test_output

SOURCES=$(wildcard $(DIR_SOURCES)/*.cpp)


# -- Installation build --------------------------------------------------

install: cppinstall pyinstall

cppinstall: measurements

pyinstall:
	echo "${INCLUDES}" > includes.txt
	pip install --user -e .
	rm includes.txt


# -- Testing build -------------------------------------------------------

test: cpptest pytest

cpptest: test_fftlog

pytest:


# -- Invididual build ----------------------------------------------------

measurements:
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_BUILD)/, $(notdir $@)) $(SOURCES) $(INCLUDES) $(LIBS) $(CLIBS)

test_fftlog: $(DIR_TESTS)/test_fftlog.cpp \
						 $(DIR_SOURCES)/fftlog.cpp \
						 $(DIR_SOURCES)/monitor.cpp $(DIR_SOURCES)/maths.cpp $(DIR_SOURCES)/arrayops.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)


# ========================================================================
# Clean
# ========================================================================

clean:
	rm -rf triumvirate/*.cpp triumvirate/*.so triumvirate/*.o triumvirate/*.a
	rm -rf core $(DIR_BUILD)/*
	rm -rf *.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

cleantest:
	rm -rf core $(DIR_TESTBUILD)/* $(DIR_TESTOUT)/*
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
