# -- Configuration ------------------------------------------------------------

# -- System detection

# Check for the NERSC computing cluster here.
SYSTYPE := $(if ${NERSC_HOST}, cluster, local)

# -- Common configuration

INCLUDES = -I./triumvirate/include
LIBS = -lgsl -lgslcblas -lfftw3

# -- System-dependent configuration

ifeq ($(strip ${SYSTYPE}), local)

CC = g++
CFLAGS = -DDBG_NORMALT# {-DTRV_USE_DISABLED_CODE, -DDBG_PARS, -DDBG_NORMALT, -DDBG_NOAC, -DDBG_BINWIDTH, ...}

endif

# Adapt for the NERSC computing cluster here.
ifeq ($(strip ${SYSTYPE}), cluster)

CC = g++
CFLAGS = -DDBG_NORMALT# {-DTRV_USE_DISABLED_CODE, -DDBG_PARS, -DDBG_NORMALT, -DDBG_NOAC, -DDBG_BINWIDTH, ...}

FFTW_DIR = ${FFTW_ROOT}

endif

ifdef GSL_DIR

INCLUDES += -I${GSL_DIR}/include
LIBS += -L${GSL_DIR}/lib

endif

ifdef FFTW_DIR

INCLUDES += -I${FFTW_DIR}/include
LIBS += -L$(FFTW_DIR)/lib

endif


# -- Build --------------------------------------------------------------------

# Installation build

install: cppinstall pyinstall

cppinstall: measurements

pyinstall:
	echo "${INCLUDES}" > includes.txt
	pip install --user -e .
	rm includes.txt

# Testing build

test: cpptest pytest

cpptest: test_monitor test_parameters test_bessel test_harmonic test_tools \
         test_particles test_field test_twopt test_threept \
				 test_fftlog

pytest:

# Invididual build

DIR_SOURCES=triumvirate/src
DIR_BUILD=build
DIR_TEST=triumvirate/tests
DIR_TESTBUILD=triumvirate/tests/test_build

SOURCES=$(wildcard $(DIR_SOURCES)/*.cpp)

measurements:
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_BUILD)/, $(notdir $@)) $(SOURCES) $(INCLUDES) $(LIBS) $(CLIBS)

test_bessel: $(DIR_TEST)/test_bessel.cpp $(DIR_SOURCES)/bessel.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_fftlog: $(DIR_TEST)/test_fftlog.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_field: $(DIR_TEST)/test_field.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_harmonic: $(DIR_TEST)/test_harmonic.cpp $(DIR_SOURCES)/harmonic.cpp \
							 $(DIR_SOURCES)/monitor.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_monitor: $(DIR_TEST)/test_monitor.cpp $(DIR_SOURCES)/monitor.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_parameters: $(DIR_TEST)/test_parameters.cpp $(DIR_SOURCES)/parameters.cpp \
						     $(DIR_SOURCES)/monitor.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_particles: $(DIR_TEST)/test_particles.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_threept: $(DIR_TEST)/test_threept.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_tools: $(DIR_TEST)/test_tools.cpp $(DIR_SOURCES)/tools.cpp \
            $(DIR_SOURCES)/monitor.cpp $(DIR_SOURCES)/parameters.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_twopt: $(DIR_TEST)/test_twopt.cpp
	$(CC) $(CFLAGS) -o $(addprefix $(DIR_TESTBUILD)/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

# Build clean-up

clean:
	rm -rf triumvirate/*.cpp triumvirate/*.so triumvirate/*.o
	rm -rf build/* core
	rm -rf *.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +

cleantest:
	rm -rf $(DIR_TESTBUILD)/* $(DIR_TEST)/test_output/*
	rm -rf core
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
