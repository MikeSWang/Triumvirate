# *****************
#   Configuration
# *****************

SYSTYPE = "cluster"

ifeq ($(SYSTYPE), "local")

CC = mpic++
CFLAGS = # -Wall

INCLUDE = -I./triumvirate/measurements/include -I./triumvirate/modelling/include

LIB = -lm

GSL_DIR = /usr/local/gsl
FFTW_DIR = /usr/local/fftw3

INCLUDE += -I${GSL_DIR}/include
INCLUDE += -I${FFTW_DIR}/include

LIB += -L${GSL_DIR}/lib -lgsl -lgslcblas
LIB += -L$(FFTW_DIR)/lib -lfftw3

endif

ifeq ($(SYSTYPE), "cluster")

CC = g++
CFLAGS = # -Wall

INCLUDE = -I./triumvirate/measurements/include -I./triumvirate/modelling/include

LIB = -lm -lgsl -lgslcblas -lfftw3

endif


# ************
#   Programs
# ************

all: measurements modelling test_common test_parameters test_bessel

measurements: triumvirate/measurements/src/measurements.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix build/, $(notdir $@)) $< $(CLIBS) $(LIB)

modelling: triumvirate/setup.py
	python $< install --user

test_bessel: triumvirate/tests/test_bessel.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_common: triumvirate/tests/test_common.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_parameters: triumvirate/tests/test_parameters.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_particles: triumvirate/tests/test_particles.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_tools: triumvirate/tests/test_tools.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

clean:
	rm -f build/* triumvirate/tests/test_build/* core

distclean:
	rm -f build/* triumvirate/tests/test_build/* *~

testclean:
	find triumvirate/tests/test_output ! -name '.gitignore' -type f -exec rm -f {} +
