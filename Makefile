# *****************
#   Configuration
# *****************

SYSTYPE = "cluster"

ifeq ($(SYSTYPE), "local")

CC = mpic++
CFLAGS = # -Wall

GSL_DIR = /usr/local/gsl
FFTW_DIR = /usr/local/fftw3

INCLUDE = -I./include
INCLUDE += -I${GSL_DIR}/include
INCLUDE += -I${FFTW_DIR}/include

LIB = -lm
LIB += -L${GSL_DIR}/lib -lgsl -lgslcblas
LIB += -L$(FFTW_DIR)/lib -lfftw3

endif

ifeq ($(SYSTYPE), "cluster")

CC = g++
CFLAGS = # -Wall

INCLUDE = -I./include

LIB = -lm -lgsl -lgslcblas -lfftw3

endif


# ************
#   Programs
# ************

all: measurements test_common test_parameters test_bessel

measurements: src/measurements.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_bessel: tests/test_bessel.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_common: tests/test_common.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_parameters: tests/test_parameters.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_particles: tests/test_particles.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_tools: tests/test_tools.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

clean:
	rm -f build/* tests/test_build/* core

distclean:
	rm -f build/* tests/test_build/* *~
