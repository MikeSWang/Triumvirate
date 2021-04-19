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
LIB += -lgsl -lgslcblas -L${GSL_DIR}/lib
LIB += -lfftw3 -L$(FFTW_DIR)/lib

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

all: trium test_common test_parameters test_bessel

trium: src/trium.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_common: tests/test_common.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_parameters: tests/test_parameters.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_bessel: tests/test_bessel.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

test_tools: tests/test_tools.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -o $(addprefix tests/test_build/, $(notdir $@)) $< $(CLIBS) $(LIB)

clean:
	rm -f build/* tests/test_build/* core

distclean:
	rm -f build/* tests/test_build/* *~
