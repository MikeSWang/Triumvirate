# -- Configuration ------------------------------------------------------------

SYSTYPE = "local"

## Common configuration

CFLAGS =# # e.g. -Wall -DTRV_USE_DISABLED_CODE -DDBG_DK

INCLUDES = -I./triumvirate/include
LIBS = -lm -lgsl -lgslcblas -lfftw3

## System-dependent configuration

ifeq (${SYSTYPE}, "local")

CC = g++

GSL_DIR = /usr/local/gsl
FFTW_DIR = /usr/local/fftw3

endif

ifeq (${SYSTYPE}, "cluster") # adapt for different systems (NERSC here)

CC = g++

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

## Installation build

install: cppinstall pyinstall

cppinstall: measurements

pyinstall:
	pip install -e .


## Testing build

test: cpptest pytest

cpptest: test_monitor test_parameters test_bessel test_harmonic test_tools \
         test_particles test_field test_twopt test_threept \
				 test_fftlog

pytest:


## Invididual build

measurements: triumvirate/src/measurements.cpp
	$(CC) $(CFLAGS) -o $(addprefix build/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_bessel: triumvirate/tests/test_bessel.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_fftlog: triumvirate/tests/test_fftlog.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_field: triumvirate/tests/test_field.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_harmonic: triumvirate/tests/test_harmonic.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_monitor: triumvirate/tests/test_monitor.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_parameters: triumvirate/tests/test_parameters.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_particles: triumvirate/tests/test_particles.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_threept: triumvirate/tests/test_threept.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_tools: triumvirate/tests/test_tools.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_twopt: triumvirate/tests/test_twopt.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)


## Build clean-up

clean:
	rm -rf triumvirate/*.cpp triumvirate/*.so triumvirate/*.o
	rm -rf build/* core
	rm -rf *.egg-info
	find . -type d -name "__pycache__" -exec rm -rf {} +

cleantest:
	rm -rf triumvirate/tests/test_build/* triumvirate/tests/test_output/*
	rm -rf core
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
