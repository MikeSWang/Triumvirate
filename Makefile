# -- Configuration ------------------------------------------------------------

SYSTYPE = "cluster"

ifeq ($(SYSTYPE), "local")

CC = mpic++
CFLAGS =# -DDBGNZ -DDBGDK -DTRIUMVIRATE_USE_DISABLED_CODE -Wall

GSL_DIR = /usr/local/gsl
FFTW_DIR = /usr/local/fftw3

INCLUDES = -I./triumvirate/include
INCLUDES += -I${GSL_DIR}/include
INCLUDES += -I${FFTW_DIR}/include

LIBS = -lm
LIBS += -L${GSL_DIR}/lib -lgsl -lgslcblas
LIBS += -L$(FFTW_DIR)/lib -lfftw3

endif

ifeq ($(SYSTYPE), "cluster")

CC = g++
CFLAGS =# -DDBGNZ -DDBGDK -DTRIUMVIRATE_USE_DISABLED_CODE -Wall

INCLUDES = -I./triumvirate/include

LIBS = -lm -lgsl -lgslcblas -lfftw3

endif


# -- Build --------------------------------------------------------------------

## Installation build

install: pyinstall cppinstall

pyinstall:
	pip install -e .

cppinstall: measurements


## Testing build

test: pytest cpptest

pytest:

cpptest: test_common test_parameters test_bessel test_harmonic test_tools \
         test_particles test_field test_twopt test_threept \
				 test_fftlog


## Invididual build

measurements: triumvirate/src/measurements.cpp
	$(CC) $(CFLAGS) -o $(addprefix build/, $(notdir $@)) $^ $(INCLUDES) $(LIBS) $(CLIBS)

test_bessel: triumvirate/tests/test_bessel.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_common: triumvirate/tests/test_common.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_fftlog: triumvirate/tests/test_fftlog.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_field: triumvirate/tests/test_field.cpp
	$(CC) $(CFLAGS) -o $(addprefix triumvirate/tests/test_build/, $(notdir $@)) $< $(INCLUDES) $(LIBS) $(CLIBS)

test_harmonic: triumvirate/tests/test_harmonic.cpp
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
	rm -rf triumvirate/*.cpp triumvirate/*.o triumvirate/*.so
	rm -rf build/*
	rm -rf *.egg-info
	rm -rf **/__pycache__/ core

cleantest:
	rm -rf triumvirate/tests/test_build/* triumvirate/tests/test_output/*
	rm -rf **/.pytest_cache/ core
