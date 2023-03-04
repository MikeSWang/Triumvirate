# ========================================================================
# Configuration
# ========================================================================

PROGNAME := triumvirate
LIBNAME := trv

# ------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------

# Repository root
DIR_ROOT := $(shell pwd)

# Package, build and test directories
DIR_PKG := ${DIR_ROOT}/src/${PROGNAME}
DIR_BUILD := ${DIR_ROOT}/build
DIR_TESTS := ${DIR_ROOT}/tests

# Package subdirectories
DIR_PKG_INCLUDE := ${DIR_PKG}/include
DIR_PKG_SRC := ${DIR_PKG}/src

# Build subdirectories
DIR_BUILDLIB := ${DIR_BUILD}/lib
DIR_BUILDOBJ := ${DIR_BUILD}/obj
DIR_BUILDBIN := ${DIR_BUILD}/bin

# Test subdirectories
DIR_TESTBUILD := ${DIR_TESTS}/test_build
DIR_TESTOUT := ${DIR_TESTS}/test_output

# Package source module subdirectory
DIR_PKG_SRCMODULES := ${DIR_PKG_SRC}/modules


# ------------------------------------------------------------------------
# OS-dependent Compilation
# ------------------------------------------------------------------------

# -- Compiler ------------------------------------------------------------

# Linux: use GNU compiler. [modify]
ifeq ($(shell uname -s), Linux)

CXX := g++

# macOS: use LLVM compiler. [modify]
else ifeq ($(shell uname -s), Darwin)

HOMEBREW_PREFIX ?= /usr/local
CXX := ${HOMEBREW_PREFIX}/opt/llvm/bin/clang++

# Default: use GNU compiler. [modify]
else

CXX := g++

endif


# -- Options -------------------------------------------------------------

INCLUDES := -I${DIR_PKG_INCLUDE}
CFLAGS := -O3 -Wall $(shell pkg-config --cflags gsl fftw3)
LDFLAGS := $(if $(shell pkg-config --libs gsl fftw3),\
                $(shell pkg-config --libs gsl fftw3),\
								$(-lgsl -lgslcblas -lm -lfftw3))


# ------------------------------------------------------------------------
# Environment-specific Settings
# ------------------------------------------------------------------------

# NERSC computer cluster: an example environment. [adapt]
ifdef NERSC_HOST

FFTW_DIR = ${FFTW_ROOT}

# GSL library
ifdef GSL_DIR
INCLUDES += -I${GSL_DIR}/include
LDFLAGS += -L${GSL_DIR}/lib
endif

# FFTW library
ifdef FFTW_DIR
INCLUDES += -I${FFTW_DIR}/include
LDFLAGS += -L${FFTW_DIR}/lib
endif

endif


# ------------------------------------------------------------------------
# Customisation
# ------------------------------------------------------------------------

# OpenMP: enabled with `useomp=[true|1]`; disabled otherwise.
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

CFLAGS += -fopenmp -DTRV_USE_OMP -DTRV_USE_FFTWOMP
LDFLAGS += -lfftw3_omp

endif
endif

# Parameter debugging: enabled with `dbgpars=[true|1]`; disabled otherwise.
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))

CFLAGS += -DDBG_MODE -DDBG_PARS

endif
endif

# Visual enhancements: enabled with `uselogo=[true|1]`; disabled otherwise.
ifdef uselogo
ifeq ($(strip ${uselogo}), $(filter $(strip ${uselogo}), true 1))

CFLAGS += -DTRV_USE_LOGO

endif
endif

# Other options:
# e.g. {-g, -DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_FLAG_NOAC, ...}.
# >>> insert <<<


# ------------------------------------------------------------------------
# Language-specific Settings
# ------------------------------------------------------------------------

# Python: export CXX compilation options as environmental variables.
export PY_CXX=${CXX}
export PY_INCLUDES=${INCLUDES}
export PY_CFLAGS=${CFLAGS}
export PY_LDFLAGS=${LDFLAGS}


# ========================================================================
# Build
# ========================================================================

.PHONY: ${PROGNAME}

PROGEXE := ${PROGNAME}
PROGLIB := ${DIR_BUILDLIB}/lib${LIBNAME}.a

PROGSRC := ${DIR_PKG_SRC}/${PROGEXE}.cpp
MODULESRC := $(wildcard ${DIR_PKG_SRCMODULES}/*.cpp)
PROGOBJ := ${DIR_BUILDOBJ}/${PROGEXE}.o
MODULEOBJ := $(patsubst ${DIR_PKG_SRCMODULES}/%.cpp,${DIR_BUILDOBJ}/%.o,${MODULESRC})

# ------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------

install: cppinstall pyinstall

cppinstall: cpplibinstall cppappbuild

pyinstall:
	@echo "Installing Triumvirate Python package (in development mode)..."
	python -m pip install --verbose --user --editable .

cppappbuild: ${PROGEXE}

cpplibinstall: ${PROGLIB}


# ------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------

unittest: cpptest pytest

cpptest:

pytest:

testit:
	@echo "Performing integration tests... (see ${DIR_TESTOUT}/$@.log)"
	bash ${DIR_TESTS}/$@.sh > ${DIR_TESTOUT}/$@.log


# ------------------------------------------------------------------------
# Components
# ------------------------------------------------------------------------

${PROGEXE}: ${PROGOBJ} ${MODULEOBJ}
	@echo "Compiling Triumvirate C++ program..."
	$(CXX) $(CFLAGS) -o $(addprefix $(DIR_BUILDBIN)/, $(notdir $@)) $^ $(LDFLAGS)

${PROGLIB}: ${MODULEOBJ}
	@echo "Installing Triumvirate C++ library..."
	if [ ! -d build/lib ]; then mkdir -p build/lib; fi
	ar -rcsv build/lib/libtrv.a $^

${PROGOBJ}: ${PROGSRC}
	if [ ! -d build/obj ]; then mkdir -p build/obj; fi
	$(CXX) $(CFLAGS) -o $@ -c $< $(INCLUDES)

${MODULEOBJ}: ${DIR_BUILDOBJ}/%.o: ${DIR_PKG_SRCMODULES}/%.cpp
	if [ ! -d build/obj ]; then mkdir -p build/obj; fi
	$(CXX) $(CFLAGS) -o $@ -c $< $(INCLUDES)


# ========================================================================
# Clean
# ========================================================================

# Ensure deletion safety by limiting the top directory.
DIR_PKG := $(or ${DIR_PKG}, '.')
DIR_BUILD := $(or ${DIR_BUILD}, '.')
DIR_TESTS := $(or ${DIR_TESTS}, '.')
DIR_TESTBUILD := $(or ${DIR_TESTBUILD}, '.')
DIR_TESTOUT := $(or ${DIR_TESTOUT}, '.')

clean:
	@echo "Cleaning up Triumvirate builds..."
	rm -rf *.egg-info
	rm -rf core
	find . -type d -name "__pycache__" -exec rm -rf {} +
	find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} +
	find ${DIR_PKG} -maxdepth 1 \( -name "*.cpp" -or -name "*.so" \) -exec rm {} +
	find ${DIR_BUILD} -mindepth 1 -maxdepth 1 ! -name ".gitignore" -exec rm -r {} +

cleantest:
	@echo "Cleaning up Triumvirate tests..."
	rm -rf ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/* ${DIR_TESTS}/*_temp*
	rm -rf core
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
