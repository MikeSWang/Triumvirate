# ========================================================================
# Configuration
# ========================================================================

PROGNAME := triumvirate

# ------------------------------------------------------------------------
# Directories
# ------------------------------------------------------------------------

# Repository root
DIR_ROOT := $(shell pwd)

# Package and build directories
DIR_PKG := ${DIR_ROOT}/${PROGNAME}
DIR_BUILD := ${DIR_ROOT}/build

# Package subdirectories
DIR_INCLUDE := ${DIR_PKG}/include
DIR_SRC := ${DIR_PKG}/src
DIR_TESTS := ${DIR_PKG}/tests

# Build object subdirectory
DIR_BUILDOBJ := ${DIR_BUILD}/obj

# Module source subdirectory
DIR_MODULESRC := ${DIR_SRC}/modules

# Test subdirectories
DIR_TESTBUILD := ${DIR_TESTS}/test_build
DIR_TESTOUT := ${DIR_TESTS}/test_output


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

INCLUDES := -I${DIR_INCLUDE}
CFLAGS := -O3 -Wall $(shell pkg-config --cflags gsl fftw3)
LDFLAGS := $(shell pkg-config --libs gsl fftw3)


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

# OpenMP: enabled with `useomp=[true,1]`; disabled otherwise.
USEOMP = 0
ifdef useomp
ifeq ($(strip ${useomp}), $(filter $(strip ${useomp}), true 1))

CFLAGS += -fopenmp -DTRV_USE_OMP -DTRV_USE_FFTWOMP
LDFLAGS += -lfftw3_omp

USEOMP = 1

endif
endif

# Parameter debugging: enabled with `dbgpars=[true,1]`; disabled otherwise.
DBGPARS = 0
ifdef dbgpars
ifeq ($(strip ${dbgpars}), $(filter $(strip ${dbgpars}), true 1))

CFLAGS += -DDBG_MODE -DDBG_PARS

DBGPARS = 1

endif
endif

# Visual enhancements: enabled with `uselogo=[true,1]`; disabled otherwise.
ifdef uselogo
ifeq ($(strip ${uselogo}), $(filter $(strip ${uselogo}), true 1))

CFLAGS += -DTRV_USE_LOGO

endif
endif

# Other options: e.g. {-g, -DTRV_USE_LEGACY_CODE, -DDBG_MODE, -DDBG_NOAC, ...}.
# >>> insert <<<


# ------------------------------------------------------------------------
# Language-specific Settings
# ------------------------------------------------------------------------

# Python: export CXX compilation options as environmental variables.
export PY_INCLUDES=${INCLUDES}
export PY_LDFLAGS=${LDFLAGS}
export PY_USEOMP=${USEOMP}
export PY_DBGPARS=${DBGPARS}


# ========================================================================
# Build
# ========================================================================

.PHONY: ${PROGNAME}

PROGSRC := ${DIR_SRC}/${PROGNAME}.cpp
MODULESRC := $(wildcard ${DIR_MODULESRC}/*.cpp)
PROGOBJ := ${DIR_BUILDOBJ}/${PROGNAME}.o
MODULEOBJ := $(patsubst ${DIR_MODULESRC}/%.cpp,${DIR_BUILDOBJ}/%.o,${MODULESRC})

# ------------------------------------------------------------------------
# Installation
# ------------------------------------------------------------------------

install: cppinstall pyinstall

cppinstall: ${PROGNAME}

pyinstall:
	@echo "Installing Triumvirate Python package..."
	pip install --user --editable .


# ------------------------------------------------------------------------
# Testing
# ------------------------------------------------------------------------

unittest: cpptest pytest

cpptest:

pytest:

testit:
	@echo "Performing integration tests... (see ${DIR_TESTOUT}/$@.log)"
	@bash ${DIR_TESTS}/$@.sh > ${DIR_TESTOUT}/$@.log


# ------------------------------------------------------------------------
# Components
# ------------------------------------------------------------------------

${PROGNAME}: ${PROGOBJ} ${MODULEOBJ}
	@echo "Building Triumvirate C++ program."
	$(CXX) $(CFLAGS) -o $(addprefix $(DIR_BUILD)/, $(notdir $@)) $^ $(LDFLAGS)

${PROGOBJ}: ${PROGSRC}
	if [ ! -d build/obj ]; then mkdir -p build/obj; fi
	$(CXX) $(CFLAGS) -o $@ -c $< $(INCLUDES)

${MODULEOBJ}: ${DIR_BUILDOBJ}/%.o: ${DIR_MODULESRC}/%.cpp
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
	find ${DIR_PKG} -maxdepth 1 -name "*.cpp" -or -name "*.so" -exec rm -rf {} +
	find ${DIR_BUILD} -mindepth 1 ! -name ".gitignore" -exec rm -rf {} +

cleantest:
	@echo "Cleaning up Triumvirate tests..."
	rm -rf ${DIR_TESTBUILD}/* ${DIR_TESTOUT}/* ${DIR_TESTS}/*_temp*
	rm -rf core
	find . -type d -name ".pytest_cache" -exec rm -rf {} +
