# @file docs/Makefile
#
# Makefile for Sphinx documentation
#
# @maintainer Mike S Wang (https://github.com/MikeSWang)
#

# ========================================================================
# Configuration
# ========================================================================

# Define executables explicitly. [MSW]
RM ?= rm -f

# Define build variables.
SPHINXBUILD ?= sphinx-build
SPHINXOPTS  ?=
SOURCEDIR   := source/
BUILDDIR    := build/

# Add separate build directories. [MSW]
BUILDDIR_SPHINX := build/sphinx
BUILDDIR_DOXY   := build/doxygen

# Exclude individual build types. [MSW]
ifdef excl
ifeq ($(strip ${excl}), $(filter $(strip ${excl}), doxy))
export excl_doxy=true
endif  # excl_doxy
ifeq ($(strip ${excl}), $(filter $(strip ${excl}), sphinx))
export excl_sphinx=true
endif  # excl_sphinx
endif  # excl*


# ========================================================================
# Recipes
# ========================================================================

.PHONY: help Makefile clean

# Put it first so that `make` without argument is like ``make help``.
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR_SPHINX)" $(SPHINXOPTS) $(O)

# Route all unknown targets to Sphinx using the new 'make mode' option.
# $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR_SPHINX)" $(SPHINXOPTS) $(O)

# Clean generated documentation artifacts. [MSW]
clean:
	if [ "${excl_doxy}" != "true" ]; then $(RM) -r ${BUILDDIR_DOXY}; fi
	if [ "${excl_sphinx}" != "true" ]; then $(RM) -r ${BUILDDIR_SPHINX}; fi

	find ${BUILDDIR} -mindepth 1 -maxdepth 1 -type d \
	  -not \( -name "doxygen" -or -name "sphinx" \) \
	  -exec $(RM) -r {} +

	for apidir in ${SOURCEDIR}/api*/; do $(RM) -r $$apidir; done
