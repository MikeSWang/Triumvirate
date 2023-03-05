# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build/sphinx

# Parse command-line arguments.
ifdef excl
ifeq ($(strip ${excl}), $(filter $(strip ${excl}), doxy))
export excl_doxy=true
endif
ifeq ($(strip ${excl}), $(filter $(strip ${excl}), sphinx))
export excl_sphinx=true
endif
endif

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	if [ "${excl_doxy}" != "true" ]; then rm -rf build/doxygen; fi
	if [ "${excl_sphinx}" != "true" ]; then rm -rf build/sphinx; fi

	find build/ -mindepth 1 -maxdepth 1 -type f \
		-not -name ".gitignore" \
		-exec rm -rf {} +
	find build/ -mindepth 1 -maxdepth 1 -type d \
	  -not \( -name "doxygen" -or -name "sphinx" \) \
		-exec rm -rf {} +

	if [ -d source/apidoc_cpp ]; then \
	  find source/apidoc_cpp -mindepth 1 -exec rm -rf {} +; \
	fi
	if [ -d source/apidoc_py ]; then \
	  find source/apidoc_py -mindepth 1 -exec rm -rf {} +; \
	fi
	if [ -d source/apiref_doxy ]; then \
	  find source/apiref_doxy -mindepth 1 -exec rm -rf {} +; \
	fi
