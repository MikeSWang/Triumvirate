# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	find build/ -mindepth 1 -maxdepth 1 \
	  -not -name ".gitignore"\
	  -exec rm -rf {} +
	find source/ -type f\
	  -not -name "index.rst"\
	  -not -name "installation.rst"\
	  -not -name "tutorials.rst"\
	  -not -name "apidoc.rst"\
	  -iname "*.rst"\
	  -exec rm -rf {} +
