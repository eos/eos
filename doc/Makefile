# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

doc-check:
	python3 doc_TEST.py

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
html: source Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@test xhtml = x$@ && touch $(BUILDDIR)/html/.nojekyll

# For gradual migration from python scripts to jinja templates:
# explicitly specify here the files that are still scripted.
# The already templated targets are covered automatically further below.
%.rst: %.rst.py %.rst.jinja
	PYTHONPATH="./:$$PYTHONPATH" python3 $< > $@

.PHONY: source
source: reference/constraints.rst reference/observables.rst reference/parameters.rst reference/bibliography.rst reference/signal-pdfs.rst reference/python.rst
