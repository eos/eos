GRAPHS = \
	graph-parameters.pdf \
	graph-observable-copy.pdf \
	graph-observable-clone.pdf
CLEANFILES = \
	*~ \
	$(GRAPHS)
MAINTAINERCLEANFILES = Makefile.in

.PHONY: manual

SUFFIXES = .dot .pdf
.dot.pdf:
	dot -Tpdf $< > $@

manual-local: $(GRAPHS)
