CLEANFILES = \
	*~ \
	parameters.tex
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = .

YAMLTEXSOURCES := $(foreach f,$(wildcard $(top_srcdir)/eos/parameters/*.yaml),parameters-$(notdir $f).tex)
TEXSOURCES = \
	parameters.tex.in

SOURCES = \
	$(TEXSOURCES) \
	$(YAMLTEXSOURCES)

manual-am: parameters.tex

parameters-%.yaml.tex: $(top_srcdir)/eos/parameters/%.yaml parameters-to-latex.py
	./parameters-to-latex.py $< > $@

parameters.tex: $(SOURCES) parameters-to-latex.py
	@cp parameters.tex.in parameters.tex.tmp
	@for f in $(YAMLTEXSOURCES) ; do \
	    echo -e "%% from file '$$f'\n\\input{appendices/$$f}\n" \
	        >> parameters.tex.tmp ; \
	done
	@mv parameters.tex.tmp parameters.tex
