CLEANFILES = \
	*~ \
	parameters.tex
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = .

PARAMSOURCES := $(foreach f,$(wildcard $(top_srcdir)/eos/parameters/*.yaml),parameters-$(notdir $f).tex)

manual-am: parameters.tex observables.tex

parameters-%.yaml.tex: $(top_srcdir)/eos/parameters/%.yaml parameters-to-latex.py
	./parameters-to-latex.py $< > $@

parameters.tex: parameters.tex.in $(PARAMSOURCES) parameters-to-latex.py
	@cp parameters.tex.in parameters.tex.tmp
	@for f in $(PARAMSOURCES) ; do \
	    echo -e "%% from file '$$f'\n\\input{appendices/$$f}\n" \
	        >> parameters.tex.tmp ; \
	done
	@mv parameters.tex.tmp parameters.tex

observables.tex: observables.tex.in observables-to-latex.py
	@cp observables.tex.in  observables.tex.tmp
	./observables-to-latex.py >> observables.tex.tmp
	@mv observables.tex.tmp observables.tex
