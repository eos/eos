CLEANFILES = \
	*~ \
	manual.acn \
	manual.acr \
	manual.alg \
	manual.aux \
	manual.bbl \
	manual.bcf \
	manual.blg \
	manual.glg \
	manual.glo \
	manual.gls \
	manual.ist \
	manual.out \
	manual.pdf \
	manual.run.xml \
	manual.tdo \
	manual.toc \
	manual.log
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = figures appendices .

TEXSOURCES = \
	acknowledgments.tex \
	manual.tex \
	glossaries.tex \
	doc/installation.tex \
	doc/usage.tex \
	doc/interface.tex \
	appendices/parameters.tex

BIBSOURCES = \
	references.bib

SOURCES = \
	$(BIBSOURCES) $(TEXSOURCES)

manual-am: manual.pdf manual.tex

manual.tex: manual.tex.in
	sed \
	    -e "s/@GITREV@/$$(git describe --tags --long --abbrev=5 --always)/g" \
	    < manual.tex.in \
	    > manual.tex

manual.pdf: $(SOURCES)
	pdflatex manual.tex
	pdflatex manual.tex
	makeglossaries manual
	biber manual
	pdflatex manual.tex
