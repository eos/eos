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
	how-to-read.tex \
	manual.tex \
	macros.tex \
	glossaries.tex \
	doc/installation.tex \
	doc/usage.tex \
	doc/interface.tex \
	doc/extending.tex \
	appendices/parameters.tex \
	physics/eft.tex

BIBSOURCES = \
	references.bib \
	references-physics.bib

SOURCES = \
	$(BIBSOURCES) $(TEXSOURCES)

manual-am: manual.pdf manual.tex

manual.tex: manual.tex.in Makefile
	sed \
	    -e "s/@EOS_VERSION@/$(PACKAGE_VERSION)/g" \
	    < manual.tex.in \
	    > manual.tex

manual.pdf: $(SOURCES)
	xelatex manual.tex
	xelatex manual.tex
	makeglossaries manual
	biber manual
	xelatex manual.tex
