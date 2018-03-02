CLEANFILES = \
	*~ \
	eos-list-references

MAINTAINERCLEANFILES = Makefile.in

if EOS_ENABLE_PYTHON

bin_SCRIPTS = \
	eos-list-references \
	eos-plot \
	eos-plot-1d \
	eos-plot-2d \
	eos-print-uncertainty

endif

EXTRA_DIST = \
	eos-list-references.in \
	eos-plot \
	eos-plot-1d \
	eos-plot-2d \
	eos-print-uncertainty

eos-list-references: eos-list-references.in Makefile
	sed \
	    -e "s#@PKGDATADIR@#$(pkgdatadir)#g" \
	    $< \
	    > $@
	chmod a+x $@
