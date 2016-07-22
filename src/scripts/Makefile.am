CLEANFILES = \
	*~ \
	eos-list-references

MAINTAINERCLEANFILES = Makefile.in

bin_SCRIPTS = \
	eos-list-references

EXTRA_DIST = \
	eos-list-references.in

eos-list-references: eos-list-references.in Makefile
	sed \
	    -e "s#@PKGDATADIR@#$(pkgdatadir)#g" \
	    $< \
	    > $@
	chmod a+x $@
