CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in configure config/* aclocal.m4 \
			config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash

SUBDIRS = test eos src doc manual

doxygen:
	$(MAKE) -C doc $@

.PHONY: manual

manual:
	$(MAKE) -C manual $@
