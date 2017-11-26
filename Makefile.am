CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in configure config/* aclocal.m4 \
			config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash

SUBDIRS = test eos python src doc manual



doxygen:
	$(MAKE) -C doc $@

.PHONY: manual deb

manual:
	$(MAKE) -C manual $@

deb:
	mkdir -p $(DESTDIR)
	make install DESTDIR=$(DESTDIR)
	make -C debian \
	    DESTDIR=$(DESTDIR) \
	    EOS_VERSION=$(VERSION) \
	    FLAVOR=$(shell lsb_release -a | sed -n -e "/^Codename:/s/^Codename:\s*//p") \
	    deb
