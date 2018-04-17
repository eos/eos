CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in configure config/* aclocal.m4 \
			config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash

SUBDIRS = test eos python src doc manual debian



doxygen:
	$(MAKE) -C doc $@

.PHONY: manual deb

manual:
	$(MAKE) -C manual $@

deb:
	mkdir -p $(DESTDIR)
	make install DESTDIR=$(DESTDIR)
	make -C debian deb

print-%:
	@echo $($*)
