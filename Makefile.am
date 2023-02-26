CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in configure config/* aclocal.m4 \
			config.h config.h.in
AUTOMAKE_OPTIONS = foreign dist-bzip2
EXTRA_DIST = autogen.bash

# MAYBE_SRC = src by default, and empty when --disable-cli is used
SUBDIRS = test eos python $(MAYBE_SRC) debian



doxygen:
	$(MAKE) -C doc $@

.PHONY: deb

deb:
	mkdir -p $(DESTDIR)
	make install DESTDIR=$(DESTDIR)
	make -C debian deb

print-%:
	@echo $($*)
