CLEANFILES =
MAINTAINERCLEANFILES = Makefile.in

EXTRA_DIST = \
    control-focal.in \
    control-jammy.in

# produce the .deb file manually
deb:
	mkdir -p $(DESTDIR)/DEBIAN
	cp control-$(OS) $(DESTDIR)/DEBIAN/control
	dpkg-deb -b $(DESTDIR)
