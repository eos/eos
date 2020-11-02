CLEANFILES = 
MAINTAINERCLEANFILES = Makefile.in

EXTRA_DIST = \
    control-xenial.in \
    control-bionic.in \
    control-disco.in  \
    control-focal.in

# produce the .deb file manually
deb:
	mkdir -p $(DESTDIR)/DEBIAN
	cp control-$(OS) $(DESTDIR)/DEBIAN/control
	dpkg-deb -b $(DESTDIR)
