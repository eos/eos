.PHONY: deb control

# produce the DEBIAN/control file
control : control-$(FLAVOR).in
	@echo "Creating DEBIAN/control for target Ubuntu $(FLAVOR)"
	@mkdir -p $(DESTDIR)/DEBIAN
	@sed \
	    -e "s/@EOS_VERSION@/$(EOS_VERSION)/g" \
	    < $< \
	    > $(DESTDIR)/DEBIAN/control

# produce the .deb file manually
deb: control
	@dpkg-deb -b $(DESTDIR)
