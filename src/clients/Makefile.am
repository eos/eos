CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -pedantic

bin_PROGRAMS = wilson-scan
wilson_scan_SOURCES = \
	wilson-scan.cc
wilson_scan_LDADD = \
	$(top_builddir)/src/rare-b-decays/libwfrarebdecays.la \
	$(top_builddir)/src/utils/libwfutils.la
