CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

bin_PROGRAMS = \
	wilson-print-polynomial \
	wilson-scan
wilson_print_polynomial_SOURCES = \
	wilson-print-polynomial.cc
wilson_print_polynomial_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la

wilson_scan_SOURCES = \
	wilson-scan.cc
wilson_scan_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la
