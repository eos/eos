CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosformfactors.la
libeosformfactors_la_SOURCES = \
	form-factor-adapter.hh \
	form-factors.hh form-factors-fwd.hh \
	mesonic.cc mesonic.hh mesonic-impl.hh
libeosformfactors_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	-lgslcblas \
	-lgsl

include_eos_rarebdecaysdir = $(includedir)/eos/form-factors
include_eos_rarebdecays_HEADERS = \
	form-factor-adapter.hh \
	form-factors.hh

TESTS = \
	mesonic_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
mesonic_TEST_SOURCES = mesonic_TEST.cc
