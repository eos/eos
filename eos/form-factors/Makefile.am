CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosformfactors.la
libeosformfactors_la_SOURCES = \
	analytic-b-to-pi.hh analytic-b-to-pi.cc \
	baryonic.cc baryonic.hh baryonic-impl.hh \
	form-factor-adapter.hh \
	form-factors.hh form-factors-fwd.hh \
	mesonic.cc mesonic.hh mesonic-impl.hh \
	pi-lcdas.cc pi-lcdas.hh
libeosformfactors_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	-lgslcblas \
	-lgsl

include_eos_rarebdecaysdir = $(includedir)/eos/form-factors
include_eos_rarebdecays_HEADERS = \
	analytic-b-to-pi.hh \
	form-factor-adapter.hh \
	form-factors.hh \
	pi-lcdas.hh

TESTS = \
	analytic-b-to-pi_TEST \
	baryonic_TEST \
	mesonic_TEST \
	pi-lcdas_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
analytic_b_to_pi_TEST_SOURCES = analytic-b-to-pi_TEST.cc

baryonic_TEST_SOURCES = baryonic_TEST.cc

mesonic_TEST_SOURCES = mesonic_TEST.cc

pi_lcdas_TEST_SOURCES = pi-lcdas_TEST.cc
