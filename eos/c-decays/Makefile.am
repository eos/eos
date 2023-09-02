CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeoscdecays.la
libeoscdecays_la_SOURCES = \
	dq-to-l-nu.cc dq-to-l-nu.hh \
	observables.cc observables.hh
libeoscdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeoscdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeoscdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lgslcblas \
	-lgsl

include_eos_cdecaysdir = $(includedir)/eos/c-decays
include_eos_cdecays_HEADERS = \
	dq-to-l-nu.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	dq-to-l-nu_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeoscdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
dq_to_l_nu_TEST_SOURCES = dq-to-l-nu_TEST.cc
