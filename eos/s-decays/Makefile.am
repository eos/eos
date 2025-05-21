CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeossdecays.la
libeossdecays_la_SOURCES = \
	k-to-l-nu.cc k-to-l-nu.hh \
	k-to-pi-l-nu.cc k-to-pi-l-nu.hh \
	observables.cc observables.hh
libeossdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeossdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeossdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lgslcblas \
	-lgsl

include_eos_sdecaysdir = $(includedir)/eos/s-decays
include_eos_sdecays_HEADERS = \
	k-to-l-nu.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	k-to-l-nu_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeossdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
k_to_l_nu_TEST_SOURCES = k-to-l-nu_TEST.cc
