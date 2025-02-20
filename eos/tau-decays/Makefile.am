CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeostaudecays.la
libeostaudecays_la_SOURCES = \
	tau-to-k-pi-nu.cc tau-to-k-pi-nu.hh \
	observables.cc observables.hh
libeostaudecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeostaudecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeostaudecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lgslcblas \
	-lgsl

include_eos_taudecaysdir = $(includedir)/eos/tau-decays
include_eos_taudecays_HEADERS = \
	tau-to-k-pi-nu.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	tau-to-k-pi-nu_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeostaudecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
tau_to_k_pi_nu_TEST_SOURCES = tau-to-k-pi-nu_TEST.cc
