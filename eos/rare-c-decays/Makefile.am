CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosrarecdecays.la
libeosrarecdecays_la_SOURCES = \
	d-to-psd-l-l-base.cc d-to-psd-l-l-base.hh \
	d-to-psd-l-l-bght2019.cc d-to-psd-l-l-bght2019.hh \
	d-to-psd-l-l.cc d-to-psd-l-l.hh \
	lambda-c-to-proton-l-l.cc lambda-c-to-proton-l-l.hh \
	observables.cc observables.hh
libeosrarecdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosrarecdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosrarecdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lboost_filesystem \
	-lgsl -lgslcblas -lm

include_eos_rarecdecaysdir = $(includedir)/eos/rare-c-decays
include_eos_rarecdecays_HEADERS = \
	d-to-psd-l-l-base.hh \
	d-to-psd-l-l-bght2019.hh \
	d-to-psd-l-l.hh \
	lambda-c-to-proton-l-l.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	d-to-psd-l-l-bght2019_TEST
	lambda-c-to-proton-l-l_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosrarecdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
d_to_psd_l_l_bght2019_TEST_SOURCES = d-to-psd-l-l-bght2019_TEST.cc
lambda_c_to_proton_l_l_TEST_SOURCES = lambda-c-to-proton-l-l_TEST.cc
