CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosbdecays.la
libeosbdecays_la_SOURCES = \
	d-to-psd-l-nu.cc d-to-psd-l-nu.hh

libeosbdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosbdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosbdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lgslcblas \
	-lgsl

include_eos_bdecaysdir = $(includedir)/eos/d-decays
include_eos_bdecays_HEADERS = \
	d-to-psd-l-nu.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = 

LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeosbdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
