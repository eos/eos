CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeosmesonmixing.la
libeosmesonmixing_la_SOURCES = \
    bq-mixing.cc bq-mixing.hh \
	observables.cc observables.hh
libeosmesonmixing_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosmesonmixing_la_LDFLAGS = $(GSL_LDFLAGS)
libeosmesonmixing_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	-lgslcblas \
	-lgsl

include_eos_mesonmixingdir = $(includedir)/eos/meson-mixing
include_eos_mesonmixing_HEADERS = \
	bq-mixing.hh \
	observables.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
    bq-mixing_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosmesonmixing.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
bq_mixing_TEST_SOURCES = bq-mixing_TEST.cc
