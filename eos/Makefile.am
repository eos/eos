CLEANFILES = \
	*~ \
	references_TEST
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = \
	constraints \
	parameters \
	utils \
	maths \
	models \
	statistics \
	form-factors \
	nonlocal-form-factors \
	nonleptonic-amplitudes \
	rare-b-decays \
	b-decays \
	c-decays \
	s-decays \
	meson-mixing \
	scattering \
	.

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	constraint.cc constraint.hh \
	observable.cc observable.hh observable-fwd.hh observable-impl.hh \
	reference.cc reference.hh \
	signal-pdf.cc signal-pdf.hh signal-pdf-fwd.hh signal-pdf-impl.hh
libeos_la_CXXFLAGS = $(AM_CXXFLAGS) \
	-DEOS_DATADIR='"$(datadir)"' \
	$(GSL_CXXFLAGS) \
	$(YAMLCPP_CXXFLAGS)
libeos_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS) $(YAMLCPP_LDFLAGS)
libeos_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/nonlocal-form-factors/libeosnonlocalformfactors.la \
	$(top_builddir)/eos/nonleptonic-amplitudes/libeosnonleptonicamplitudes.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/c-decays/libeoscdecays.la \
	$(top_builddir)/eos/s-decays/libeossdecays.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/meson-mixing/libeosmesonmixing.la \
	$(top_builddir)/eos/scattering/libeosscattering.la

include_eosdir = $(includedir)/eos
include_eos_HEADERS = \
	constraint.hh \
	observable.hh \
	reference.hh \
	signal-pdf.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
	export EOS_TESTS_REFERENCES="$(top_srcdir)/eos/";

TESTS = \
	constraint_TEST \
	observable_TEST \
	reference_TEST \
	signal-pdf_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/c-decays/libeoscdecays.la \
	$(top_builddir)/eos/s-decays/libeossdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/nonleptonic-amplitudes/libeosnonleptonicamplitudes.la \
	$(top_builddir)/eos/nonlocal-form-factors/libeosnonlocalformfactors.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/scattering/libeosscattering.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = \
	constraint_TEST \
	observable_TEST \
	reference_TEST \
	signal-pdf_TEST

constraint_TEST_SOURCES = constraint_TEST.cc
constraint_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
constraint_TEST_LDADD = $(LDADD) -lyaml-cpp

observable_TEST_SOURCES = observable_TEST.cc
observable_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
observable_TEST_LDADD = $(LDADD) -lyaml-cpp

reference_TEST_SOURCES = reference_TEST.cc
reference_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
reference_TEST_LDADD = $(LDADD) -lyaml-cpp

signal_pdf_TEST_SOURCES = signal-pdf_TEST.cc
signal_pdf_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
signal_pdf_TEST_LDADD = $(LDADD) -lyaml-cpp

pkgdata_DATA = references.yaml
EXTRA_DIST = \
	references.yaml
