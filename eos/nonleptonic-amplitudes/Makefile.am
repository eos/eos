CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosnonleptonicamplitudes.la
libeosnonleptonicamplitudes_la_SOURCES = \
	nonleptonic-amplitudes.cc nonleptonic-amplitudes.hh \
	nonleptonic-amplitudes-fwd.hh nonleptonic-amplitudes-adapter.hh \
	observables.cc observables.hh \
	qcdf-amplitudes.cc qcdf-amplitudes.hh \
	su3f-amplitudes.cc su3f-amplitudes.hh \
	topological-amplitudes.cc topological-amplitudes.hh \
	qcdf_coefs.cc qcdf_coefs.hh

libeosnonleptonicamplitudes_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosnonleptonicamplitudes_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosnonleptonicamplitudes_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	-lgslcblas \
	-lgsl

include_eos_nonleptonicamplitudesdir = $(includedir)/eos/nonleptonic-amplitudes
include_eos_nonleptonicamplitudes_HEADERS = \
	nonleptonic-amplitudes.hh \
	nonleptonic-amplitudes-fwd.hh \
	nonleptonic-amplitudes-adapter.hh \
	observables.hh \
	qcdf-amplitudes.hh \
	su3f-amplitudes.hh \
	topological-amplitudes.hh \
	qcdf_coefs.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	qcdf-amplitudes_TEST \
	su3f-amplitudes_TEST \
	topological-amplitudes_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosnonleptonicamplitudes.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
qcdf_amplitudes_TEST_SOURCES = qcdf-amplitudes_TEST.cc

su3f_amplitudes_TEST_SOURCES = su3f-amplitudes_TEST.cc

topological_amplitudes_TEST_SOURCES = topological-amplitudes_TEST.cc
