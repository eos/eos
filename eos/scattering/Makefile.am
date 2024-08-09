CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeosscattering.la
libeosscattering_la_SOURCES = \
    ee-to-ccbar.cc ee-to-ccbar.hh \
    observables.cc observables.hh \
    parametric-gmkprdey2011.hh parametric-gmkprdey2011.cc \
    scattering-amplitudes.cc scattering-amplitudes.hh scattering-amplitudes-fwd.hh \
    single-channel.cc single-channel.hh \
    single-channel-processes.hh
libeosscattering_la_LIBADD = \
    $(top_builddir)/eos/utils/libeosutils.la \
    -lgslcblas \
    -lgsl

include_eos_scatteringdir = $(includedir)/eos/scattering
include_eos_scattering_HEADERS = \
    ee-to-ccbar.hh \
    parametric-gmkprdey2011.hh \
    scattering-amplitudes.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
    export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
    ee-to-ccbar_TEST \
    parametric-gmkprdey2011_TEST \
    single-channel_TEST

LDADD = \
    $(top_builddir)/test/libeostest.la \
    libeosscattering.la \
    $(top_builddir)/eos/utils/libeosutils.la \
    $(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
ee_to_ccbar_TEST_SOURCES = ee-to-ccbar_TEST.cc

parametric_gmkprdey2011_TEST_SOURCES = parametric-gmkprdey2011_TEST.cc

single_channel_TEST_SOURCES = single-channel_TEST.cc
