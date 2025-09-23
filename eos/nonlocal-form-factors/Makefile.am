CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosnonlocalformfactors.la
libeosnonlocalformfactors_la_SOURCES = \
	charm-loops.cc charm-loops.hh \
	charm-loops-analytic.cc charm-loops-impl.hh \
	charm-loops-analytic-F279b-helper.cc \
	charm-loops-analytic-F27b.cc \
	charm-loops-analytic-F29b.cc \
	charm-loops-analytic-F27d-helper.cc \
	charm-loops-analytic-F29d-helper.cc \
	charm-loops-analytic-weight4.cc \
	charm-loops-analytic-weight4wx3wx4.cc \
	charm-loops-analytic-weight4w4w5w7.cc \
	charm-loops-analytic-F279d-helper.cc \
	charm-loops-analytic-F27d.cc \
	charm-loops-analytic-F29d.cc \
	hard-scattering.cc hard-scattering.hh \
	long-distance.cc long-distance.hh \
	observables.cc observables.hh \
	nonlocal-formfactors.hh nonlocal-formfactors.cc nonlocal-formfactors-fwd.hh \
	nonlocal-formfactors-p-to-p.cc \
	nonlocal-formfactors-p-to-v.cc
libeosnonlocalformfactors_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosnonlocalformfactors_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosnonlocalformfactors_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lboost_filesystem \
	-lgsl -lgslcblas -lm

include_eos_nonlocalformfactorsdir = $(includedir)/eos/nonlocal-form-factors
include_eos_nonlocalformfactors_HEADERS = \
	charm-loops.hh \
	hard-scattering.hh \
	long-distance.hh \
	nonlocal-formfactors.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	charm-loops_TEST \
	charm-loops-analytic_TEST \
	hard-scattering_TEST \
	long-distance_TEST \
	nonlocal-formfactors-p-to-p_TEST \
	nonlocal-formfactors-p-to-v_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosnonlocalformfactors.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
charm_loops_TEST_SOURCES = charm-loops_TEST.cc

charm_loops_analytic_TEST_SOURCES = charm-loops-analytic_TEST.cc

hard_scattering_TEST_SOURCES = hard-scattering_TEST.cc

long_distance_TEST_SOURCES = long-distance_TEST.cc

nonlocal_formfactors_p_to_p_TEST_SOURCES = nonlocal-formfactors-p-to-p_TEST.cc

nonlocal_formfactors_p_to_v_TEST_SOURCES = nonlocal-formfactors-p-to-v_TEST.cc
