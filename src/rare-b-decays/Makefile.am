CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

noinst_LTLIBRARIES = libeosrarebdecays.la
libeosrarebdecays_la_SOURCES = \
	bremsstrahlung.cc bremsstrahlung.hh \
	charm-loops.cc charm-loops.hh \
	decays.hh \
	exclusive-b-to-s-dilepton.cc exclusive-b-to-s-dilepton.hh \
	exclusive-b-to-s-gamma.cc exclusive-b-to-s-gamma.hh \
	inclusive-b-to-s-dilepton.cc inclusive-b-to-s-dilepton.hh \
	inclusive-b-to-s-gamma.cc inclusive-b-to-s-gamma.hh \
	factory.cc factory.hh \
	form-factors.cc form-factors.hh \
	hard-scattering.cc hard-scattering.hh \
	long-distance.cc long-distance.hh
libeosrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libeosutils.la -lgslcblas -lgsl

noinst_PROGRAMS = amplitudes scan2 scan2-raw scan3 integrated observables
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la
integrated_SOURCES = \
	integrated.cc
integrated_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la
observables_SOURCES = \
	observables.cc
observables_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la
scan2_SOURCES = \
	scan2.cc
scan2_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la
scan2_raw_SOURCES = \
	scan2-raw.cc
scan2_raw_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la
scan3_SOURCES = \
	scan3.cc
scan3_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la

TESTS = \
	bremsstrahlung_TEST \
	charm-loops_TEST \
	exclusive-b-to-s-dilepton_TEST \
	form-factors_TEST \
	inclusive-b-to-s-dilepton_TEST \
	inclusive-b-to-s-gamma_TEST \
	hard-scattering_TEST \
	long-distance_TEST
TESTS_LD_ADD = $(top_builddir)/test/libeostest.a $(top_builddir)/src/utils/libeosutils.la libeosrarebdecays.la

check_PROGRAMS = $(TESTS)
bremsstrahlung_TEST_SOURCES = bremsstrahlung_TEST.cc
bremsstrahlung_TEST_LDADD = $(TESTS_LD_ADD)

charm_loops_TEST_SOURCES = charm-loops_TEST.cc
charm_loops_TEST_LDADD = $(TESTS_LD_ADD)

exclusive_b_to_s_dilepton_TEST_SOURCES = exclusive-b-to-s-dilepton_TEST.cc
exclusive_b_to_s_dilepton_TEST_LDADD = $(TESTS_LD_ADD)

form_factors_TEST_SOURCES = form-factors_TEST.cc
form_factors_TEST_LDADD = $(TESTS_LD_ADD)

inclusive_b_to_s_dilepton_TEST_SOURCES = inclusive-b-to-s-dilepton_TEST.cc
inclusive_b_to_s_dilepton_TEST_LDADD = $(TESTS_LD_ADD)

inclusive_b_to_s_gamma_TEST_SOURCES = inclusive-b-to-s-gamma_TEST.cc
inclusive_b_to_s_gamma_TEST_LDADD = $(TESTS_LD_ADD)

hard_scattering_TEST_SOURCES = hard-scattering_TEST.cc
hard_scattering_TEST_LDADD = $(TESTS_LD_ADD)

long_distance_TEST_SOURCES = long-distance_TEST.cc
long_distance_TEST_LDADD = $(TESTS_LD_ADD)
