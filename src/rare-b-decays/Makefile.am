CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosrarebdecays.la
libeosrarebdecays_la_SOURCES = \
	bremsstrahlung.cc bremsstrahlung.hh \
	charm-loops.cc charm-loops.hh \
	decays.hh \
	em-contributions.hh em-contributions.cc \
	exclusive-b-to-s-dilepton-low-recoil.cc exclusive-b-to-s-dilepton-large-recoil.cc \
	exclusive-b-to-s-dilepton-low-recoil.hh exclusive-b-to-s-dilepton-large-recoil.hh \
	exclusive-b-to-s-gamma.cc exclusive-b-to-s-gamma.hh \
	inclusive-b-to-s-dilepton.cc inclusive-b-to-s-dilepton.hh \
	inclusive-b-to-s-gamma.cc inclusive-b-to-s-gamma.hh \
	form-factors.cc form-factors.hh form-factors-impl.hh \
	hard-scattering.cc hard-scattering.hh \
	long-distance.cc long-distance.hh
libeosrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libeosutils.la -lgslcblas -lgsl

noinst_PROGRAMS = amplitudes
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libeosrarebdecays.la $(top_builddir)/src/utils/libeosutils.la

TESTS = \
	bremsstrahlung_TEST \
	charm-loops_TEST \
	em-contributions_TEST \
	exclusive-b-to-s-dilepton-low-recoil_TEST \
	exclusive-b-to-s-dilepton-large-recoil_TEST \
	exclusive-b-to-s-gamma_TEST \
	form-factors_TEST \
	inclusive-b-to-s-dilepton_TEST \
	inclusive-b-to-s-gamma_TEST \
	hard-scattering_TEST \
	long-distance_TEST
TESTS_LD_ADD = \
	$(top_builddir)/test/libeostest.a \
	$(top_builddir)/src/libeos.la \
	$(top_builddir)/src/utils/libeosutils.la \
	libeosrarebdecays.la

check_PROGRAMS = $(TESTS)
bremsstrahlung_TEST_SOURCES = bremsstrahlung_TEST.cc
bremsstrahlung_TEST_LDADD = $(TESTS_LD_ADD)

charm_loops_TEST_SOURCES = charm-loops_TEST.cc
charm_loops_TEST_LDADD = $(TESTS_LD_ADD)

em_contributions_TEST_SOURCES = em-contributions_TEST.cc
em_contributions_TEST_LDADD = $(TESTS_LD_ADD)

exclusive_b_to_s_dilepton_low_recoil_TEST_SOURCES = exclusive-b-to-s-dilepton-low-recoil_TEST.cc
exclusive_b_to_s_dilepton_low_recoil_TEST_LDADD = $(TESTS_LD_ADD)

exclusive_b_to_s_dilepton_large_recoil_TEST_SOURCES = exclusive-b-to-s-dilepton-large-recoil_TEST.cc
exclusive_b_to_s_dilepton_large_recoil_TEST_LDADD = $(TESTS_LD_ADD)

exclusive_b_to_s_gamma_TEST_SOURCES = exclusive-b-to-s-gamma_TEST.cc
exclusive_b_to_s_gamma_TEST_LDADD = $(TESTS_LD_ADD)

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
