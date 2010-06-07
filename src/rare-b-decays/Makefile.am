CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -pedantic

noinst_LTLIBRARIES = libwfrarebdecays.la
libwfrarebdecays_la_SOURCES = \
	bremsstrahlung.cc bremsstrahlung.hh \
	charm-loops.cc charm-loops.hh \
	decays.hh \
	exclusive-b-to-s-dilepton.cc exclusive-b-to-s-dilepton.hh \
	exclusive-b-to-s-gamma.cc exclusive-b-to-s-gamma.hh \
	inclusive-b-to-s-dilepton.cc inclusive-b-to-s-dilepton.hh \
	factory.cc factory.hh \
	form_factors.cc form_factors.hh
libwfrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libwfutils.la -lgslcblas -lgsl

noinst_PROGRAMS = amplitudes ff scan2 scan2-raw scan3 integrated observables
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
ff_SOURCES = \
	ff.cc
ff_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
integrated_SOURCES = \
	integrated.cc
integrated_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
observables_SOURCES = \
	observables.cc
observables_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan2_SOURCES = \
	scan2.cc
scan2_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan2_raw_SOURCES = \
	scan2-raw.cc
scan2_raw_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan3_SOURCES = \
	scan3.cc
scan3_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la

TESTS = \
	bremsstrahlung_TEST \
	charm-loops_TEST \
	exclusive-b-to-s-dilepton_TEST
TESTS_LD_ADD = $(top_builddir)/test/libwftest.a $(top_builddir)/src/utils/libwfutils.la libwfrarebdecays.la

check_PROGRAMS = $(TESTS)
bremsstrahlung_TEST_SOURCES = bremsstrahlung_TEST.cc
bremsstrahlung_TEST_LDADD = $(TESTS_LD_ADD)

charm_loops_TEST_SOURCES = charm-loops_TEST.cc
charm_loops_TEST_LDADD = $(TESTS_LD_ADD)

exclusive_b_to_s_dilepton_TEST_SOURCES = exclusive-b-to-s-dilepton_TEST.cc
exclusive_b_to_s_dilepton_TEST_LDADD = $(TESTS_LD_ADD)
