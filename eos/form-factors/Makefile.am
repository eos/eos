CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeosformfactors.la
libeosformfactors_la_SOURCES = \
	analytic-b-to-pi.hh analytic-b-to-pi.cc \
	analytic-b-to-pi-pi.hh analytic-b-to-pi-pi.cc \
	analytic-b-to-p-lcsr.hh analytic-b-to-p-lcsr-impl.hh \
	analytic-b-to-pi-lcsr.cc  analytic-b-to-k-lcsr.cc  analytic-b-to-d-lcsr.cc \
	analytic-bs-to-k-lcsr.cc  analytic-bs-to-ds-lcsr.cc \
	analytic-b-to-v-lcsr.hh analytic-b-to-v-lcsr-impl.hh \
	analytic-bs-to-kstar-lcsr.cc analytic-bs-to-dsstar-lcsr.cc analytic-bs-to-phi-lcsr.cc \
	analytic-b-to-rho-lcsr.cc analytic-b-to-kstar-lcsr.cc analytic-b-to-dstar-lcsr.cc \
	baryonic.cc baryonic.hh baryonic-impl.hh \
	b-lcdas.cc b-lcdas.hh \
	form-factor-adapter.hh \
	form-factors.cc form-factors.hh form-factors-fwd.hh \
	hqet-b-to-c.hh hqet-b-to-c.cc \
	kstar-lcdas.cc kstar-lcdas.hh \
	mesonic.cc mesonic.hh mesonic-impl.hh \
	mesonic-hqet.hh \
	mesonic-processes.hh \
	observables.cc observables.hh \
	pi-lcdas.cc pi-lcdas.hh \
	unitarity-bounds.cc unitarity-bounds.hh \
	zero-recoil-sum-rule.cc zero-recoil-sum-rule.hh
libeosformfactors_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosformfactors_la_LDFLAGS = $(GSL_LDFLAGS)
libeosformfactors_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	-lgslcblas \
	-lgsl

include_eos_rarebdecaysdir = $(includedir)/eos/form-factors
include_eos_rarebdecays_HEADERS = \
	analytic-b-to-pi.hh \
	analytic-b-to-pi-pi.hh \
	b-lcdas.hh \
	form-factor-adapter.hh \
	form-factors.hh \
	observables.hh \
	pi-lcdas.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	analytic-b-to-pi_TEST \
	analytic-b-to-pi-pi_TEST \
	analytic-b-to-p-lcsr_TEST \
	analytic-b-to-v-lcsr_TEST \
	baryonic_TEST \
	b-lcdas_TEST \
	hqet-b-to-c_TEST \
	kstar-lcdas_TEST \
	mesonic_TEST \
	mesonic-hqet_TEST \
	pi-lcdas_TEST \
	unitarity-bounds_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeosformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

analytic_b_to_pi_TEST_SOURCES = analytic-b-to-pi_TEST.cc

analytic_b_to_pi_pi_TEST_SOURCES = analytic-b-to-pi-pi_TEST.cc
analytic_b_to_pi_pi_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)

analytic_b_to_p_lcsr_TEST_SOURCES = analytic-b-to-p-lcsr_TEST.cc

analytic_b_to_v_lcsr_TEST_SOURCES = analytic-b-to-v-lcsr_TEST.cc

baryonic_TEST_SOURCES = baryonic_TEST.cc

b_lcdas_TEST_SOURCES = b-lcdas_TEST.cc

hqet_b_to_c_TEST_SOURCES = hqet-b-to-c_TEST.cc

kstar_lcdas_TEST_SOURCES = kstar-lcdas_TEST.cc

mesonic_TEST_SOURCES = mesonic_TEST.cc

mesonic_hqet_TEST_SOURCES = mesonic-hqet_TEST.cc

pi_lcdas_TEST_SOURCES = pi-lcdas_TEST.cc

unitarity_bounds_TEST_SOURCES = unitarity-bounds_TEST.cc

