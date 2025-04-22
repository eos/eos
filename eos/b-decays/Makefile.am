CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosbdecays.la
libeosbdecays_la_SOURCES = \
	b-to-d-l-x-nu.cc b-to-d-l-x-nu.hh \
	b-to-gamma-l-nu.cc b-to-gamma-l-nu.hh \
	b-to-l-nu.cc b-to-l-nu.hh \
	b-to-3l-nu.cc b-to-3l-nu.hh \
	b-to-pi-pi-l-nu.cc b-to-pi-pi-l-nu.hh \
	b-to-pi-l-x-nu.cc b-to-pi-l-x-nu.hh \
	b-to-psd-l-nu.cc b-to-psd-l-nu.hh \
	b-to-psd-psd.cc b-to-psd-psd.hh \
	b-to-psd-psd-l-nu.cc b-to-psd-psd-l-nu.hh b-to-psd-psd-l-nu-impl.hh \
	b-to-v-l-nu.hh \
	b-to-vec-l-nu.cc b-to-vec-l-nu.hh b-to-vec-l-nu-impl.hh \
	bq-to-dq-psd.cc bq-to-dq-psd.hh \
	bq-to-dstarq-psd.cc bq-to-dstarq-psd.hh \
	inclusive-b-to-u.cc inclusive-b-to-u.hh \
	lambdab-to-lambdac-l-nu.cc lambdab-to-lambdac-l-nu.hh \
	lambdab-to-lambdac2595-l-nu.cc lambdab-to-lambdac2595-l-nu.hh \
	lambdab-to-lambdac2625-l-nu.cc lambdab-to-lambdac2625-l-nu.hh \
	lifetime.cc lifetime.hh \
	observables.cc observables.hh \
	properties.cc properties.hh \
	signal-pdfs.cc signal-pdfs.hh

libeosbdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosbdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosbdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/nonleptonic-amplitudes/libeosnonleptonicamplitudes.la \
	-lgslcblas \
	-lgsl

include_eos_bdecaysdir = $(includedir)/eos/b-decays
include_eos_bdecays_HEADERS = \
	b-to-d-l-x-nu.hh \
	b-to-gamma-l-nu.hh \
	b-to-l-nu.hh \
	b-to-3l-nu.hh \
	b-to-pi-pi-l-nu.hh \
	b-to-pi-l-x-nu.hh \
	b-to-psd-l-nu.hh \
	b-to-psd-psd.hh \
	b-to-psd-psd-l-nu.hh \
	b-to-v-l-nu.hh \
	b-to-vec-l-nu.hh \
	bq-to-dq-psd.hh \
	bq-to-dstarq-psd.hh \
	inclusive-b-to-u.hh \
	lambdab-to-lambdac-l-nu.hh \
	lifetime.hh \
	observables.hh \
	properties.hh \
	signal-pdfs.hh

EXTRA_DIST =

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	b-to-l-nu_TEST \
	b-to-3l-nu_TEST \
	b-to-d-l-nu_TEST \
	b-to-dstar-l-nu_TEST \
	b-to-eta-l-nu_TEST \
	b-to-gamma-l-nu_TEST \
	b-to-pi-l-nu_TEST \
	b-to-pi-pi-l-nu_TEST \
	b-to-pi-pi-pw-l-nu_TEST \
	b-to-psd-psd_TEST \
	b-to-rho-l-nu_TEST \
	bq-to-dq-psd_TEST \
	bq-to-dstarq-psd_TEST \
	bs-to-kstar-l-nu_TEST \
	lifetime_TEST \
	lambdab-to-lambdac-l-nu_TEST \
	lambdab-to-lambdac2595-l-nu_TEST \
	lambdab-to-lambdac2625-l-nu_TEST \
	properties_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosbdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/nonleptonic-amplitudes/libeosnonleptonicamplitudes.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
b_to_l_nu_TEST_SOURCES = b-to-l-nu_TEST.cc

b_to_3l_nu_TEST_SOURCES = b-to-3l-nu_TEST.cc

b_to_d_l_nu_TEST_SOURCES = b-to-d-l-nu_TEST.cc

b_to_dstar_l_nu_TEST_SOURCES = b-to-dstar-l-nu_TEST.cc

b_to_eta_l_nu_TEST_SOURCES = b-to-eta-l-nu_TEST.cc

b_to_gamma_l_nu_TEST_SOURCES = b-to-gamma-l-nu_TEST.cc

b_to_pi_l_nu_TEST_SOURCES = b-to-pi-l-nu_TEST.cc

b_to_pi_pi_l_nu_TEST_SOURCES = b-to-pi-pi-l-nu_TEST.cc

b_to_pi_pi_pw_l_nu_TEST_SOURCES = b-to-pi-pi-pw-l-nu_TEST.cc

b_to_psd_psd_TEST_SOURCES = b-to-psd-psd_TEST.cc

b_to_rho_l_nu_TEST_SOURCES = b-to-rho-l-nu_TEST.cc

bq_to_dq_psd_TEST_SOURCES = bq-to-dq-psd_TEST.cc

bq_to_dstarq_psd_TEST_SOURCES = bq-to-dstarq-psd_TEST.cc

bs_to_kstar_l_nu_TEST_SOURCES = bs-to-kstar-l-nu_TEST.cc

lifetime_TEST_SOURCES = lifetime_TEST.cc

lambdab_to_lambdac_l_nu_TEST_SOURCES = lambdab-to-lambdac-l-nu_TEST.cc

lambdab_to_lambdac2595_l_nu_TEST_SOURCES = lambdab-to-lambdac2595-l-nu_TEST.cc

lambdab_to_lambdac2625_l_nu_TEST_SOURCES = lambdab-to-lambdac2625-l-nu_TEST.cc

properties_TEST_SOURCES = properties_TEST.cc
