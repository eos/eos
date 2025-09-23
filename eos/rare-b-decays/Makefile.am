CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosrarebdecays.la
libeosrarebdecays_la_SOURCES = \
	b-to-k-charmonium.cc b-to-k-charmonium.hh \
	b-to-k-ll.cc b-to-k-ll.hh \
	b-to-k-ll-base.cc b-to-k-ll-base.hh \
	b-to-k-ll-bfs2004.cc b-to-k-ll-bfs2004.hh \
	b-to-k-ll-gp2004.cc b-to-k-ll-gp2004.hh \
	b-to-k-ll-gvdv2020.cc b-to-k-ll-gvdv2020.hh \
	b-to-k-ll-naive.cc b-to-k-ll-naive.hh \
	b-to-kstar-charmonium.cc b-to-kstar-charmonium.hh \
	b-to-kstar-gamma.cc b-to-kstar-gamma.hh \
	b-to-kstar-gamma-base.cc b-to-kstar-gamma-base.hh \
	b-to-kstar-gamma-bfs2004.cc b-to-kstar-gamma-bfs2004.hh \
	b-to-kstar-ll.cc b-to-kstar-ll.hh b-to-kstar-ll-impl.hh \
	b-to-kstar-ll-base.cc b-to-kstar-ll-base.hh \
	b-to-kstar-ll-bfs2004.cc b-to-kstar-ll-bfs2004.hh \
	b-to-kstar-ll-gp2004.cc b-to-kstar-ll-gp2004.hh \
	b-to-kstar-ll-gvdv2020.cc b-to-kstar-ll-gvdv2020.hh \
	b-to-kstar-ll-naive.cc b-to-kstar-ll-naive.hh \
	b-to-ll.cc b-to-ll.hh \
	b-to-psd-nu-nu.cc b-to-psd-nu-nu.hh \
	b-to-vec-nu-nu.cc b-to-vec-nu-nu.hh \
	bremsstrahlung.cc bremsstrahlung.hh \
	bs-to-phi-charmonium.cc bs-to-phi-charmonium.hh \
	bs-to-phi-ll.cc bs-to-phi-ll.hh \
	bs-to-phi-ll-base.cc bs-to-phi-ll-base.hh \
	bs-to-phi-ll-bfs2004.cc bs-to-phi-ll-bfs2004.hh \
	bs-to-phi-ll-gvdv2020.cc bs-to-phi-ll-gvdv2020.hh \
	bs-to-phi-ll-naive.cc bs-to-phi-ll-naive.hh \
	decays.hh \
	em-contributions.hh em-contributions.cc \
	inclusive-b-to-s-dilepton.cc inclusive-b-to-s-dilepton.hh \
	inclusive-b-to-s-gamma.cc inclusive-b-to-s-gamma.hh \
	lambda-b-to-lambda-dilepton.cc lambda-b-to-lambda-dilepton.hh \
	lambda-b-to-lambda-nu-nu.cc lambda-b-to-lambda-nu-nu.hh lambda-b-to-lambda-nu-nu-impl.hh\
	lambda-b-to-lambda1520-gamma.cc lambda-b-to-lambda1520-gamma.hh \
	lambda-b-to-lambda1520-gamma-base.cc lambda-b-to-lambda1520-gamma-base.hh \
	lambda-b-to-lambda1520-gamma-naive.cc lambda-b-to-lambda1520-gamma-naive.hh \
	lambda-b-to-lambda1520-ll.cc lambda-b-to-lambda1520-ll.hh \
	lambda-b-to-lambda1520-ll-base.cc lambda-b-to-lambda1520-ll-base.hh \
	lambda-b-to-lambda1520-ll-naive.cc lambda-b-to-lambda1520-ll-naive.hh \
	observables.cc observables.hh \
	qcdf-integrals.cc qcdf-integrals.hh qcdf-integrals-impl.hh\
	qcdf-integrals-analytical.cc \
	qcdf-integrals-mixed.cc \
	qcdf-integrals-numerical.cc \
	signal-pdfs.cc signal-pdfs.hh
libeosrarebdecays_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosrarebdecays_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosrarebdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/nonlocal-form-factors/libeosnonlocalformfactors.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	-lboost_filesystem \
	-lgsl -lgslcblas -lm

include_eos_rarebdecaysdir = $(includedir)/eos/rare-b-decays
include_eos_rarebdecays_HEADERS = \
	b-to-psd-nu-nu.hh \
	b-to-vec-nu-nu.hh \
	bremsstrahlung.hh \
	decays.hh \
	em-contributions.hh \
	lambda-b-to-lambda-dilepton.hh \
	lambda-b-to-lambda-nu-nu.hh \
	qcdf-integrals.hh

EXTRA_DIST = \
	exclusive-b-to-s-dilepton-low-recoil_TEST-btokll.data \
	exclusive-b-to-s-dilepton-low-recoil_TEST-btokstarll.data \
	exclusive-b-to-s-gamma_TEST-btokstargamma.data

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	bremsstrahlung_TEST \
	b-to-k-charmonium_TEST \
	b-to-k-ll-bfs2004_TEST \
	b-to-k-ll-gp2004_TEST \
	b-to-k-ll-gvdv2020_TEST \
	b-to-k-nu-nu_TEST \
	b-to-kstar-charmonium_TEST \
	b-to-kstar-gamma-bfs2004_TEST \
	b-to-kstar-ll-bfs2004_TEST \
	b-to-kstar-ll-gp2004_TEST \
	b-to-kstar-ll-gvdv2020_TEST \
	b-to-kstar-nu-nu_TEST \
	b-to-ll_TEST \
	bs-to-phi-charmonium_TEST \
	bs-to-phi-ll-bfs2004_TEST \
	bs-to-phi-ll-gvdv2020_TEST \
	em-contributions_TEST \
	inclusive-b-to-s-dilepton_TEST \
	inclusive-b-to-s-gamma_TEST \
	lambda-b-to-lambda-dilepton_TEST \
	lambda-b-to-lambda1520-gamma-naive_TEST \
	lambda-b-to-lambda1520-ll-naive_TEST \
	qcdf-integrals_TEST

LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosrarebdecays.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/nonlocal-form-factors/libeosnonlocalformfactors.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
bremsstrahlung_TEST_SOURCES = bremsstrahlung_TEST.cc

b_to_k_charmonium_TEST_SOURCES = b-to-k-charmonium_TEST.cc

b_to_k_ll_bfs2004_TEST_SOURCES = b-to-k-ll-bfs2004_TEST.cc

b_to_k_ll_gp2004_TEST_SOURCES = b-to-k-ll-gp2004_TEST.cc

b_to_k_ll_gvdv2020_TEST_SOURCES = b-to-k-ll-gvdv2020_TEST.cc

b_to_k_nu_nu_TEST_SOURCES = b-to-k-nu-nu_TEST.cc

b_to_kstar_charmonium_TEST_SOURCES = b-to-kstar-charmonium_TEST.cc

b_to_kstar_gamma_bfs2004_TEST_SOURCES = b-to-kstar-gamma-bfs2004_TEST.cc

b_to_kstar_ll_bfs2004_TEST_SOURCES = b-to-kstar-ll-bfs2004_TEST.cc

b_to_kstar_ll_gp2004_TEST_SOURCES = b-to-kstar-ll-gp2004_TEST.cc

b_to_kstar_ll_gvdv2020_TEST_SOURCES = b-to-kstar-ll-gvdv2020_TEST.cc

b_to_kstar_nu_nu_TEST_SOURCES = b-to-kstar-nu-nu_TEST.cc

b_to_ll_TEST_SOURCES = b-to-ll_TEST.cc

bs_to_phi_charmonium_TEST_SOURCES = bs-to-phi-charmonium_TEST.cc

bs_to_phi_ll_bfs2004_TEST_SOURCES = bs-to-phi-ll-bfs2004_TEST.cc

bs_to_phi_ll_gvdv2020_TEST_SOURCES = bs-to-phi-ll-gvdv2020_TEST.cc

em_contributions_TEST_SOURCES = em-contributions_TEST.cc

inclusive_b_to_s_dilepton_TEST_SOURCES = inclusive-b-to-s-dilepton_TEST.cc

inclusive_b_to_s_gamma_TEST_SOURCES = inclusive-b-to-s-gamma_TEST.cc

lambda_b_to_lambda_dilepton_TEST_SOURCES = lambda-b-to-lambda-dilepton_TEST.cc

lambda_b_to_lambda1520_gamma_naive_TEST_SOURCES = lambda-b-to-lambda1520-gamma-naive_TEST.cc

lambda_b_to_lambda1520_ll_naive_TEST_SOURCES = lambda-b-to-lambda1520-ll-naive_TEST.cc

qcdf_integrals_TEST_SOURCES = qcdf-integrals_TEST.cc
