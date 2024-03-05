CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosformfactors.la
libeosformfactors_la_SOURCES = \
	analytic-b-to-gamma-qcdf.cc analytic-b-to-gamma-qcdf.hh \
	analytic-b-to-psd-dkmmo2008.hh analytic-b-to-psd-dkmmo2008-impl.hh \
	analytic-b-to-pi-dkmmo2008.cc \
	analytic-bs-to-k-dkmmo2008.cc \
	analytic-b-to-pi-pi.hh analytic-b-to-pi-pi.cc \
	analytic-b-to-p-lcsr.hh analytic-b-to-p-lcsr-impl.hh \
	analytic-b-to-pi-lcsr.cc  analytic-b-to-k-lcsr.cc  analytic-b-to-d-lcsr.cc \
	analytic-bs-to-k-lcsr.cc  analytic-bs-to-ds-lcsr.cc \
	analytic-b-to-v-lcsr.hh analytic-b-to-v-lcsr-impl.hh \
	analytic-bs-to-kstar-lcsr.cc analytic-bs-to-dsstar-lcsr.cc analytic-bs-to-phi-lcsr.cc \
	analytic-b-to-rho-lcsr.cc analytic-b-to-kstar-lcsr.cc analytic-b-to-dstar-lcsr.cc \
	baryonic.cc baryonic.hh \
	baryonic-processes.hh \
	form-factor-adapter.hh \
	form-factors.cc form-factors.hh form-factors-fwd.hh \
	heavy-meson-lcdas.cc heavy-meson-lcdas.hh \
	heavy-meson-lcdas-exponential.cc heavy-meson-lcdas-exponential.hh \
	heavy-meson-lcdas-flvd2022.cc heavy-meson-lcdas-flvd2022.hh \
	hqet-b-to-c.hh hqet-b-to-c.cc \
	k-lcdas.cc k-lcdas.hh \
	k-star-lcdas.cc k-star-lcdas.hh \
	lcda.cc lcda.hh \
	mesonic.cc mesonic.hh \
	mesonic-processes.hh \
	observables.cc observables.hh \
	parametric-abr2022.hh parametric-abr2022-impl.hh parametric-abr2022.cc \
	parametric-bbgorvd2018.hh parametric-bbgorvd2018-impl.hh parametric-bbgorvd2018.cc \
	parametric-bcl2008.hh parametric-bcl2008-impl.hh parametric-bcl2008.cc \
	parametric-bfvd2014.hh parametric-bfvd2014.cc \
	parametric-bfw2010.hh parametric-bfw2010-impl.hh parametric-bfw2010.cc \
	parametric-bgl1997.hh parametric-bgl1997-impl.hh parametric-bgl1997.cc \
	parametric-bgjvd2019.hh parametric-bgjvd2019-impl.hh parametric-bgjvd2019.cc \
	parametric-bmrvd2022.hh parametric-bmrvd2022-impl.hh parametric-bmrvd2022.cc \
	parametric-bsz2015.hh parametric-bsz2015-impl.hh parametric-bsz2015.cc \
	parametric-dkmr2017.hh parametric-dkmr2017-impl.hh parametric-dkmr2017.cc \
	parametric-dm2016.hh parametric-dm2016-impl.hh parametric-dm2016.cc \
	parametric-fvdv2018.hh parametric-fvdv2018-impl.hh parametric-fvdv2018.cc \
	parametric-kkrvd2024.hh parametric-kkrvd2024.cc \
	parametric-kkvdz2022.hh parametric-kkvdz2022.cc \
	parametric-kmpw2010.hh parametric-kmpw2010-impl.hh parametric-kmpw2010.cc \
	pi-lcdas.cc pi-lcdas.hh \
	psd-lcdas.cc psd-lcdas.hh \
	rho-lcdas.cc rho-lcdas.hh \
	unitarity-bounds.cc unitarity-bounds.hh \
	vec-lcdas.cc vec-lcdas.hh \
	zero-recoil-sum-rule.cc zero-recoil-sum-rule.hh
libeosformfactors_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
libeosformfactors_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)
libeosformfactors_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	-lgslcblas \
	-lgsl

include_eos_rarebdecaysdir = $(includedir)/eos/form-factors
include_eos_rarebdecays_HEADERS = \
	analytic-b-to-pi-pi.hh \
	analytic-b-to-psd-dkmmo2008.hh \
	form-factor-adapter.hh \
	form-factors.hh \
	heavy-meson-lcdas.hh \
	heavy-meson-lcdas-exponential.hh \
	heavy-meson-lcdas-flvd2022.hh \
	observables.hh \
	parametric-abr2022.hh \
	parametric-bbgorvd2018.hh \
	parametric-bcl2008.hh \
	parametric-bfw2010.hh \
	parametric-bfvd2014.hh \
	parametric-bgl1997.hh \
	parametric-bgjvd2019.hh \
	parametric-bsz2015.hh \
	parametric-dkmr2017.hh \
	parametric-dm2016.hh \
	parametric-fvdv2018.hh \
	parametric-kkrvd2024.hh \
	parametric-kmpw2010.hh \
	pi-lcdas.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	analytic-b-to-gamma-qcdf_TEST \
	analytic-b-to-pi_TEST \
	analytic-bs-to-k_TEST \
	analytic-b-to-pi-pi_TEST \
	analytic-b-to-p-lcsr_TEST \
	analytic-b-to-v-lcsr_TEST \
	baryonic_TEST \
	heavy-meson-lcdas_TEST \
	heavy-meson-lcdas-exponential_TEST \
	heavy-meson-lcdas-flvd2022_TEST \
	hqet-b-to-c_TEST \
	k-lcdas_TEST \
	kstar-lcdas_TEST \
	lcda_TEST \
	mesonic_TEST \
	parametric-abr2022_TEST \
	parametric-bbgorvd2018_TEST \
	parametric-bcl2008_TEST \
	parametric-bfvd2014_TEST \
	parametric-bfw2010_TEST \
	parametric-bgl1997_TEST \
	parametric-bgjvd2019_TEST \
	parametric-bmrvd2022_TEST \
	parametric-bsz2015_TEST \
	parametric-dkmr2017_TEST \
	parametric-dm2016_TEST \
	parametric-fvdv2018_TEST \
	parametric-kkrvd2024_TEST \
	parametric-kkvdz2022_TEST \
	parametric-kmpw2010_TEST \
	pi-lcdas_TEST \
	psd-lcdas_TEST \
	rho-lcdas_TEST \
	unitarity-bounds_TEST \
	vec-lcdas_TEST
LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosformfactors.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

analytic_b_to_gamma_qcdf_TEST_SOURCES = analytic-b-to-gamma-qcdf_TEST.cc

analytic_b_to_pi_TEST_SOURCES = analytic-b-to-pi_TEST.cc

analytic_bs_to_k_TEST_SOURCES = analytic-bs-to-k_TEST.cc

analytic_b_to_pi_pi_TEST_SOURCES = analytic-b-to-pi-pi_TEST.cc
analytic_b_to_pi_pi_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)

analytic_b_to_p_lcsr_TEST_SOURCES = analytic-b-to-p-lcsr_TEST.cc

analytic_b_to_v_lcsr_TEST_SOURCES = analytic-b-to-v-lcsr_TEST.cc

baryonic_TEST_SOURCES = baryonic_TEST.cc

heavy_meson_lcdas_TEST_SOURCES = heavy-meson-lcdas_TEST.cc

heavy_meson_lcdas_exponential_TEST_SOURCES = heavy-meson-lcdas-exponential_TEST.cc

heavy_meson_lcdas_flvd2022_TEST_SOURCES = heavy-meson-lcdas-flvd2022_TEST.cc

hqet_b_to_c_TEST_SOURCES = hqet-b-to-c_TEST.cc

k_lcdas_TEST_SOURCES = k-lcdas_TEST.cc

kstar_lcdas_TEST_SOURCES = kstar-lcdas_TEST.cc

lcda_TEST_SOURCES = lcda_TEST.cc

mesonic_TEST_SOURCES = mesonic_TEST.cc

parametric_abr2022_TEST_SOURCES = parametric-abr2022_TEST.cc

parametric_bbgorvd2018_TEST_SOURCES = parametric-bbgorvd2018_TEST.cc

parametric_bcl2008_TEST_SOURCES = parametric-bcl2008_TEST.cc

parametric_bfvd2014_TEST_SOURCES = parametric-bfvd2014_TEST.cc

parametric_bfw2010_TEST_SOURCES = parametric-bfw2010_TEST.cc

parametric_bgl1997_TEST_SOURCES = parametric-bgl1997_TEST.cc

parametric_bgjvd2019_TEST_SOURCES = parametric-bgjvd2019_TEST.cc

parametric_bsz2015_TEST_SOURCES = parametric-bsz2015_TEST.cc

parametric_bmrvd2022_TEST_SOURCES = parametric-bmrvd2022_TEST.cc

parametric_dkmr2017_TEST_SOURCES = parametric-dkmr2017_TEST.cc

parametric_dm2016_TEST_SOURCES = parametric-dm2016_TEST.cc

parametric_fvdv2018_TEST_SOURCES = parametric-fvdv2018_TEST.cc

parametric_kkrvd2024_TEST_SOURCES = parametric-kkrvd2024_TEST.cc

parametric_kkvdz2022_TEST_SOURCES = parametric-kkvdz2022_TEST.cc

parametric_kmpw2010_TEST_SOURCES = parametric-kmpw2010_TEST.cc

pi_lcdas_TEST_SOURCES = pi-lcdas_TEST.cc

psd_lcdas_TEST_SOURCES = psd-lcdas_TEST.cc

rho_lcdas_TEST_SOURCES = rho-lcdas_TEST.cc

unitarity_bounds_TEST_SOURCES = unitarity-bounds_TEST.cc

vec_lcdas_TEST_SOURCES = vec-lcdas_TEST.cc
