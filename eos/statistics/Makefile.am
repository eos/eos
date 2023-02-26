CLEANFILES = \
	*~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosstatistics.la
libeosstatistics_la_SOURCES = \
	goodness-of-fit.cc goodness-of-fit.hh \
	log-likelihood.cc log-likelihood.hh log-likelihood-fwd.hh \
	log-posterior.cc log-posterior.hh log-posterior-fwd.hh \
	log-prior.cc log-prior.hh log-prior-fwd.hh \
	test-statistic.cc test-statistic.hh test-statistic-impl.hh
libeosstatistics_la_LIBADD = -lpthread -lgsl -lgslcblas -lm -lyaml-cpp
libeosstatistics_la_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS) $(YAMLCPP_CXXFLAGS)
libeosstatistics_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS) $(YAMLCPP_LDFLAGS)

include_eos_statisticsdir = $(includedir)/eos/statistics
include_eos_statistics_HEADERS = \
	goodness-of-fit.hh \
	log-likelihood.hh log-likelihood-fwd.hh \
	log-posterior.hh log-posterior-fwd.hh \
	log-prior.hh log-prior-fwd.hh \
	test-statistic.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	log-likelihood_TEST \
	log-posterior_TEST \
	log-prior_TEST
LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosstatistics.la \
	$(top_builddir)/eos/libeos.la \
	$(top_builddir)/eos/utils/libeosutils.la

check_PROGRAMS = $(TESTS)

log_likelihood_TEST_SOURCES = log-likelihood_TEST.cc
log_likelihood_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
log_likelihood_TEST_LDFLAGS = $(GSL_LDFLAGS)

log_posterior_TEST_SOURCES = log-posterior_TEST.cc log-posterior_TEST.hh
log_posterior_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)

log_prior_TEST_SOURCES = log-prior_TEST.cc
log_prior_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
log_prior_TEST_LDFLAGS = $(GSL_LDFLAGS)
