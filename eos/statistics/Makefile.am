CLEANFILES = \
	*~ \
	markov_chain_sampler_TEST.hdf5 \
	prior_sampler_TEST.hdf5 \
	proposal_functions_TEST-rdwr.hdf5 \
	proposal_functions_TEST-block-decomposition.hdf5
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosstatistics.la
libeosstatistics_la_SOURCES = \
	analysis.cc analysis.hh analysis-fwd.hh \
	chi-squared.hh chi-squared.cc \
	cluster.cc cluster.hh \
	hierarchical-clustering.cc hierarchical-clustering.hh \
	histogram.cc histogram.hh \
	log_likelihood.cc log_likelihood.hh log_likelihood-fwd.hh \
	log_prior.cc log_prior.hh log_prior-fwd.hh \
	markov_chain.cc markov_chain.hh \
	markov_chain_sampler.cc markov_chain_sampler.hh \
	prior_sampler.cc prior_sampler.hh \
	proposal_functions.cc proposal_functions.hh \
	rvalue.cc rvalue.hh \
	welford.cc welford.hh
libeosstatistics_la_LIBADD = -lpthread -lgsl -lgslcblas -lm -lMinuit2

if EOS_ENABLE_PMC
libeosstatistics_la_SOURCES += \
	population_monte_carlo_sampler.cc population_monte_carlo_sampler.hh
libeosstatistics_la_LIBADD += -lpmc -ldl
endif

include_eos_statisticsdir = $(includedir)/eos/utils
include_eos_statistics_HEADERS = \
	analysis.hh analysis-fwd.hh \
	chi-squared.hh \
	hierarchical-clustering.hh \
	histogram.hh \
	log_likelihood.hh log_likelihood-fwd.hh \
	log_prior.hh log_prior-fwd.hh \
	markov_chain.hh \
	markov_chain_sampler.hh \
	prior_sampler.hh \
	proposal_functions.hh \
	rvalue.hh \
	welford.hh

TESTS = \
	analysis_TEST \
	chi-squared_TEST \
	hierarchical-clustering_TEST \
	histogram_TEST \
	log_likelihood_TEST \
	log_prior_TEST \
	markov_chain_TEST \
	markov_chain_sampler_TEST \
	prior_sampler_TEST \
	proposal_functions_TEST \
	rvalue_TEST \
	welford_TEST
LDADD = \
	-lgsl \
	-lhdf5 \
	-lMinuit2 \
	$(top_builddir)/test/libeostest.a \
	libeosstatistics.la \
	$(top_builddir)/eos/libeos.la \
	$(top_builddir)/eos/utils/libeosutils.la

if EOS_ENABLE_PMC
TESTS += population_monte_carlo_sampler_TEST
LDADD += -lpmc -ldl
endif

check_PROGRAMS = $(TESTS)

analysis_TEST_SOURCES = analysis_TEST.cc analysis_TEST.hh

chi_squared_TEST_SOURCES = chi-squared_TEST.cc

hierarchical_clustering_TEST_SOURCES = hierarchical-clustering_TEST.cc

histogram_TEST_SOURCES = histogram_TEST.cc

log_likelihood_TEST_SOURCES = log_likelihood_TEST.cc

log_prior_TEST_SOURCES = log_prior_TEST.cc

markov_chain_TEST_SOURCES = markov_chain_TEST.cc

markov_chain_sampler_TEST_SOURCES = markov_chain_sampler_TEST.cc

if EOS_ENABLE_PMC
population_monte_carlo_sampler_TEST_SOURCES = population_monte_carlo_sampler_TEST.cc
endif

prior_sampler_TEST_SOURCES = prior_sampler_TEST.cc

proposal_functions_TEST_SOURCES = proposal_functions_TEST.cc

rvalue_TEST_SOURCES = rvalue_TEST.cc

welford_TEST_SOURCES = welford_TEST.cc
