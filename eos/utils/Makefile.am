CLEANFILES = \
	*~ \
	markov_chain_sampler_TEST.hdf5 \
	scan_file_TEST.hdf5
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosutils.la
libeosutils_la_SOURCES = \
	accumulator.cc accumulator.hh \
	analysis.cc analysis.hh analysis-fwd.hh \
	apply.hh \
	cartesian-product.hh \
	chi-squared.hh chi-squared.cc \
	complex.hh \
	concrete_observable.hh \
	condition_variable.cc condition_variable.hh \
	destringify.hh \
	exception.cc exception.hh \
	histogram.cc histogram.hh \
	integrate.cc integrate.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	kinematic.cc kinematic.hh \
	lcda.cc lcda.hh \
	lock.cc lock.hh \
	log.cc log.hh \
	log_likelihood.cc log_likelihood.hh log_likelihood-fwd.hh \
	log_prior.cc log_prior.hh log_prior-fwd.hh \
	markov_chain.cc markov_chain.hh \
	markov_chain_sampler.cc markov_chain_sampler.hh \
	matrix.hh \
	memoise.hh \
	model.cc model.hh \
	mutex.cc mutex.hh \
	observable_cache.cc observable_cache.hh \
	observable_set.cc observable_set.hh \
	one-of.hh \
	options.cc options.hh \
	parameters.cc parameters.hh \
	power_of.hh \
	prior_sampler.cc prior_sampler.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	random_number_generator.cc random_number_generator.hh \
	rvalue.cc rvalue.hh \
	save.hh \
	scan_file.cc scan_file.hh \
	standard-model.cc standard-model.hh \
	stringify.hh \
	thread.cc thread.hh \
	thread_pool.cc thread_pool.hh \
	ticket.cc ticket.hh \
	top-loops.hh top-loops.cc \
	verify.cc verify.hh \
	wilson_coefficients.cc wilson_coefficients.hh \
	wilson-polynomial.cc wilson-polynomial.hh \
	wilson_scan_model.cc wilson_scan_model.hh \
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh
libeosutils_la_LIBADD = -lpthread -lgslcblas -lgsl -lhdf5 -lhdf5_hl

TESTS = \
	analysis_TEST \
	apply_TEST \
	cartesian-product_TEST \
	chi-squared_TEST \
	histogram_TEST \
	integrate_TEST \
	kinematic_TEST \
	lcda_TEST \
	log_TEST \
	log_likelihood_TEST \
	log_prior_TEST \
	matrix_TEST \
	markov_chain_TEST \
	markov_chain_sampler_TEST \
	memoise_TEST \
	observable_set_TEST \
	options_TEST \
	one-of_TEST \
	parameters_TEST \
	power_of_TEST \
	prior_sampler_TEST \
	qcd_TEST \
	random_number_generator_TEST \
	rvalue_TEST \
	save_TEST \
	scan_file_TEST \
	standard_model_TEST \
	top-loops_TEST \
	stringify_TEST \
	verify_TEST \
	wilson_coefficients_TEST \
	wilson-polynomial_TEST \
	wilson_scan_model_TEST
TESTS_LD_ADD = $(top_builddir)/test/libeostest.a $(top_builddir)/eos/libeos.la libeosutils.la

check_PROGRAMS = $(TESTS)

analysis_TEST_SOURCES = analysis_TEST.cc analysis_TEST.hh
analysis_TEST_LDADD = $(TESTS_LD_ADD)

apply_TEST_SOURCES = apply_TEST.cc
apply_TEST_LDADD = $(TESTS_LD_ADD)

cartesian_product_TEST_SOURCES = cartesian-product_TEST.cc
cartesian_product_TEST_LDADD = $(TESTS_LD_ADD)

chi_squared_TEST_SOURCES = chi-squared_TEST.cc
chi_squared_TEST_LDADD = $(TESTS_LD_ADD)

histogram_TEST_SOURCES = histogram_TEST.cc
histogram_TEST_LDADD = $(TESTS_LD_ADD)

integrate_TEST_SOURCES = integrate_TEST.cc
integrate_TEST_LDADD = $(TESTS_LD_ADD)

kinematic_TEST_SOURCES = kinematic_TEST.cc
kinematic_TEST_LDADD = $(TESTS_LD_ADD)

lcda_TEST_SOURCES = lcda_TEST.cc
lcda_TEST_LDADD = $(TESTS_LD_ADD)

log_TEST_SOURCES = log_TEST.cc
log_TEST_LDADD = $(TESTS_LD_ADD)

log_likelihood_TEST_SOURCES = log_likelihood_TEST.cc
log_likelihood_TEST_LDADD = $(TESTS_LD_ADD)

log_prior_TEST_SOURCES = log_prior_TEST.cc
log_prior_TEST_LDADD = $(TESTS_LD_ADD)

matrix_TEST_SOURCES = matrix_TEST.cc
matrix_TEST_LDADD = $(TESTS_LD_ADD)

markov_chain_TEST_SOURCES = markov_chain_TEST.cc
markov_chain_TEST_LDADD = $(TESTS_LD_ADD)

markov_chain_sampler_TEST_SOURCES = markov_chain_sampler_TEST.cc
markov_chain_sampler_TEST_LDADD = $(TESTS_LD_ADD) $(top_builddir)/eos/libeos.la

memoise_TEST_SOURCES = memoise_TEST.cc
memoise_TEST_LDADD = $(TESTS_LD_ADD)

observable_set_TEST_SOURCES = observable_set_TEST.cc
observable_set_TEST_LDADD = $(TESTS_LD_ADD) 

one_of_TEST_SOURCES = one-of_TEST.cc
one_of_TEST_LDADD = $(TESTS_LD_ADD)

options_TEST_SOURCES = options_TEST.cc
options_TEST_LDADD = $(TESTS_LD_ADD)

parameters_TEST_SOURCES = parameters_TEST.cc
parameters_TEST_LDADD = $(TESTS_LD_ADD)

power_of_TEST_SOURCES = power_of_TEST.cc
power_of_TEST_LDADD = $(TESTS_LD_ADD)

prior_sampler_TEST_SOURCES = prior_sampler_TEST.cc
prior_sampler_TEST_LDADD = $(TESTS_LD_ADD)

qcd_TEST_SOURCES = qcd_TEST.cc
qcd_TEST_LDADD = $(TESTS_LD_ADD)

random_number_generator_TEST_SOURCES = random_number_generator_TEST.cc
random_number_generator_TEST_LDADD = $(TESTS_LD_ADD)

rvalue_TEST_SOURCES = rvalue_TEST.cc
rvalue_TEST_LDADD = $(TESTS_LD_ADD)

save_TEST_SOURCES = save_TEST.cc
save_TEST_LDADD = $(TESTS_LD_ADD)

scan_file_TEST_SOURCES = scan_file_TEST.cc
scan_file_TEST_LDADD = $(TESTS_LD_ADD)

stringify_TEST_SOURCES = stringify_TEST.cc
stringify_TEST_LDADD = $(TESTS_LD_ADD)

standard_model_TEST_SOURCES = standard_model_TEST.cc
standard_model_TEST_LDADD = $(TESTS_LD_ADD)

top_loops_TEST_SOURCES = top-loops_TEST.cc
top_loops_TEST_LDADD = $(TESTS_LD_ADD)

verify_TEST_SOURCES = verify_TEST.cc
verify_TEST_LDADD = $(TESTS_LD_ADD)

wilson_coefficients_TEST_SOURCES = wilson_coefficients_TEST.cc
wilson_coefficients_TEST_LDADD = $(TESTS_LD_ADD)

wilson_polynomial_TEST_SOURCES = wilson-polynomial_TEST.cc
wilson_polynomial_TEST_LDADD = $(TESTS_LD_ADD)

wilson_scan_model_TEST_SOURCES = wilson_scan_model_TEST.cc
wilson_scan_model_TEST_LDADD = $(TESTS_LD_ADD)
