CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosutils.la
libeosutils_la_SOURCES = \
	accumulator.cc accumulator.hh \
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
	likelihood.cc likelihood.hh \
	lcda.cc lcda.hh \
	lock.cc lock.hh \
	matrix.hh \
	memoise.hh \
	model.cc model.hh \
	mutex.cc mutex.hh \
	one-of.hh \
	options.cc options.hh \
	parameters.cc parameters.hh \
	power_of.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	random_number_engine.cc random_number_engine.hh \
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
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh
libeosutils_la_LIBADD = -lpthread -lgslcblas -lgsl -lhdf5 -lhdf5_hl

TESTS = \
	apply_TEST \
	cartesian-product_TEST \
	chi-squared_TEST \
	histogram_TEST \
	integrate_TEST \
	lcda_TEST \
	likelihood_TEST \
	matrix_TEST \
	memoise_TEST \
	one-of_TEST \
	parameters_TEST \
	power_of_TEST \
	qcd_TEST \
	save_TEST \
	scan_file_TEST \
	standard_model_TEST \
	top-loops_TEST \
	stringify_TEST \
	verify_TEST \
	wilson_coefficients_TEST \
	wilson-polynomial_TEST
TESTS_LD_ADD = $(top_builddir)/test/libeostest.a $(top_builddir)/src/libeos.la libeosutils.la

check_PROGRAMS = $(TESTS)

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

lcda_TEST_SOURCES = lcda_TEST.cc
lcda_TEST_LDADD = $(TESTS_LD_ADD)

likelihood_TEST_SOURCES = likelihood_TEST.cc
likelihood_TEST_LDADD = $(TESTS_LD_ADD)

matrix_TEST_SOURCES = matrix_TEST.cc
matrix_TEST_LDADD = $(TESTS_LD_ADD)

memoise_TEST_SOURCES = memoise_TEST.cc
memoise_TEST_LDADD = $(TESTS_LD_ADD)

one_of_TEST_SOURCES = one-of_TEST.cc
one_of_TEST_LDADD = $(TESTS_LD_ADD)

parameters_TEST_SOURCES = parameters_TEST.cc
parameters_TEST_LDADD = $(TESTS_LD_ADD)

power_of_TEST_SOURCES = power_of_TEST.cc
power_of_TEST_LDADD = $(TESTS_LD_ADD)

qcd_TEST_SOURCES = qcd_TEST.cc
qcd_TEST_LDADD = $(TESTS_LD_ADD)

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
