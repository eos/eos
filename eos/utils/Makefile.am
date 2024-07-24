CLEANFILES = \
	*~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosutils.la
libeosutils_la_SOURCES = \
	cartesian-product.hh \
	concrete_observable.cc concrete_observable.hh \
	concrete-cacheable-observable.hh \
	concrete-signal-pdf.hh \
	condition_variable.cc condition_variable.hh \
	density.cc density.hh density-fwd.hh density-impl.hh \
	destringify.cc destringify.hh \
	diagnostics.cc diagnostics.hh \
	exception.cc exception.hh \
	expression.cc expression.hh expression-fwd.hh \
	expression-cacher.hh \
	expression-cloner.hh \
	expression-evaluator.hh \
	expression-kinematic-reader.hh \
	expression-maker.hh \
	expression-observable.cc expression-observable.hh \
	expression-parser.cc expression-parser.hh expression-parser-impl.hh \
	expression-printer.hh \
	expression-used-parameter-reader.hh \
	expression-visitors.cc \
	gsl-hacks.cc \
	indirect-iterator.hh indirect-iterator-fwd.hh indirect-iterator-impl.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	iterator-range.hh \
	join.hh \
	kinematic.cc kinematic.hh \
	kmatrix.cc kmatrix.hh kmatrix-impl.hh\
	lock.cc lock.hh \
	log.cc log.hh \
	memoise.cc memoise.hh \
	mutable.cc mutable.hh mutable-fwd.hh \
	mutex.cc mutex.hh \
	named-value.hh named-value-fwd.hh \
	observable_cache.cc observable_cache.hh \
	observable_set.cc observable_set.hh \
	observable_stub.cc observable_stub.hh \
	omnes-factor.hh omnes-factor-impl.hh \
	one-of.hh \
	options.cc options.hh options-impl.hh \
	parameters.cc parameters.hh parameters-fwd.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	qualified-name.cc qualified-name.hh \
	quantum-numbers.cc quantum-numbers.hh \
	reference-name.cc reference-name.hh \
	stringify.hh \
	test-observable.cc test-observable.hh \
	thread.cc thread.hh \
	thread_pool.cc thread_pool.hh \
	ticket.cc ticket.hh \
	transitions.hh \
	tuple-maker.hh \
	type-list.hh type-list-fwd.hh \
	units.cc units.hh \
	verify.cc verify.hh \
	visitor.hh visitor-fwd.hh \
	wilson-polynomial.cc wilson-polynomial.hh \
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh

libeosutils_la_LIBADD = \
	-lboost_filesystem -lboost_system \
	-lgsl -lgslcblas -lm \
	-lpthread \
	-lyaml-cpp
libeosutils_la_CXXFLAGS = $(AM_CXXFLAGS) \
	-DEOS_DATADIR='"$(datadir)"' \
	$(GSL_CXXFLAGS) \
	$(YAMLCPP_CXXFLAGS)
libeosutils_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS) $(YAMLCPP_LDFLAGS)

include_eos_utilsdir = $(includedir)/eos/utils
include_eos_utils_HEADERS = \
	cartesian-product.hh \
	concrete_observable.hh \
	concrete-signal-pdf.hh \
	condition_variable.hh \
	density.hh density-fwd.hh \
	destringify.hh \
	exception.hh \
	expression.hh expression-fwd.hh \
	expression-parser.hh expression-parser-impl.hh \
	indirect-iterator.hh indirect-iterator-fwd.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	iterator-range.hh \
	join.hh \
	kinematic.hh \
	kmatrix.hh kmatrix-impl.hh\
	lock.hh \
	log.hh \
	memoise.hh \
	mutex.hh \
	mutable.hh mutable-fwd.hh \
	observable_cache.hh \
	observable_set.hh \
	omnes-factor.hh omnes-factor-impl.hh \
	one-of.hh \
	options.hh \
	parameters.hh parameters-fwd.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.hh \
	qualified-name.hh \
	quantum-numbers.hh \
	reference-name.hh \
	rge.hh rge-impl.hh \
	stringify.hh \
	thread.hh \
	thread_pool.hh \
	ticket.hh \
	transitions.hh \
	tuple-maker.hh \
	units.hh \
	verify.hh \
	wilson-polynomial.hh \
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	cacheable-observable_TEST \
	cartesian-product_TEST \
	expression-parser_TEST \
	gsl-hacks_TEST \
	indirect-iterator_TEST \
	join_TEST \
	kinematic_TEST \
	kmatrix_TEST \
	log_TEST \
	memoise_TEST \
	mutable_TEST \
	observable_set_TEST \
	observable_stub_TEST \
	options_TEST \
	omnes-factor_TEST \
	one-of_TEST \
	parameters_TEST \
	qcd_TEST \
	qualified-name_TEST \
	quantum-numbers_TEST \
	reference-name_TEST \
	rge_TEST \
	stringify_TEST \
	verify_TEST \
	wilson-polynomial_TEST
LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosutils.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/models/libeosmodels.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

cacheable_observable_TEST_SOURCES = cacheable_observable_TEST.cc

cartesian_product_TEST_SOURCES = cartesian-product_TEST.cc

expression_parser_TEST_SOURCES = expression-parser_TEST.cc

gsl_hacks_TEST_SOURCES = gsl-hacks_TEST.cc
gsl_hacks_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
gsl_hacks_TEST_LDFLAGS = $(GSL_LDFLAGS)

indirect_iterator_TEST_SOURCES = indirect-iterator_TEST.cc

join_TEST_SOURCES = join_TEST.cc

kinematic_TEST_SOURCES = kinematic_TEST.cc

kmatrix_TEST_SOURCES = kmatrix_TEST.cc

log_TEST_SOURCES = log_TEST.cc

memoise_TEST_SOURCES = memoise_TEST.cc

mutable_TEST_SOURCES = mutable_TEST.cc

observable_set_TEST_SOURCES = observable_set_TEST.cc

observable_stub_TEST_SOURCES = observable_stub_TEST.cc

omnes_factor_TEST_SOURCES = omnes-factor_TEST.cc

one_of_TEST_SOURCES = one-of_TEST.cc

options_TEST_SOURCES = options_TEST.cc

parameters_TEST_SOURCES = parameters_TEST.cc

qcd_TEST_SOURCES = qcd_TEST.cc

qualified_name_TEST_SOURCES = qualified-name_TEST.cc

quantum_numbers_TEST_SOURCES = quantum-numbers_TEST.cc

reference_name_TEST_SOURCES = reference-name_TEST.cc

rge_TEST_SOURCES = rge_TEST.cc

stringify_TEST_SOURCES = stringify_TEST.cc

verify_TEST_SOURCES = verify_TEST.cc

wilson_polynomial_TEST_SOURCES = wilson-polynomial_TEST.cc
