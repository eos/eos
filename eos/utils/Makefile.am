CLEANFILES = \
	*~ \
	hdf5_TEST-attribute.hdf5 \
	hdf5_TEST-file.hdf5 \
	hdf5_TEST-copy.hdf5
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeosutils.la
libeosutils_la_SOURCES = \
	accumulator.cc accumulator.hh \
	apply.hh \
	cartesian-product.hh \
	ckm_scan_model.cc ckm_scan_model.hh \
	complex.hh \
	concrete_observable.hh \
	concrete-signal-pdf.hh \
	condition_variable.cc condition_variable.hh \
	density.cc density.hh density-fwd.hh density-impl.hh \
	derivative.cc derivative.hh \
	destringify.cc destringify.hh \
	diagnostics.cc diagnostics.hh \
	equation_solver.cc equation_solver.hh \
	exception.cc exception.hh \
	gsl-cblas-hack.cc \
	hdf5.cc hdf5.hh hdf5-fwd.hh \
	indirect-iterator.hh indirect-iterator-fwd.hh indirect-iterator-impl.hh \
	integrate.cc integrate.hh integrate-impl.hh \
	integrate-cubature.hh integrate-cubature.cc \
	instantiation_policy.hh instantiation_policy-impl.hh \
	iterator-range.hh \
	join.hh \
	kinematic.cc kinematic.hh \
	lcda.cc lcda.hh \
	lock.cc lock.hh \
	log.cc log.hh \
	matrix.hh \
	memoise.cc memoise.hh \
	model.cc model.hh \
	mutable.cc mutable.hh mutable-fwd.hh \
	mutex.cc mutex.hh \
	named-value.hh named-value-fwd.hh \
	observable_cache.cc observable_cache.hh \
	observable_set.cc observable_set.hh \
	observable_stub.cc observable_stub.hh \
	one-of.hh \
	options.cc options.hh options-impl.hh \
	parameters.cc parameters.hh parameters-fwd.hh \
	polylog.cc polylog.hh \
	power_of.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	qualified-name.cc qualified-name.hh \
	random_number_generator.cc random_number_generator.hh \
	save.hh \
	standard-model.cc standard-model.hh \
	stringify.hh \
	thread.cc thread.hh \
	thread_pool.cc thread_pool.hh \
	ticket.cc ticket.hh \
	top-loops.hh top-loops.cc \
	tuple-maker.hh \
	type-list.hh type-list-fwd.hh \
	verify.cc verify.hh \
	visitor.hh visitor-fwd.hh \
	wilson_coefficients.cc wilson_coefficients.hh \
	wilson-polynomial.cc wilson-polynomial.hh \
	wilson_scan_model.cc wilson_scan_model.hh \
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh

libeosutils_la_LIBADD = \
	-lboost_filesystem -lboost_system \
	-lgsl -lgslcblas -lm \
	-lhdf5 -lhdf5_hl \
	-lMinuit2 \
	-lpthread \
	-lyaml-cpp
libeosutils_la_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS) $(MINUIT2_CXXFLAGS)
libeosutils_la_LDFLAGS = $(AM_CXXFLAGS) $(HDF5_LDFLAGS) $(MINUIT2_LDFLAGS)

EXTRA_DIST = \
	polylog_TEST_dilog.bin \
	polylog_TEST_trilog.bin \
	polylog_TEST_z.bin

include_eos_utilsdir = $(includedir)/eos/utils
include_eos_utils_HEADERS = \
	accumulator.hh \
	apply.hh \
	cartesian-product.hh \
	ckm_scan_model.hh \
	complex.hh \
	concrete_observable.hh \
	concrete-signal-pdf.hh \
	condition_variable.hh \
	density.hh density-fwd.hh \
	derivative.hh \
	destringify.hh \
	equation_solver.hh \
	exception.hh \
	hdf5.hh hdf5-fwd.hh \
	indirect-iterator.hh indirect-iterator-fwd.hh \
	integrate.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	iterator-range.hh \
	join.hh \
	kinematic.hh \
	lcda.hh \
	lock.hh \
	log.hh \
	matrix.hh \
	memoise.hh \
	model.hh \
	mutex.hh \
	mutable.hh mutable-fwd.hh \
	observable_cache.hh \
	observable_set.hh \
	one-of.hh \
	options.hh \
	parameters.hh parameters-fwd.hh \
	power_of.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.hh \
	qualified-name.hh \
	random_number_generator.cc random_number_generator.hh \
	save.hh \
	standard-model.hh \
	stringify.hh \
	thread.hh \
	thread_pool.hh \
	ticket.hh \
	top-loops.hh \
	tuple-maker.hh \
	verify.hh \
	wilson_coefficients.hh \
	wilson-polynomial.hh \
	wilson_scan_model.hh \
	wrapped_forward_iterator.hh wrapped_forward_iterator-fwd.hh wrapped_forward_iterator-impl.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	apply_TEST \
	cartesian-product_TEST \
	ckm_scan_model_TEST \
	derivative_TEST \
	equation_solver_TEST \
	hdf5_TEST \
	indirect-iterator_TEST \
	integrate_TEST \
	join_TEST \
	kinematic_TEST \
	lcda_TEST \
	log_TEST \
	matrix_TEST \
	memoise_TEST \
	mutable_TEST \
	observable_set_TEST \
	observable_stub_TEST \
	options_TEST \
	one-of_TEST \
	parameters_TEST \
	polylog_TEST \
	power_of_TEST \
	qcd_TEST \
	qualified-name_TEST \
	random_number_generator_TEST \
	save_TEST \
	standard_model_TEST \
	top-loops_TEST \
	stringify_TEST \
	verify_TEST \
	wilson_coefficients_TEST \
	wilson-polynomial_TEST \
	wilson_scan_model_TEST
LDADD = \
	-lyaml-cpp \
	$(top_builddir)/test/libeostest.a \
	libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

apply_TEST_SOURCES = apply_TEST.cc

cartesian_product_TEST_SOURCES = cartesian-product_TEST.cc

ckm_scan_model_TEST_SOURCES = ckm_scan_model_TEST.cc

derivative_TEST_SOURCES = derivative_TEST.cc

hdf5_TEST_SOURCES = hdf5_TEST.cc
hdf5_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS)
hdf5_TEST_LDFLAGS = $(AM_CXXFLAGS) $(HDF5_LDFLAGS)

equation_solver_TEST_SOURCES = equation_solver_TEST.cc
equation_solver_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(MINUIT2_CXXFLAGS)
equation_solver_TEST_LDFLAGS = $(MINUIT2_LDFLAGS)

indirect_iterator_TEST_SOURCES = indirect-iterator_TEST.cc

integrate_TEST_SOURCES = integrate_TEST.cc

join_TEST_SOURCES = join_TEST.cc

kinematic_TEST_SOURCES = kinematic_TEST.cc

lcda_TEST_SOURCES = lcda_TEST.cc

log_TEST_SOURCES = log_TEST.cc

matrix_TEST_SOURCES = matrix_TEST.cc

memoise_TEST_SOURCES = memoise_TEST.cc

mutable_TEST_SOURCES = mutable_TEST.cc

observable_set_TEST_SOURCES = observable_set_TEST.cc

observable_stub_TEST_SOURCES = observable_stub_TEST.cc

one_of_TEST_SOURCES = one-of_TEST.cc

options_TEST_SOURCES = options_TEST.cc

parameters_TEST_SOURCES = parameters_TEST.cc

polylog_TEST_SOURCES = polylog_TEST.cc

power_of_TEST_SOURCES = power_of_TEST.cc

qcd_TEST_SOURCES = qcd_TEST.cc

qualified_name_TEST_SOURCES = qualified-name_TEST.cc

random_number_generator_TEST_SOURCES = random_number_generator_TEST.cc

save_TEST_SOURCES = save_TEST.cc

stringify_TEST_SOURCES = stringify_TEST.cc

standard_model_TEST_SOURCES = standard_model_TEST.cc

top_loops_TEST_SOURCES = top-loops_TEST.cc

verify_TEST_SOURCES = verify_TEST.cc

wilson_coefficients_TEST_SOURCES = wilson_coefficients_TEST.cc

wilson_polynomial_TEST_SOURCES = wilson-polynomial_TEST.cc

wilson_scan_model_TEST_SOURCES = wilson_scan_model_TEST.cc
