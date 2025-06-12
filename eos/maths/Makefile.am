CLEANFILES = \
	*~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

lib_LTLIBRARIES = libeosmaths.la
libeosmaths_la_SOURCES = \
	angular-integrals.cc angular-integrals.hh \
	complex.hh \
	derivative.cc derivative.hh \
	gegenbauer-polynomial.cc gegenbauer-polynomial.hh \
	gsl-interface.hh \
	integrate.cc integrate.hh integrate-impl.hh \
	integrate-cubature.cc integrate-cubature.hh \
	interpolation.cc interpolation.hh \
	lagrange-polynomial.hh \
	legendre-polynomial-vector.hh \
	matrix.hh \
	multiplepolylog-li22.cc multiplepolylog-li22.hh \
	multiplepolylog-li22-const.cc multiplepolylog-li22-const.hh \
	omnes-factor.hh omnes-factor-impl.hh \
	outer-function.hh outer-function.cc \
	polylog.cc polylog.hh \
	power-of.hh \
	szego-polynomial.hh

libeosmaths_la_LIBADD = \
	-lgsl -lgslcblas -lm
libeosmaths_la_CXXFLAGS = $(AM_CXXFLAGS) \
	-DEOS_DATADIR='"$(datadir)"' \
	$(GSL_CXXFLAGS)
libeosmaths_la_LDFLAGS = $(AM_LDFLAGS) $(GSL_LDFLAGS)

EXTRA_DIST = \
	polylog_TEST_dilog.bin \
	polylog_TEST_quadlog.bin \
	polylog_TEST_trilog.bin \
	polylog_TEST_z.bin

include_eos_utilsdir = $(includedir)/eos/utils
include_eos_utils_HEADERS = \
    angular-integrals.hh \
	complex.hh \
	derivative.hh \
	gegenbauer-polynomial.hh \
	gsl-interface.hh \
	integrate.hh \
	integrate-cubature.hh \
	interpolation.hh \
	lagrange-polynomial.hh \
	legendre-polynomial-vector.hh \
	matrix.hh \
	multiplepolylog-li22.hh \
	multiplepolylog-li22-const.hh \
	omnes-factor.hh omnes-factor-impl.hh \
	outer-function.hh \
	polylog.hh \
	power-of.hh \
	szego-polynomial.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
    angular-integrals_TEST \
	derivative_TEST \
	gegenbauer-polynomial_TEST \
	gsl-interface_TEST \
	integrate_TEST \
	interpolation_TEST \
	lagrange-polynomial_TEST \
	legendre-polynomial-vector_TEST \
	matrix_TEST \
	multiplepolylog_TEST \
	omnes-factor_TEST \
	outer-function_TEST \
	polylog_TEST \
	power_of_TEST \
	szego-polynomial_TEST
LDADD = \
	$(top_builddir)/test/libeostest.la \
	libeosmaths.la \
	$(top_builddir)/eos/maths/libeosmaths.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

angular_integrals_TEST_SOURCES = angular-integrals_TEST.cc

derivative_TEST_SOURCES = derivative_TEST.cc

gegenbauer_polynomial_TEST_SOURCES = gegenbauer-polynomial_TEST.cc

gsl_interface_TEST_SOURCES = gsl-interface_TEST.cc

integrate_TEST_SOURCES = integrate_TEST.cc
integrate_TEST_CXXFLAGS = $(AM_CXXFLAGS) $(GSL_CXXFLAGS)
integrate_TEST_LDFLAGS = $(GSL_LDFLAGS)

interpolation_TEST_SOURCES = interpolation_TEST.cc

lagrange_polynomial_TEST_SOURCES = lagrange-polynomial_TEST.cc

legendre_polynomial_vector_TEST_SOURCES = legendre-polynomial-vector_TEST.cc

matrix_TEST_SOURCES = matrix_TEST.cc

multiplepolylog_TEST_SOURCES = multiplepolylog_TEST.cc

omnes_factor_TEST_SOURCES = omnes-factor_TEST.cc

outer_function_TEST_SOURCES = outer-function_TEST.cc

polylog_TEST_SOURCES = polylog_TEST.cc

power_of_TEST_SOURCES = power-of_TEST.cc

szego_polynomial_TEST_SOURCES = szego-polynomial_TEST.cc
