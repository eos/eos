CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

lib_LTLIBRARIES = libwfutils.la
libwfutils_la_SOURCES = \
	complex.cc complex.hh complex-impl.hh \
	exception.cc exception.hh \
	kinematics.hh \
	parameters.cc parameters.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh

TESTS = \
	complex_TEST
TESTS_LD_ADD = $(top_builddir)/test/libwftest.a libwfutils.la

check_PROGRAMS = $(TESTS)
complex_TEST_SOURCES = complex_TEST.cc
complex_TEST_LDADD = $(TESTS_LD_ADD)
