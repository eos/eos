MAINTAINERCLEANFILES = Makefile.in

lib_LTLIBRARIES = libeosoptimize.la
libeosoptimize_la_SOURCES = \
	optimizer.cc optimizer.hh \
	optimizer-gsl.cc optimizer-gsl.hh
libeosoptimize_la_CXXFLAGS = $(AM_CXXFLAGS)
libeosoptimize_la_LIBADD =

include_eos_optimizedir = $(includedir)/eos/optimize
include_eos_optimize_HEADERS = \
	optimizer.hh

TESTS = \
	optimizer-gsl_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	libeosoptimize.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

optimizer_gsl_TEST_SOURCES = optimizer-gsl_TEST.cc
