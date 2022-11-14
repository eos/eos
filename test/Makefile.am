CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = @AM_CXXFLAGS@
AM_LDFLAGS = @AM_LDFLAGS@

noinst_LTLIBRARIES = libeostest.la

libeostest_la_SOURCES = \
	test.cc test.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
	export EOS_TESTS_REFERENCES="$(top_srcdir)/eos/";

TESTS = \
	test_TEST
LDADD = libeostest.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cc
