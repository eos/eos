CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

noinst_LTLIBRARIES = libeostest.la

libeostest_la_SOURCES = \
	test.cc test.hh

TESTS = \
	test_TEST
LDADD = libeostest.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cc
