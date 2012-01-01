CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

noinst_LIBRARIES = libeostest.a

libeostest_a_SOURCES = \
	test.cc test.hh

TESTS = \
	test_TEST
LDADD = libeostest.a \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cc
