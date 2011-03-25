CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

noinst_LIBRARIES = libeostest.a

libeostest_a_SOURCES = \
	test.cc test.hh

TESTS = \
	test_TEST
TESTS_LD_ADD = libeostest.a -lgslcblas

check_PROGRAMS = $(TESTS)

test_TEST_SOURCES = test_TEST.cc
test_TEST_LDADD = $(TESTS_LD_ADD)
