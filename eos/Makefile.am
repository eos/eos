CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = utils rare-b-decays .

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	constraint.cc constraint.hh \
	observable.cc observable.hh
libeos_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	-lpthread -lgslcblas -lgsl

include_eosdir = $(includedir)/eos
include_eos_HEADERS = \
	constraint.hh \
	observable.hh

TESTS = \
	constraint_TEST
TESTS_LD_ADD = $(top_builddir)/test/libeostest.a $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la

check_PROGRAMS = $(TESTS)

constraint_TEST_SOURCES = constraint_TEST.cc
constraint_TEST_LDADD = $(TESTS_LD_ADD)
