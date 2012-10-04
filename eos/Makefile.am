CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = utils form-factors rare-b-decays b-decays .

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	constraint.cc constraint.hh \
	observable.cc observable.hh
libeos_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la

include_eosdir = $(includedir)/eos
include_eos_HEADERS = \
	constraint.hh \
	observable.hh

TESTS = \
	constraint_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)

constraint_TEST_SOURCES = constraint_TEST.cc
