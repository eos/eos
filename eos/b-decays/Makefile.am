CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeosbdecays.la
libeosbdecays_la_SOURCES = \
	b-to-l-nu.cc b-to-l-nu.hh \
	b-to-pi-l-nu.cc b-to-pi-l-nu.hh \
	b-to-v-l-nu.hh \
	bs-to-kstar-l-nu.cc bs-to-kstar-l-nu.hh \
	properties.cc properties.hh
libeosbdecays_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	-lgslcblas \
	-lgsl

include_eos_bdecaysdir = $(includedir)/eos/b-decays
include_eos_bdecays_HEADERS = \
	properties.hh

EXTRA_DIST =

TESTS = \
	b-to-pi-l-nu_TEST \
	bs-to-kstar-l-nu_TEST \
	properties_TEST

LDADD = \
	$(top_builddir)/test/libeostest.a \
	libeosbdecays.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = $(TESTS)
b_to_pi_l_nu_TEST_SOURCES = b-to-pi-l-nu_TEST.cc

bs_to_kstar_l_nu_TEST_SOURCES = bs-to-kstar-l-nu_TEST.cc

properties_TEST_SOURCES = properties_TEST.cc
