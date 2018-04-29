CLEANFILES = \
	*~ \
	references_TEST
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = constraints parameters utils optimize statistics form-factors rare-b-decays b-decays .

AM_CXXFLAGS = @AM_CXXFLAGS@

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	constraint.cc constraint.hh \
	decays.hh \
	observable.cc observable.hh \
	signal-pdf.cc signal-pdf.hh
libeos_la_LIBADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la

include_eosdir = $(includedir)/eos
include_eos_HEADERS = \
	constraint.hh \
	decays.hh \
	observable.hh \
	signal-pdf.hh

AM_TESTS_ENVIRONMENT = \
	export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
	export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters";

TESTS = \
	constraint_TEST \
	references_TEST
LDADD = \
	$(top_builddir)/test/libeostest.a \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/libeos.la

check_PROGRAMS = \
	constraint_TEST \
	references_TEST

constraint_TEST_SOURCES = constraint_TEST.cc

references_TEST_SOURCES = references_TEST.in

pkgdata_DATA = references.yaml
EXTRA_DIST = \
	references.yaml \
	references_TEST.in

references_TEST$(EXEEXT): references_TEST.in Makefile
	sed \
	    -e "s#@TOPSRCDIR@#$(top_srcdir)#g" \
	    $< \
	    > $@
	chmod a+x $@
