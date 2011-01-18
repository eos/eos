CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

bin_PROGRAMS = \
	eos-evaluate \
	eos-list-parameters \
	eos-print-polynomial \
	eos-scan-polynomial \
	eos-scan

eos_evaluate_SOURCES = \
	eos-evaluate.cc
eos_evaluate_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la

eos_list_parameters_SOURCES = \
	eos-list-parameters.cc
eos_list_parameters_LDADD = \
	$(top_builddir)/src/utils/libeosutils.la

eos_print_polynomial_SOURCES = \
	eos-print-polynomial.cc
eos_print_polynomial_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la

eos_scan_polynomial_SOURCES = \
	eos-scan-polynomial.cc
eos_scan_polynomial_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la

eos_scan_SOURCES = \
	eos-scan.cc
eos_scan_LDADD = \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/src/utils/libeosutils.la
