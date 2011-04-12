CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

bin_PROGRAMS = \
	eos-evaluate \
	eos-inspect \
	eos-list-parameters \
	eos-merge \
	eos-marginalise \
	eos-print-polynomial \
	eos-realify \
	eos-scan-polynomial \
	eos-scan
noinst_PROGRAMS = \
	integrated \
	observables

EOS_CLIENTS_LIBS = $(top_builddir)/src/libeos.la $(top_builddir)/src/utils/libeosutils.la

eos_evaluate_SOURCES = \
	eos-evaluate.cc
eos_evaluate_LDADD = $(EOS_CLIENTS_LIBS)

eos_inspect_SOURCES = \
	eos-inspect.cc
eos_inspect_LDADD = $(EOS_CLIENTS_LIBS)

eos_list_parameters_SOURCES = \
	eos-list-parameters.cc
eos_list_parameters_LDADD = $(EOS_CLIENTS_LIBS)

eos_marginalise_SOURCES = \
	eos-marginalise.cc
eos_marginalise_LDADD = $(EOS_CLIENTS_LIBS)

eos_realify_SOURCES = \
	eos-realify.cc
eos_realify_LDADD = $(EOS_CLIENTS_LIBS)

eos_merge_SOURCES = \
	eos-merge.cc
eos_merge_LDADD = $(EOS_CLIENTS_LIBS)

eos_print_polynomial_SOURCES = \
	eos-print-polynomial.cc
eos_print_polynomial_LDADD = $(EOS_CLIENTS_LIBS)

eos_scan_polynomial_SOURCES = \
	eos-scan-polynomial.cc
eos_scan_polynomial_LDADD = $(EOS_CLIENTS_LIBS)

eos_scan_SOURCES = \
	eos-scan.cc
eos_scan_LDADD = $(EOS_CLIENTS_LIBS)

integrated_SOURCES = \
	integrated.cc
integrated_LDADD = $(EOS_CLIENTS_LIBS)

observables_SOURCES = \
	observables.cc
observables_LDADD = $(EOS_CLIENTS_LIBS)
