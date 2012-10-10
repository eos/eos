CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

bin_PROGRAMS = \
	eos-evaluate \
	eos-find-confidence-region \
	eos-inspect \
	eos-list-parameters \
	eos-merge \
	eos-marginalise \
	eos-print-polynomial \
	eos-propagate-uncertainty \
	eos-realify \
	eos-scan-mc \
	eos-scan-polynomial \
	eos-scan
noinst_PROGRAMS = \
	integrated \
	observables

LDADD = \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/libeos.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la

eos_evaluate_SOURCES = eos-evaluate.cc

eos_find_confidence_region_SOURCES = eos-find-confidence-region.cc

eos_inspect_SOURCES = eos-inspect.cc

eos_list_parameters_SOURCES = eos-list-parameters.cc

eos_marginalise_SOURCES = eos-marginalise.cc

eos_realify_SOURCES = eos-realify.cc

eos_merge_SOURCES = eos-merge.cc

eos_print_polynomial_SOURCES = eos-print-polynomial.cc

eos_propagate_uncertainty_SOURCES = eos-propagate-uncertainty.cc

eos_scan_mc_SOURCES = eos-scan-mc.cc
eos_scan_mc_LDADD = $(LDADD) -lgsl -lMinuit2

eos_scan_polynomial_SOURCES = eos-scan-polynomial.cc

eos_scan_SOURCES = eos-scan.cc

integrated_SOURCES = integrated.cc

observables_SOURCES = observables.cc
