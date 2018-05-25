CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++14 -Wall -Wextra -pedantic

noinst_LIBRARIES = libcli.a

libcli_a_SOURCES = \
	cli_dumper.cc cli_dumper.hh \
	cli_error.cc cli_error.hh \
	cli_group.cc cli_group.hh \
	cli_handler.cc cli_handler.hh \
	cli_section.cc cli_section.hh \
	cli_option.cc cli_option.hh \
	cli_visitor.cc cli_visitor.hh

bin_PROGRAMS = \
	eos-evaluate \
	eos-find-mode \
	eos-list-constraints \
	eos-list-observables \
	eos-list-parameters \
	eos-list-signal-pdfs \
	eos-print-polynomial \
	eos-propagate-uncertainty \
	eos-sample-mcmc \
	eos-sample-events-mcmc \
	eos-scan-mc
noinst_PROGRAMS = \
	integrated \
	observables

LDADD = \
	$(top_builddir)/eos/statistics/libeosstatistics.la \
	$(top_builddir)/eos/optimize/libeosoptimize.la \
	$(top_builddir)/eos/utils/libeosutils.la \
	$(top_builddir)/eos/form-factors/libeosformfactors.la \
	$(top_builddir)/eos/b-decays/libeosbdecays.la \
	$(top_builddir)/eos/rare-b-decays/libeosrarebdecays.la \
	$(top_builddir)/eos/libeos.la \
	libcli.a \
	-lboost_filesystem -lboost_system \
	-lMinuit2 \
	-lyaml-cpp

if EOS_ENABLE_PMC
bin_PROGRAMS += \
	eos-sample-pmc

LDADD += \
	-lpmc -ldl

eos_sample_pmc_SOURCES = eos-sample-pmc.cc
eos_sample_pmc_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS)
eos_sample_pmc_LDADD = $(LDADD) $(HDF5_LDFLAGS) -lpmc -ldl
endif

eos_evaluate_SOURCES = eos-evaluate.cc

eos_find_mode_SOURCES = eos-find-mode.cc
eos_find_mode_CXXFLAGS = $(AM_CXXFLAGS) $(MINUIT2_CXXFLAGS)
eos_find_mode_LDADD = $(LDADD) $(MINUIT2_LDFLAGS)

eos_list_constraints_SOURCES = eos-list-constraints.cc

eos_list_observables_SOURCES = eos-list-observables.cc

eos_list_parameters_SOURCES = eos-list-parameters.cc

eos_list_signal_pdfs_SOURCES = eos-list-signal-pdfs.cc

eos_print_polynomial_SOURCES = eos-print-polynomial.cc

eos_propagate_uncertainty_SOURCES = eos-propagate-uncertainty.cc
eos_propagate_uncertainty_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS)
eos_propagate_uncertainty_LDADD = $(LDADD) $(HDF5_LDFLAGS)

eos_sample_mcmc_SOURCES = eos-sample-mcmc.cc

eos_sample_events_mcmc_SOURCES = eos-sample-events-mcmc.cc
eos_sample_events_mcmc_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS)
eos_sample_events_mcmc_LDADD = $(LDADD) $(HDF5_LDFLAGS)

eos_scan_mc_SOURCES = eos-scan-mc.cc
eos_scan_mc_CXXFLAGS = $(AM_CXXFLAGS) $(HDF5_CXXFLAGS) $(MINUIT2_CXXFLAGS)
eos_scan_mc_LDADD = $(LDADD) $(HDF5_LDFLAGS) $(MINUIT2_LDFLAGS)

integrated_SOURCES = integrated.cc

observables_SOURCES = observables.cc
