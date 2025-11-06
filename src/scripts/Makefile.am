CLEANFILES = \
	*~

MAINTAINERCLEANFILES = Makefile.in

if EOS_ENABLE_PYTHON

bin_SCRIPTS = \
	eos-analysis \
	eos-data \
	eos-figure \
	eos-list-observables \
	eos-list-references \
	eos-make-constraint \
	eos-merge-mcmc \
	eos-plot \
	eos-plot-1d \
	eos-plot-2d \
	eos-print-modes \
	eos-print-uncertainty

AM_TESTS_ENVIRONMENT = \
			 export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
			 export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
			 export PYTHONPATH="$(top_builddir)/python/.libs/:$(top_srcdir)/python"; \
			 export PYTHON="$(PYTHON)"; \
			 export SOURCE_DIR="$(abs_srcdir)";

LOG_COMPILER="/bin/bash"
TESTS = \
	eos-analysis_TEST

endif

EXTRA_DIST = \
	eos-analysis \
	eos-data \
	eos-figure \
	eos-list-observables \
	eos-list-references \
	eos-make-constraint \
	eos-merge-mcmc \
	eos-plot \
	eos-plot-1d \
	eos-plot-2d \
	eos-print-modes \
	eos-print-uncertainty \
	eos-analysis_TEST \
	eos-analysis_TEST.d/analysis.yaml
