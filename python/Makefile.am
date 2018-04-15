CLEANFILES = *~

EXTRA_DIST = \
	eos/__init__.py \
	eos/ipython.py \
	eos/data/__init__.py \
	eos/plot/__init__.py \
	eos/plot/config.py \
	eos/plot/plotter.py

AM_TESTS_ENVIRONMENT = \
			 export EOS_TESTS_CONSTRAINTS="$(top_srcdir)/eos/constraints"; \
			 export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
			 export PYTHONPATH="$(top_builddir)/python/.libs/:$(top_builddir)/python";

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS =

if EOS_ENABLE_PYTHON

eosdir = $(pkgpythondir)
eos_PYTHON =
eos_SCRIPTS = \
	eos/__init__.py \
	eos/ipython.py

eosdatadir = $(pkgpythondir)/data
eosdata_SCRIPTS = \
	eos/data/__init__.py

eosplotdir = $(pkgpythondir)/plot
eosplot_SCRIPTS = \
	eos/plot/__init__.py \
	eos/plot/config.py \
	eos/plot/plotter.py

_eosdir = $(pythondir)
_eos_LTLIBRARIES = _eos.la
_eos_la_SOURCES = _eos.cc
_eos_la_CXXFLAGS = $(PYTHON_CXXFLAGS) @AM_CXXFLAGS@
_eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex PyInit__eos
_eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python$(BOOST_PYTHON_SUFFIX)

TESTS = \
	eos_TEST.py

EXTRA_DIST += $(TESTS)

endif
