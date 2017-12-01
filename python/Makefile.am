CLEANFILES = *~

EXTRA_DIST = eosdata.py

AM_TESTS_ENVIRONMENT = \
       export EOS_TESTS_PARAMETERS="$(top_srcdir)/eos/parameters"; \
       export PYTHONPATH="$(top_builddir)/python/.libs/";

TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = 

if EOS_ENABLE_PYTHON

eosdir = $(pkgpythondir)
eos_PYTHON = 

pyexec_LTLIBRARIES = _eos.la
pyexec_SCRIPTS = eos.py eosdata.py

_eos_la_SOURCES = _eos.cc
_eos_la_CXXFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir) -std=c++0x
_eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex initeos
_eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python$(BOOST_PYTHON_SUFFIX)

TESTS = \
	eos_TEST.py

EXTRA_DIST += $(TESTS)

endif
