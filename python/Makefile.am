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

pyexec_LTLIBRARIES = eos.la
pyexec_SCRIPTS = eosdata.py

eos_la_SOURCES = eos.cc
eos_la_CXXFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir) -std=c++0x
eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex initeos
eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python$(BOOST_PYTHON_SUFFIX)

TESTS = \
	eos_TEST.py

EXTRA_DIST += $(TESTS)

endif
