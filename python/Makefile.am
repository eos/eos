CLEANFILES = *~

EXTRA_DIST = eosdata.py

if EOS_ENABLE_PYTHON

eosdir = $(pkgpythondir)
eos_PYTHON = 

pyexec_LTLIBRARIES = eos.la
pyexec_SCRIPTS = eosdata.py

eos_la_SOURCES = eos.cc
eos_la_CXXFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir) -std=c++0x
eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex initeos
eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python$(BOOST_PYTHON_SUFFIX)

AM_TESTS_ENVIRONMENT = PYTHONPATH="$(top_builddir)/python/.libs/"; export PYTHONPATH;
TEST_EXTENSIONS = .py
PY_LOG_COMPILER = $(PYTHON)
AM_PY_LOG_FLAGS = 

TESTS = \
	eos_TEST.py

EXTRA_DIST += $(TESTS)
endif
