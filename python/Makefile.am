CLEANFILES = *~

if EOS_ENABLE_PYTHON

eosdir = $(pkgpythondir)
eos_PYTHON = 

pyexec_LTLIBRARIES = eos.la

eos_la_SOURCES = eos.cc
eos_la_CXXFLAGS = $(PYTHON_CFLAGS) -I$(top_srcdir) -std=c++0x
eos_la_LDFLAGS = -module -avoid-version -export-symbols-regex initeos
eos_la_LIBADD = $(top_builddir)/eos/libeos.la $(top_builddir)/eos/utils/libeosutils.la -lboost_python

EXTRA_DIST = 

endif
