CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

SUBDIRS = utils rare-b-decays . clients

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -Wall -Wextra -pedantic

lib_LTLIBRARIES = libeos.la
libeos_la_SOURCES = \
	factory.cc factory.hh \
	observable.cc observable.hh
libeos_la_LIBADD = \
	$(top_builddir)/src/utils/libeosutils.la \
	$(top_builddir)/src/rare-b-decays/libeosrarebdecays.la \
	-lpthread -lgslcblas -lgsl
