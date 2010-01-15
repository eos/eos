CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

noinst_LTLIBRARIES = libwflowrecoil.la
libwflowrecoil_la_SOURCES = \
	low-recoil.cc low-recoil.hh

noinst_PROGRAMS = observables
observables_SOURCES = \
	observables.cc
observables_LDADD = libwflowrecoil.la $(top_builddir)/src/utils/libwfutils.la
