CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

noinst_LTLIBRARIES = libwfrarebdecays.la
libwfrarebdecays_la_SOURCES = \
	exclusive-b-to-s-dilepton.cc exclusive-b-to-s-dilepton.hh
libwfrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libwfutils.la

noinst_PROGRAMS = amplitudes scan observables
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
observables_SOURCES = \
	observables.cc
observables_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan_SOURCES = \
	scan.cc
scan_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
