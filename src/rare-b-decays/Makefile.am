CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x -ggdb

noinst_LTLIBRARIES = libwfrarebdecays.la
libwfrarebdecays_la_SOURCES = \
	exclusive-b-to-s-dilepton.cc exclusive-b-to-s-dilepton.hh \
	factory.cc factory.hh \
	form_factors.cc form_factors.hh
libwfrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libwfutils.la

noinst_PROGRAMS = amplitudes ff scan2 scan3 observables
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
ff_SOURCES = \
	ff.cc
ff_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
observables_SOURCES = \
	observables.cc
observables_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan2_SOURCES = \
	scan2.cc
scan2_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan3_SOURCES = \
	scan3.cc
scan3_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
