CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

noinst_LTLIBRARIES = libwfrarebdecays.la
libwfrarebdecays_la_SOURCES = \
	exclusive-b-to-s-dilepton.cc exclusive-b-to-s-dilepton.hh \
	form_factors.cc form_factors.hh
libwfrarebdecays_la_LIBADD = $(top_builddir)/src/utils/libwfutils.la

noinst_PROGRAMS = amplitudes ff scan-sm scan-rh observables
amplitudes_SOURCES = \
	amplitudes.cc
amplitudes_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
ff_SOURCES = \
	ff.cc
ff_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
observables_SOURCES = \
	observables.cc
observables_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan_sm_SOURCES = \
	scan-sm.cc
scan_sm_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
scan_rh_SOURCES = \
	scan-rh.cc
scan_rh_LDADD = libwfrarebdecays.la $(top_builddir)/src/utils/libwfutils.la
