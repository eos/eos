CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

lib_LTLIBRARIES = libwfutils.la
libwfutils_la_SOURCES = \
	accumulator.cc accumulator.hh \
	complex.cc complex.hh complex-impl.hh \
	concrete_observable.hh \
	condition_variable.cc condition_variable.hh \
	exception.cc exception.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	kinematic.cc kinematic.hh \
	lock.cc lock.hh \
	mutex.cc mutex.hh \
	parameters.cc parameters.hh \
	observable.cc observable.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	thread.cc thread.hh \
	thread_pool.cc thread_pool.hh \
	ticket.cc ticket.hh
libwfutils_la_LIBADD = -lpthread

TESTS = \
	complex_TEST
TESTS_LD_ADD = $(top_builddir)/test/libwftest.a libwfutils.la

check_PROGRAMS = $(TESTS)
complex_TEST_SOURCES = complex_TEST.cc
complex_TEST_LDADD = $(TESTS_LD_ADD)
