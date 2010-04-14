CLEANFILES = *~
MAINTAINERCLEANFILES = Makefile.in

AM_CXXFLAGS = -I$(top_srcdir) -std=c++0x

lib_LTLIBRARIES = libwfutils.la
libwfutils_la_SOURCES = \
	accumulator.cc accumulator.hh \
	complex.hh \
	concrete_observable.hh \
	condition_variable.cc condition_variable.hh \
	destringify.hh \
	exception.cc exception.hh \
	integrate.cc integrate.hh \
	instantiation_policy.hh instantiation_policy-impl.hh \
	kinematic.cc kinematic.hh \
	lock.cc lock.hh \
	mutex.cc mutex.hh \
	parameters.cc parameters.hh \
	observable.cc observable.hh \
	private_implementation_pattern.hh private_implementation_pattern-impl.hh \
	qcd.cc qcd.hh \
	random_number_engine.cc random_number_engine.hh \
	thread.cc thread.hh \
	thread_pool.cc thread_pool.hh \
	ticket.cc ticket.hh
libwfutils_la_LIBADD = -lpthread

TESTS = \
	integrate_TEST \
	qcd_TEST
TESTS_LD_ADD = $(top_builddir)/test/libwftest.a libwfutils.la

check_PROGRAMS = $(TESTS)
integrate_TEST_SOURCES = integrate_TEST.cc
integrate_TEST_LDADD = $(TESTS_LD_ADD)

qcd_TEST_SOURCES = qcd_TEST.cc
qcd_TEST_LDADD = $(TESTS_LD_ADD)
