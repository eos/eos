dnl vim: set sw=8 sts=8 noet ft=config foldmethod=marker foldmarker={{{,}}}:

AC_INIT([EOS], [0.0])
AC_PREREQ(2.5)
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

dnl {{{ check for required programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_LIBTOOL
AC_PROG_MAKE_SET
dnl }}}

dnl {{{ use HDF5 path if supplied via --with-hdf5=
AC_ARG_WITH([hdf5],
	[AS_HELP_STRING([--with-hdf5], [used for I/O to scan data])]
	[],
	[])
if test "x$with_hdf5" == xno ; then
	AC_MSG_FAILURE(["Cannot build EOS without HDF5"])
fi
if test "x$with_hdf5" != x ; then
	if ! test -d "$with_hdf5" ; then
		AC_MSG_FAILURE(["$with_hdf5 is not a directory"])
	fi
	CXXFLAGS="-I$with_hdf5/include $CXXFLAGS"
	LDFLAGS="-L$with_hdf5/lib $LDFLAGS"
fi
dnl }}}
dnl {{{ use Minuit2 path if supplied via --with-minuit2=
AC_ARG_WITH([minuit2],
	[AS_HELP_STRING([--with-minuit2], [used to find starting points for population montecarlo sampling])]
	[],
	[])
if test "x$with_minuit2" == xno ; then
	AC_MSG_FAILURE(["Cannot build EOS without Minuit2"])
fi
if test "x$with_minuit2" != x ; then
	if ! test -d "$with_minuit2" ; then
		AC_MSG_FAILURE(["$with_minuit2 is not a directory"])
	fi
	CXXFLAGS="-I$with_minuit2/include $CXXFLAGS"
	LDFLAGS="-L$with_minuit2/lib $LDFLAGS"
fi
dnl }}}
dnl {{{ enable PMC path if supplied via --with-pmc=
AC_ARG_WITH([pmc],
	[AS_HELP_STRING([--with-pmc], [used to find perform population montecarlo sampling])]
	[],
	[])
if test "x$with_pmc" == xno -o "x$with_pmc" == x ; then
	AC_MSG_NOTICE(["Building EOS without PMC support"])
else 
	if ! test -d "$with_pmc" ; then
		AC_MSG_FAILURE(["$with_pmc is not a directory"])
	fi
	CXXFLAGS="-I$with_pmc/include $CXXFLAGS"
	LDFLAGS="-Wl,-rpath,$with_pm/lib -L$with_pmc/lib $LDFLAGS"
	AC_DEFINE([EOS_ENABLE_PMC], 1, [Whether we use PMC])
fi
AM_CONDITIONAL([EOS_ENABLE_PMC], [test "x$with_pmc" != "xno" -a "x$with_pmc" != x])
dnl }}}

dnl {{{ check for git HEAD
if test -d "${GIT_DIR:-${ac_top_srcdir:-./}/.git}" ; then
	GITHEAD=`git describe 2> /dev/null`
	if test -z ${GITHEAD} ; then
		GITHEAD=`git rev-parse --short=7 HEAD`
	fi
	if test -n "`git diff-index -m --name-only HEAD`" ; then
		GITHEAD=${GITHEAD}-dirty
	fi
else
	GITHEAD=$PACKAGE_VERSION
fi
AC_SUBST([GITHEAD])
AC_DEFINE_UNQUOTED([EOS_GITHEAD], "$GITHEAD", [GIT revision of the sources])
dnl }}}

dnl {{{ define CHECK_CXXFLAG
AC_DEFUN([CHECK_CXXFLAG], [
	save_CXXFLAGS=$CXXFLAGS
	CXXFLAGS="$CXXFLAGS $1 -Werror"
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([
#include <string>
#include <iostream>
template <typename T_> struct S { };
int main(int, char **)
{
	std::string s("test");
	std::cout << s << std::endl;
}
	])],
		[AC_MSG_RESULT([yes])],
		[AC_MSG_RESULT([no])
		 AC_MSG_ERROR([your compiler does not support CXXFLAGS=]$1)
		]
	)
	CXXFLAGS=$save_CXXFLAGS
])
dnl }}}

dnl {{{ check for compiler requirements
dnl {{{ check for -std=c++0x
AC_MSG_CHECKING([for C++0x support])
CHECK_CXXFLAG([-std=c++0x])
dnl }}}
dnl {{{ check for braced initializers
AC_MSG_CHECKING([for braced initializers])
save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -std=c++0x"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	struct Foo { double x, y; };
	int main(int, char **)
	{
		Foo foo{ 17.24, 23.08 };
	}
])],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([your compiler does not support braced initializers])
	]
)
CXXFLAGS=$save_CXXFLAGS
dnl }}}
dnl }}}

dnl {{{ check for libraries
save_LIBS=$LIBS
LIBS="$LIBS -lgslcblas" # work around Debian/Ubuntu oddities
dnl {{{ check for -lgsl
AC_CHECK_LIB(gsl, gsl_sf_zeta,
	[],
	[AC_MSG_ERROR([you do seem to be lacking the GNU Scientific Library (http://www.gnu.org/software/gsl/).])])
dnl }}}
dnl {{{ check for -lhdf5
AC_CHECK_LIB(hdf5, H5open,
	[],
	[AC_MSG_ERROR([you do seem to be lacking libhdf5 (http://hdf.ncsa.uiuc.edu/HDF5/).])])
dnl }}}
LIBS=$save_LIBS
dnl }}}

dnl {{{ define relevant EOS_* macros
AC_DEFINE_UNQUOTED([EOS_BUILDDIR], "$ac_pwd", [Absolute path to the build directory])
AC_DEFINE_UNQUOTED([EOS_SRCDIR], "`readlink -f $ac_pwd/$srcdir`", [Absolute path to the top-level source directory])
dnl }}}

dnl {{{ output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	doc/Makefile
	eos/Makefile
	eos/form-factors/Makefile
	eos/rare-b-decays/Makefile
        eos/b-decays/Makefile
	eos/utils/Makefile
	src/Makefile
	src/clients/Makefile
	test/Makefile
	)
dnl }}}
