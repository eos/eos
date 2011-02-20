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
	AC_COMPILE_IFELSE([
#include <string>
#include <iostream>
template <typename T_> struct S { };
int main(int, char **)
{
	std::string s("test");
	std::cout << s << std::endl;
}
	],
		[AC_MSG_RESULT([yes])],
		[AC_MSG_RESULT([no])
		 AC_MSG_ERROR([your compiler does not support CXXFLAGS=]$1)
		])
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
AC_COMPILE_IFELSE([
	struct Foo { double x, y; };
	int main(int, char **)
	{
		Foo foo{ 17.24, 23.08 };
	}
],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([your compiler does not support braced initializers])
	]
)
dnl }}}
dnl }}}

dnl {{{ check for libraries
save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -lgslcblas" # work around Debian/Ubuntu oddities
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
CXXFLAGS=$save_CXXFLAGS
dnl }}}

dnl {{{ output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	doc/Makefile
	src/Makefile
	src/clients/Makefile
	src/rare-b-decays/Makefile
	src/utils/Makefile
	test/Makefile
	)
dnl }}}
