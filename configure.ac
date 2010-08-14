dnl vim: set sw=8 sts=8 noet ft=config foldmethod=marker foldmarker={{{,}}}:

AC_INIT([wfitter], [0.0])
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

dnl {{{ output
AM_CONFIG_HEADER(config.h)
AC_OUTPUT(
	Makefile
	src/Makefile
	src/clients/Makefile
	src/rare-b-decays/Makefile
	src/utils/Makefile
	test/Makefile
	)
dnl }}}
