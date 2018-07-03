dnl vim: set sw=8 sts=8 noet ft=config foldmethod=marker foldmarker={{{,}}}:

AC_INIT([EOS], [0.2.1])
AC_PREREQ(2.5)
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(1.9)
AC_LANG([C++])

dnl {{{ check for required programs
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_LN_S
LT_INIT([disable-static])
AC_PROG_MAKE_SET
PKG_PROG_PKG_CONFIG
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
dnl {{{ check for -std=c++14
AC_MSG_CHECKING([for C++14 support])
CHECK_CXXFLAG([-std=c++14])
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
dnl {{{ check for constexpr const
AC_MSG_CHECKING([for constexpr const])
save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -std=c++0x"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	struct X
	{
		static constexpr const int foo = 42;
	};
	int main(int, char **)
	{
		return 42 == X::foo;
	}
])],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([your compiler does not support constexpr const])
	]
)
CXXFLAGS=$save_CXXFLAGS
dnl }}}
CXXFLAGS="$CXXFLAGS -std=c++14"
dnl }}}

dnl {{{ check for required packages
dnl {{{ boost filesystem
AC_CHECK_LIB([boost_filesystem], main,
	[],
	[AC_MSG_ERROR([you do seem to be lacking 'libboost_filesystem'.])]
	)
dnl }}}
dnl {{{ boost system
AC_CHECK_LIB([boost_system], main,
	[],
	[AC_MSG_ERROR([you do seem to be lacking 'libboost_system'.])]
	)
dnl }}}
dnl {{{ yaml-cpp
AC_CHECK_LIB([yaml-cpp], main)
AC_CHECK_HEADERS([yaml-cpp/yaml.h],
        [eos_found_yaml_cpp_headers=yes; break;])

AS_IF([test "x$eos_found_yaml_cpp_headers" != "xyes"],
        [AC_MSG_ERROR([you do seem to be lacking 'libyaml-cpp-dev'.])])
dnl }}}
dnl }}}

dnl {{{ use HDF5 path if supplied via --with-hdf5=, otherwise use pkg-config; check for HDF5
AC_ARG_WITH([hdf5],
	[AS_HELP_STRING([--with-hdf5], [used for I/O to scan data])]
	[],
	[])
if test "x$with_hdf5" == xno ; then
	AC_MSG_FAILURE(["Cannot build EOS without HDF5"])
fi
if test "x$with_hdf5" != x -a "x$with_hdf5" != xyes ; then
	if ! test -d "$with_hdf5" ; then
		AC_MSG_FAILURE(["$with_hdf5 is not a directory. Please check the manual for how to compile with HDF5"])
	fi
	HDF5_CXXFLAGS="-I$with_hdf5/include"
	HDF5_LDFLAGS="-L$with_hdf5/lib"
else # use pkg-config
	PKG_CHECK_MODULES([HDF5],
		[hdf5],
		[
			HDF5_CXXFLAGS="$HDF5_CFLAGS"
			HDF5_LDFLAGS="$HDF5_LIBS"
			AC_MSG_NOTICE([using pkg-config for HDF5 location])
		],
		[
			HDF5_CXXFLAGS=""
			HDF5_LDFLAGS="-lhdf5"
			AC_MSG_NOTICE([falling back to default for HDF5 location])
		])
fi

AC_SUBST([HDF5_LDFLAGS])
AC_SUBST([HDF5_CXXFLAGS])

_CXXFLAGS=$CXXFLAGS
_LDFLAGS=$LDFLAGS
_LIBS=$LIBS

CXXFLAGS="$HDF5_CXXFLAGS $CXXFLAGS"
LDFLAGS="$HDF5_LDFLAGS $LDFLAGS"
LIBS=""

AC_CHECK_LIB(hdf5, H5open,
	[],
	[AC_MSG_ERROR([you do seem to be lacking libhdf5 (http://hdf.ncsa.uiuc.edu/HDF5/).])],
	[]
	)

CXXFLAGS=$_CXXFLAGS
LDFLAGS=$_LDFLAGS
LIBS="$_LIBS"
dnl }}}

dnl {{{ use GSL path if supplied via --with-gsl=, otherwise use pkg-config; check for GSL
AC_ARG_WITH([gsl],
	[AS_HELP_STRING([--with-gsl], [used for e.g. RNG, special functions, linear algebra])]
	[],
	[])
if test "x$with_gsl" == xno ; then
	AC_MSG_FAILURE(["Cannot build EOS without GSL"])
fi
if test "x$with_gsl" != x -a "x$with_gsl" != xyes ; then
	if ! test -d "$with_gsl" ; then
	AC_MSG_FAILURE(["$with_gsl is not a directory. Please check the manual for how to compile with GSL."])
	fi
	GSL_CXXFLAGS="-I$with_gsl/include"
	GSL_LDFLAGS="-L$with_gsl/lib"
else # use pkg-config
	PKG_CHECK_MODULES([GSL],
		[gsl],
		[
			GSL_CXXFLAGS="$GSL_CFLAGS"
			GSL_LDFLAGS="`pkg-config --libs-only-L gsl`"
			AC_MSG_NOTICE([using pkg-config for GSL location])
		],
		[
			GSL_CXXFLAGS=""
			GSL_LDFLAGS=""
			AC_MSG_NOTICE([falling back to default for GSL location])
		])
fi
_CXXFLAGS=$CXXFLAGS
_LDFLAGS=$LDFLAGS
_LIBS=$LIBS

CXXFLAGS="$GSL_CXXFLAGS $CXXFLAGS"
LDFLAGS="$GSL_LDFLAGS $LDFLAGS"
LIBS=""

AC_CHECK_LIB(gsl, gsl_sf_zeta,
	[],
	[AC_MSG_ERROR([you do seem to be lacking the GNU Scientific Library (http://www.gnu.org/software/gsl/).])],
	[-lgslcblas -lm])
CXXFLAGS=$_CXXFLAGS
LDFLAGS=$_LDFLAGS
LIBS="$LIBS $_LIBS"
dnl }}}

dnl {{{ use Minuit2 path if supplied via --with-minuit2=
AC_ARG_WITH([minuit2],
	[AS_HELP_STRING([--with-minuit2], [used to find starting points for population Monte Carlo sampling])]
	[],
	[])
if test "x$with_minuit2" == xroot ; then
	AC_CHECK_PROG(ROOTCONFIG, root-config, yes)
	if test "x$ROOTCONFIG" == xyes ; then
		MINUIT2_CXXFLAGS="-I`root-config --incdir` $CXXFLAGS"
		MINUIT2_LDFLAGS="-L`root-config --libdir` $LDFLAGS"
	else
		AC_MSG_FAILURE(["Cannot find root-config"])
	fi
elif test "x$with_minuit2" != xno -a "x$with_minuit2" != x ; then
	if ! test -d "$with_minuit2" ; then
	AC_MSG_FAILURE(["$with_minuit2 is not a directory. Please check the manual for how to compile with Minuit2"])
	fi
	CXXFLAGS="-I$with_minuit2/include $CXXFLAGS"
	LDFLAGS="-L$with_minuit2/lib $LDFLAGS"
fi
AC_SUBST([MINUIT2_LDFLAGS])
AC_SUBST([MINUIT2_CXXFLAGS])
dnl }}}
dnl {{{ use PMC path if supplied via --with-pmc=
AC_ARG_WITH([pmc],
	[AS_HELP_STRING([--with-pmc], [used to perform population Monte Carlo sampling])]
	[],
	[])
if test "x$with_pmc" != x ; then
	if ! test -d "$with_pmc" ; then
	AC_MSG_FAILURE(["$with_pmc is not a directory. Please check the manual for how to compile with pmclib"])
	else
		CXXFLAGS="-I$with_pmc/include $CXXFLAGS"
		LDFLAGS="-Wl,-rpath,$with_pmc/lib -L$with_pmc/lib $LDFLAGS"
	fi
fi
dnl }}}

dnl {{{ check if we should enable PMC
AC_ARG_ENABLE([pmc],
	[AS_HELP_STRING([--enable-pmc], [enables usage of pmclib by EOS])],
	[],
	[])
if test "x$enable_pmc" == "xyes" ; then
	AC_DEFINE([EOS_ENABLE_PMC], [1], [Define if PMC is enabled.])
fi
AM_CONDITIONAL([EOS_ENABLE_PMC], [test "x$enable_pmc" == "xyes"])
dnl }}}
dnl {{{ check if we should enable Python
AC_ARG_ENABLE([python],
	[AS_HELP_STRING([--enable-python], [enables interfacing EOS within python])],
	[],
	[])
AM_CONDITIONAL([EOS_ENABLE_PYTHON], [test "x$enable_python" == "xyes"])
AC_ARG_WITH([boost-python-suffix],
	[AS_HELP_STRING([--with-boost-python-suffix], [suffix used to specify the libboost-python shared object])]
	[],
	[])
dnl }}}
dnl {{{ python checks (you can change the required python version below)
if test "x$enable_python" == xyes ; then
	AM_PATH_PYTHON([3.5])
	PY_PREFIX=`$PYTHON -c 'import sys ; print(sys.prefix)'`
	PYTHON_LIBS="-lpython$PYTHON_VERSION"
	PYTHON_LDFLAGS="`python$PYTHON_VERSION-config --ldflags`"
	PYTHON_CXXFLAGS="`python$PYTHON_VERSION-config --includes`"
	os_id="`. /etc/os-release ; echo $ID`"
	if test "x$os_id" == xubuntu -o "x$os_id" == xdebian ; then
		BOOST_PYTHON_SUFFIX="-py`echo $PYTHON_VERSION | sed -e 's/\.//'`"
	elif test "x$os_id" == xcentos ; then
		BOOST_PYTHON_SUFFIX=""
	elif test "x`uname`" == xDarwin ; then
		BOOST_PYTHON_SUFFIX=$with_boost_python_suffix
	fi
	AC_SUBST([PYTHON_LIBS])
	AC_SUBST([PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_CXXFLAGS])
	AC_SUBST([BOOST_PYTHON_SUFFIX])

	_CXXFLAGS=$CXXFLAGS
	_LDFLAGS=$LDFLAGS
	_LIBS=$LIBS

	CXXFLAGS=$PYTHON_CXXFLAGS
	LDFLAGS=$PYTHON_LDFLAGS
	LIBS=$PYTHON_LIBS

	AC_CHECK_LIB([boost_python$BOOST_PYTHON_SUFFIX], main,
	        [],
	        [AC_MSG_ERROR([you do seem to be lacking 'libboost_python'.])]
	)

	CXXFLAGS=$_CXXFLAGS
	LDFLAGS=$_LDFLAGS
	LIBS=$_LIBS
fi
dnl }}}
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

dnl {{{ define relevant EOS_* macros
AC_DEFINE_UNQUOTED([EOS_BUILDDIR], "$ac_pwd", [Absolute path to the build directory])
AC_DEFINE_UNQUOTED([EOS_SRCDIR], "$ac_abs_confdir", [Absolute path to the top-level source directory])
eval ac_datadir=`eval x=$datadir ; echo $x`
AC_DEFINE_UNQUOTED([EOS_DATADIR], "$ac_datadir", [Absolute path to the installation data directory])
dnl }}}

dnl {{{ output
AM_CONFIG_HEADER(config.h)
AM_EXTRA_RECURSIVE_TARGETS([manual])
AC_SUBST([AM_CXXFLAGS], '-I$(top_srcdir) -std=c++14 -Wall -Wextra -pedantic')
AC_OUTPUT(
	Makefile
	debian/control-bionic
	debian/control-xenial
	debian/Makefile
	doc/Makefile
	eos/Makefile
	eos/constraints/Makefile
	eos/form-factors/Makefile
	eos/rare-b-decays/Makefile
	eos/b-decays/Makefile
	eos/parameters/Makefile
	eos/utils/Makefile
	eos/optimize/Makefile
	eos/statistics/Makefile
	manual/Makefile
	manual/appendices/Makefile
	manual/figures/Makefile
	src/Makefile
	src/clients/Makefile
	src/scripts/Makefile
	test/Makefile
	python/Makefile
	)
dnl }}}
