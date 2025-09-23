dnl vim: set sw=8 sts=8 noet ft=config foldmethod=marker foldmarker={{{,}}}:

AC_INIT([EOS],[1.0.17])
AC_PREREQ([2.69])
AC_CONFIG_SRCDIR([])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE([1.9 subdir-objects])
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
dnl {{{ check for -std=c++20
AC_MSG_CHECKING([for C++20 support])
CHECK_CXXFLAG([-std=c++20])
dnl }}}
dnl {{{ check for braced initializers
AC_MSG_CHECKING([for braced initializers])
save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -std=c++20"
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
CXXFLAGS="$CXXFLAGS -std=c++20"
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
CXXFLAGS="$CXXFLAGS -std=c++20"
dnl {{{ check for <source_location>
AC_MSG_CHECKING([for <source_location>])
save_CXXFLAGS=$CXXFLAGS
CXXFLAGS="$CXXFLAGS -std=c++20"
AC_COMPILE_IFELSE([AC_LANG_SOURCE([
	#ifdef __clang__
	#  include <experimental/source_location>
	#else
	#  include <source_location>
	#endif

	#ifdef __clang__
		using source_location = std::experimental::source_location;
	#else
		using source_location = std::source_location;
	#endif

	int main(int, char **)
	{
		auto loc = source_location::current();
		return 0;
	}
])],
	[AC_MSG_RESULT([yes])],
	[AC_MSG_RESULT([no])
	 AC_MSG_ERROR([your compiler does not support <source_location>])
	]
)
CXXFLAGS=$save_CXXFLAGS
dnl }}}
dnl }}}

dnl {{{ determine which pieces of EOS to build
AC_ARG_ENABLE([cli],
	[AS_HELP_STRING([--enable-cli],
	    [build command-line utilities (default is yes)])]
)
if test "x$enable_cli" = xno ; then
    MAYBE_SRC=
else
    MAYBE_SRC=src
fi
AC_SUBST([MAYBE_SRC])
dnl }}}

dnl {{{ check for required packages
dnl {{{ boost filesystem
AC_CHECK_LIB([boost_filesystem], main,
	[],
	[AC_MSG_ERROR([you do seem to be lacking 'libboost_filesystem'.])]
	)
dnl }}}
dnl {{{ yaml-cpp
PKG_CHECK_MODULES([YAMLCPP],
	[yaml-cpp],
	[
		YAMLCPP_CXXFLAGS="$YAMLCPP_CFLAGS"
		YAMLCPP_LDFLAGS="$YAMLCPP_LIBS"
		AC_MSG_NOTICE([using pkg-config for yaml-cpp location])
	],
	[
		YAMLCPP_CXXFLAGS=""
		YAMLCPP_LDFLAGS="-lyaml-cpp"
		AC_MSG_NOTICE([falling back to default for yaml-cpp location])
	])

AC_SUBST([YAMLCPP_CXXFLAGS])
AC_SUBST([YAMLCPP_LDFLAGS])

_CXXFLAGS=$CXXFLAGS
_LDFLAGS=$LDFLAGS
_LIBS=$LIBS

CXXFLAGS="$YAMLCPP_CXXFLAGS $CXXFLAGS"
LDFLAGS="$YAMLCPP_LDFLAGS $LDFLAGS"
LIBS=""

AC_CHECK_LIB([yaml-cpp], main)
AC_CHECK_HEADERS([yaml-cpp/yaml.h],
        [eos_found_yaml_cpp_headers=yes; break;])

AS_IF([test "x$eos_found_yaml_cpp_headers" != "xyes"],
        [AC_MSG_ERROR([you do seem to be lacking 'libyaml-cpp-dev'.])])

CXXFLAGS=$_CXXFLAGS
LDFLAGS=$_LDFLAGS
LIBS=$_LIBS
dnl }}}
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
AC_SUBST([GSL_LDFLAGS])
AC_SUBST([GSL_CXXFLAGS])

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

AC_CHECK_LIB(gsl, gsl_linalg_cholesky_decomp1,
	[
		EOS_USE_GSL_LINALG_CHOLESKY_DECOMP="0"
	],
	[
		EOS_USE_GSL_LINALG_CHOLESKY_DECOMP="1"
		AC_MSG_NOTICE([falling back to legace GSL cholesky decomposition])
	],
	[])

AC_DEFINE_UNQUOTED([EOS_USE_GSL_LINALG_CHOLESKY_DECOMP], $EOS_USE_GSL_LINALG_CHOLESKY_DECOMP, [Whether to use legacy gsl_linalg_cholesky_decomp])

CXXFLAGS=$_CXXFLAGS
LDFLAGS=$_LDFLAGS
LIBS="$LIBS $_LIBS"
AC_SUBST([GSL_LDFLAGS])
AC_SUBST([GSL_CXXFLAGS])
dnl }}}

dnl {{{ use custom pythondir if supplied via --with-custom-pythondir=
AC_ARG_WITH([custom-pythondir],
	[AS_HELP_STRING([--with-custom-pythondir], [used to override Debian/Ubuntu default choice for \${pythondir} and \${pyexecdir}])],
	[],
	[])
dnl }}}
dnl {{{ check if we should enable Python
AC_ARG_ENABLE([python],
	[AS_HELP_STRING([--enable-python], [enables interfacing EOS within python])],
	[],
	[enable_python=yes])
AM_CONDITIONAL([EOS_ENABLE_PYTHON], [test "x$enable_python" == "xyes"])
has_boost_python_suffix=no
AC_ARG_WITH([boost-python-suffix],
	[AS_HELP_STRING([--with-boost-python-suffix], [suffix used to specify the libboost-python shared object])],
	[has_boost_python_suffix=yes],
	[with_boost_python_suffix=""])
dnl }}}
dnl {{{ python checks (you can change the required python version below)
if test "x$enable_python" == xyes ; then
	AM_PATH_PYTHON([3.10])
	PY_PREFIX=`$PYTHON -c 'import sys ; print(sys.prefix)'`
	PYTHON_LIBS="`python$PYTHON_VERSION-config --libs`"
	if python$PYTHON_VERSION-config --libs --embed 2>&1 > /dev/null ; then
		PYTHON_LIBS_CONFTEST=`python$PYTHON_VERSION-config --libs --embed`
	else
		PYTHON_LIBS_CONFTEST=$PYTHON_LIBS
	fi
	PYTHON_LDFLAGS="`python$PYTHON_VERSION-config --ldflags`"
	PYTHON_CXXFLAGS="`python$PYTHON_VERSION-config --includes`"
	if test "x$has_boost_python_suffix" == "xno"; then
		os_id="`. /etc/os-release ; echo $ID`"
		if test "x$os_id" == xubuntu -o "x$os_id" == xdebian ; then
			codename="`. /etc/os-release ; echo $VERSION_CODENAME`"
			case $codename in
				xenial|bionic)
					BOOST_PYTHON_SUFFIX="-py`echo $PYTHON_VERSION | sed -e 's/\.//'`"
					;;
				disco)
					BOOST_PYTHON_SUFFIX="`echo $PYTHON_VERSION | sed -e 's/\.//'`"
					;;
				*)
					AC_MSG_NOTICE([assuming Ubuntu Disco or more recent, BOOST_PYTHON_SUFFIX might be wrong!])
					BOOST_PYTHON_SUFFIX="`echo $PYTHON_VERSION | sed -e 's/\.//'`"
					;;
			esac
		elif test "x$os_id" == xcentos ; then
			BOOST_PYTHON_SUFFIX=""
		fi
	else
		BOOST_PYTHON_SUFFIX=$with_boost_python_suffix
		PYTHON_CXXFLAGS="$PYTHON_CXXFLAGS $BOOST_PYTHON_CXXFLAGS"
		PYTHON_LDFLAGS="$PYTHON_LDFLAGS $BOOST_PYTHON_LDFLAGS"
	fi
	AC_SUBST([PYTHON_LIBS])
	AC_SUBST([PYTHON_LDFLAGS])
	AC_SUBST([PYTHON_CXXFLAGS])
	AC_SUBST([BOOST_PYTHON_SUFFIX])
	AC_SUBST([BOOST_PYTHON_LDFLAGS])
	AC_SUBST([BOOST_PYTHON_CXXFLAGS])
	if test "x$with_custom_pythondir" != "x"; then
		# override: Debian/Ubuntu use 'dist-packages' rather than 'site-packages'
		AC_MSG_NOTICE([replacing python script directory with \${prefix}/lib/python$PYTHON_VERSION/$with_custom_pythondir])
		pythondir="\${prefix}/lib/python$PYTHON_VERSION/$with_custom_pythondir"
		AC_MSG_NOTICE([replacing python extension directory with \${exec_prefix}/lib/python$PYTHON_VERSION/$with_custom_pythondir])
		pyexecdir="\${exec_prefix}/lib/python$PYTHON_VERSION/$with_custom_pythondir"
		AC_SUBST([pythondir])
		AC_SUBST([pyexecdir])
	fi

	_CXXFLAGS=$CXXFLAGS
	_LDFLAGS=$LDFLAGS
	_LIBS=$LIBS

	CXXFLAGS=$PYTHON_CXXFLAGS
	LDFLAGS=$PYTHON_LDFLAGS
	LIBS=$PYTHON_LIBS_CONFTEST

	if test "x$BOOST_PYTHON_SKIP" == "xyes"; then
	AC_MSG_NOTICE([skipping checks for -lboost_python$BOOST_PYTHON_SUFFIX])
	PYTHON_LDFLAGS="$BOOST_PYTHON_LDFLAGS"
	PYTHON_LIBS="$BOOST_PYTHON_LIBS"
	else
	AC_CHECK_LIB([boost_python$BOOST_PYTHON_SUFFIX], main,
	        [],
	        [AC_MSG_ERROR([you do seem to be lacking 'libboost_python'.])]
	)

	AC_MSG_CHECKING([for working boost_python])
	AC_COMPILE_IFELSE([AC_LANG_SOURCE([
		#include <boost/python.hpp>
		char const* greet() {
			return "Hello, world!";
		}
		BOOST_PYTHON_MODULE(hello)
		{
			boost::python::def("greet", greet);
		}
	])],
		[AC_MSG_RESULT([yes])],
		[AC_MSG_RESULT([no])
		 AC_MSG_ERROR([boost_python example does not compile])
		]
	)
	fi

	CXXFLAGS=$_CXXFLAGS
	LDFLAGS=$_LDFLAGS
	LIBS=$_LIBS
fi
dnl }}}
dnl }}}

dnl {{{
AC_ARG_ENABLE([prerelease],
	[AS_HELP_STRING([--enable-prerelease], [enable prerelease and adjust version information])],
	[
		AC_MSG_CHECKING([for PYPI_VERSION])
		if test "x$enableval" == "xno" ; then
			PYPI_VERSION=$PACKAGE_VERSION
		else
			pypi_version_main=`python3 -c "v = [[int(n) for n in '${PACKAGE_VERSION}'.split('.')]] ; v += [[0]] * (3 - [len(v)]) ; v[[-1]] += 1 ; v = [[str(n) for n in v]] ; print('.'.join(v))"`
			PYPI_VERSION=$pypi_version_main.dev$enableval
		fi
		AC_MSG_RESULT([$PYPI_VERSION])
	],
	[PYPI_VERSION=$PACKAGE_VERSION]
	)
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
dnl }}}

dnl {{{ use link-time optimization (LTO) to reduce the size of the binaries
AC_ARG_ENABLE([lto],
	[AS_HELP_STRING([--enable-lto], [enables link-time optimiztation to reduce the size of the binaries])]
	[],
	[])
dnl }}}

dnl {{{ define default CXXFLAGS and LDFLAGS
AM_CXXFLAGS="-I$ac_abs_confdir -std=c++20 -Wall -Wextra -pedantic"
AM_LDFLAGS=""
if test "x$enable_lto" == xyes ; then
    AM_CXXFLAGS="$AM_CXXFLAGS -flto"
    AM_LDFLAGS="$AM_LDFLAGS -flto"
fi
dnl }}}

dnl {{{ output
AC_CONFIG_HEADERS(config.h)
AM_EXTRA_RECURSIVE_TARGETS([manual])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])
AC_SUBST([PYPI_VERSION])
AC_CONFIG_FILES([Makefile
	debian/control-focal
	debian/control-jammy
	debian/Makefile
	eos/Makefile
	eos/b-decays/Makefile
	eos/c-decays/Makefile
	eos/s-decays/Makefile
	eos/constraints/Makefile
	eos/form-factors/Makefile
	eos/maths/Makefile
	eos/meson-mixing/Makefile
	eos/models/Makefile
	eos/nonlocal-form-factors/Makefile
	eos/nonleptonic-amplitudes/Makefile
	eos/parameters/Makefile
	eos/rare-b-decays/Makefile
	eos/scattering/Makefile
	eos/statistics/Makefile
	eos/utils/Makefile
	src/Makefile
	src/clients/Makefile
	src/scripts/Makefile
	test/Makefile
	python/Makefile
	python/setup.py
	])
AC_OUTPUT
dnl }}}
