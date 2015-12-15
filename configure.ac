AC_INIT([trng], [4.19], [trng@mail.de])
AM_INIT_AUTOMAKE()

AC_DEFINE([TRNG_VERSION], [4.19], [TRNG version number.])

AC_LANG_CPLUSPLUS
AC_LANG(C++)

# Checks for programs.
AC_PROG_CXX
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PROG_LIBTOOL

# Check for libraries and functions
AC_CXX_LONG_LONG(,)
AC_PATH_BOOST(103000, 
              AC_DEFINE(TRNG_HAVE_BOOST, 1, [Define to 1 if you have BOOST library.]),
	      )
AC_SUBST(BOOST_CXXFLAGS)
AC_LANG_CPLUSPLUS
ACX_MPI(AC_DEFINE(TRNG_HAVE_MPI, 1, [Define to 1 if you have a MPI library.]))
AX_OPENMP
AC_SUBST(OPENMP_CXXFLAGS)
AC_LANG_PUSH([C++])
AC_CHECK_FUNCS(erff, AC_DEFINE(TRNG_HAVE_ERFF, 1, [Define to 1 if you have the `errff' function.]))
#AC_TRY_COMPILE([#include <cmath>
#extern "C" float erff(float) throw(); ], [erff(1.f); ], AC_MSG_RESULT(ok), AC_MSG_RESULT(fail) )
AC_CHECK_FUNCS(erf, AC_DEFINE(TRNG_HAVE_ERF, 1, [Define to 1 if you have the `errf' function.]))
AC_CHECK_FUNCS(erfl, AC_DEFINE(TRNG_HAVE_ERFL, 1, [Define to 1 if you have the `errfl' function.]))
AC_CHECK_FUNCS(erfcf, AC_DEFINE(TRNG_HAVE_ERFCF, 1, [Define to 1 if you have the `erfcf' function.]))
AC_CHECK_FUNCS(erfc, AC_DEFINE(TRNG_HAVE_ERFC, 1, [Define to 1 if you have the `erfc' function.]))
AC_CHECK_FUNCS(erfcl, AC_DEFINE(TRNG_HAVE_ERFCL, 1, [Define to 1 if you have the `erfcl' function.]))
AC_CHECK_FUNCS(tgammaf, AC_DEFINE(TRNG_HAVE_TGAMMAF, 1, [Define to 1 if you have the `tgammaf' function.]))
AC_CHECK_FUNCS(tgamma, AC_DEFINE(TRNG_HAVE_TGAMMA, 1, [Define to 1 if you have the `tgamma' function.]))
AC_CHECK_FUNCS(tgammal, AC_DEFINE(TRNG_HAVE_TGAMMAL, 1, [Define to 1 if you have the `tgammal' function.]))
AC_CHECK_FUNCS(lgammaf, AC_DEFINE(TRNG_HAVE_LGAMMAF, 1, [Define to 1 if you have the `lgammaf' function.]))
AC_CHECK_FUNCS(lgamma, AC_DEFINE(TRNG_HAVE_LGAMMA, 1, [Define to 1 if you have the `lgamma' function.]))
AC_CHECK_FUNCS(lgammal, AC_DEFINE(TRNG_HAVE_LGAMMAL, 1, [Define to 1 if you have the `lgammal' function.]))
AC_LANG_POP
# check for  Intel Threading Building Blocks
AC_CHECK_HEADER([tbb/parallel_for.h],	
                [ LIB_TBB="-ltbb" 
                  AC_SUBST([LIB_TBB])
                  AC_DEFINE([TRNG_HAVE_TBB], [1], [Define to 1 if Intel Threading Building Blocks are available.])
                  AC_MSG_RESULT(Will generate parallel code using Intel Threading Building Blocks.)], 
                  AC_MSG_RESULT(Cannot find Intel Threading Building Blocks.))

AM_CONFIG_HEADER(src/config.hpp)

AC_OUTPUT(Makefile trng/Makefile src/Makefile examples/Makefile)
