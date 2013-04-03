dnl AC_PATH_BOOST([MINIMUM-VERSION], [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
dnl Test for the Boost C++ libraries of a particular version (or newer)
dnl Defines:
dnl   BOOST_CXXFLAGS to the set of flags required to compiled Boost
AC_DEFUN([AC_PATH_BOOST], 
[
  BOOST_CXXFLAGS=""
  path_given="no"

dnl Extract the path name from a --with-boost=PATH argument
  AC_ARG_WITH(boost,
    [  --with-boost=PATH absolute path name where the Boost C++ libraries
    reside. Alternatively, the BOOST_ROOT environment variable will be used],
    if test "$withval" = no ; then
	path_given="no"
	BOOST_CXXFLAGS=""
    else
      if test "$withval" != yes ; then
        path_given="yes"
        BOOST_CXXFLAGS="-I$withval"
        BOOST_ROOT="$withval"
      fi
    fi    
  )

dnl If no path with given and there is a BOOST_ROOT environment variable,
dnl use it
  if test "$path_given" = "no"; then
    if test "x$BOOST_ROOT" = "x"; then
      BOOST_CXXFLAGS=""
    else
      BOOST_CXXFLAGS="-I$BOOST_ROOT"
    fi
  fi

  boost_min_version=ifelse([$1], ,102000,$1)

  AC_LANG_SAVE
  AC_LANG_CPLUSPLUS
  OLD_CXXFLAGS=$CXXFLAGS
  CXXFLAGS="$CXXFLAGS $BOOST_CXXFLAGS"
  AC_MSG_CHECKING([for the Boost C++ libraries, version $boost_min_version or newer])
  AC_TRY_COMPILE(
    [
#include <boost/version.hpp>
],
    [],
    [
      have_boost="yes"
    ],
    [
      AC_MSG_RESULT(no)
      have_boost="no"
      ifelse([$3], , :, [$3])
    ])

  if test "$have_boost" = "yes"; then
    WANT_BOOST_VERSION=$boost_min_version
    AC_TRY_COMPILE(
      [
#include <boost/version.hpp>
],
      [
#if BOOST_VERSION >= $WANT_BOOST_VERSION
// Everything is okay
#else
#  error Boost version is too old
#endif

],
      [
        AC_MSG_RESULT(yes)
        ifelse([$2], , :, [$2])
      ],
      [
        AC_MSG_RESULT([no, version of installed Boost libraries is too old])
        ifelse([$3], , :, [$3])
      ])
  fi
  CXXFLAGS=$OLD_CXXFLAGS
  AC_LANG_RESTORE
])


AC_DEFUN([AC_CXX_LONG_LONG],
[AC_CACHE_CHECK(whether the compiler recognizes long long as a built-in type,
ac_cv_cxx_long_long,
[AC_LANG_SAVE
 AC_LANG_CPLUSPLUS
 AC_TRY_COMPILE([
int f(int  x){return 1;}
int f(char x){return 1;}
int f(long long x){return 1;}
],[long long b = true; return f(b);],
 ac_cv_cxx_long_long=yes, ac_cv_cxx_long_long=no)
 AC_LANG_RESTORE
])
if test "$ac_cv_cxx_long_long" = yes; then
   AC_DEFINE(TRNG_HAVE_LONG_LONG,,[Define to 1 if long long is a built-in type.])
 else
   echo ""
   echo "Sorry, your compiler does not support the built-in type long long."
   exit 1
fi
])


AC_DEFUN([ACX_MPI], [
AC_PREREQ(2.50) dnl for AC_LANG_CASE

AC_LANG_CASE([C], [
        AC_REQUIRE([AC_PROG_CC])
        AC_ARG_VAR(MPICC,[MPI C compiler command])
        AC_CHECK_PROGS(MPICC, mpicc hcc mpcc mpcc_r mpxlc cmpicc, $CC)
        acx_mpi_save_CC="$CC"
        CC="$MPICC"
        AC_SUBST(MPICC)
],
[C++], [
        AC_REQUIRE([AC_PROG_CXX])
        AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
        AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, $CXX)
        acx_mpi_save_CXX="$CXX"
        CXX="$MPICXX"
        AC_SUBST(MPICXX)
],
[Fortran 77], [
        AC_REQUIRE([AC_PROG_F77])
        AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
        AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, $F77)
        acx_mpi_save_F77="$F77"
        F77="$MPIF77"
        AC_SUBST(MPIF77)
])

if test x = x"$MPILIBS"; then
        AC_LANG_CASE([C], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [C++], [AC_CHECK_FUNC(MPI_Init, [MPILIBS=" "])],
                [Fortran 77], [AC_MSG_CHECKING([for MPI_Init])
                        AC_TRY_LINK([],[      call MPI_Init], [MPILIBS=" "
                                AC_MSG_RESULT(yes)], [AC_MSG_RESULT(no)])])
fi
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
fi
if test x = x"$MPILIBS"; then
        AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

dnl We have to use AC_TRY_COMPILE and not AC_CHECK_HEADER because the
dnl latter uses $CPP, not $CC (which may be mpicc).
AC_LANG_CASE([C], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi],
[C++], [if test x != x"$MPILIBS"; then
        AC_MSG_CHECKING([for mpi.h])
        AC_TRY_COMPILE([#include <mpi.h>],[],[AC_MSG_RESULT(yes)], [MPILIBS=""
                AC_MSG_RESULT(no)])
fi])

AC_LANG_CASE([C], [CC="$acx_mpi_save_CC"],
        [C++], [CXX="$acx_mpi_save_CXX"],
        [Fortran 77], [F77="$acx_mpi_save_F77"])

AC_SUBST(MPILIBS)

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x = x"$MPILIBS"; then
        $2
        :
else
        ifelse([$1],,[AC_DEFINE(HAVE_MPI,1,[Define if you have the MPI library.])],[$1])
        :
fi
])dnl ACX_MPI







AC_DEFUN([AX_OPENMP], [
AC_PREREQ(2.59) dnl for _AC_LANG_PREFIX

AC_CACHE_CHECK([for OpenMP flag of _AC_LANG compiler], ax_cv_[]_AC_LANG_ABBREV[]_openmp, [save[]_AC_LANG_PREFIX[]FLAGS=$[]_AC_LANG_PREFIX[]FLAGS
ax_cv_[]_AC_LANG_ABBREV[]_openmp=unknown
# Flags to try:  -fopenmp (gcc), -openmp (icc), -mp (SGI & PGI),
#                -xopenmp (Sun), -omp (Tru64), -qsmp=omp (AIX), none
ax_openmp_flags="-fopenmp -openmp -mp -xopenmp -omp -qsmp=omp none"
if test "x$OPENMP_[]_AC_LANG_PREFIX[]FLAGS" != x; then
  ax_openmp_flags="$OPENMP_[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flags"
fi
for ax_openmp_flag in $ax_openmp_flags; do
  case $ax_openmp_flag in
    none) []_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[] ;;
    *) []_AC_LANG_PREFIX[]FLAGS="$save[]_AC_LANG_PREFIX[]FLAGS $ax_openmp_flag" ;;
  esac
  AC_TRY_LINK([#ifdef __cplusplus
extern "C"
#endif
void omp_set_num_threads(int);], [const int N = 100000;
  int i, arr[N];

  omp_set_num_threads(2);

  #pragma omp parallel for
  for (i = 0; i < N; i++) {
    arr[i] = i;
  }], [ax_cv_[]_AC_LANG_ABBREV[]_openmp=$ax_openmp_flag; break])
done
[]_AC_LANG_PREFIX[]FLAGS=$save[]_AC_LANG_PREFIX[]FLAGS
])
if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" = "xunknown"; then
  m4_default([$2],:)
else
  if test "x$ax_cv_[]_AC_LANG_ABBREV[]_openmp" != "xnone"; then
    OPENMP_[]_AC_LANG_PREFIX[]FLAGS=$ax_cv_[]_AC_LANG_ABBREV[]_openmp
  fi
  m4_default([$1], [AC_DEFINE(TRNG_HAVE_OPENMP,1,[Define if OpenMP is enabled])])
fi
])dnl AX_OPENMP
