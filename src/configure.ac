#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([TelSeq], [0.0.3], [tnakao@broadinstitute.org])
AM_INIT_AUTOMAKE(foreign)
AC_CONFIG_SRCDIR([Telseq/telseq.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX
AC_PROG_RANLIB

# Checks for header files.
AC_LANG([C++])

# Set compiler flags.
AC_SUBST(AM_CXXFLAGS, "-Wall -Wextra -Wno-unknown-pragmas -g -O3 -pg")

AC_SUBST(CXXFLAGS, "-std=c++11")
AC_SUBST(CFLAGS, "-std=c++11")

AC_SUBST(LIBS, "$LIBS -lz -lpthread")

# use huslib instead of bamtools
AC_ARG_WITH([htslib],
  [AS_HELP_STRING([--with-htslib=PATH], [Path to htslib directory])],
  [htslib_prefix="$withval"],
  [htslib_prefix="/deps/htslib-1.18"])

AC_MSG_CHECKING([for htslib library and headers])
PKG_CHECK_MODULES([HTSLIB], [htslib >= 1.18], [], [AC_MSG_ERROR([htslib >= 1.18 not found.])])
AC_MSG_RESULT([found])
AC_SUBST([HTSLIB_CFLAGS])
AC_SUBST([HTSLIB_LIBS], [-L${htslib_prefix}/lib -lhts])

AC_CONFIG_FILES([Makefile
		Util/Makefile
    Telseq/Makefile])

AC_OUTPUT