#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT([e-antic], [0.1], [vincent.delecroix@math.cnrs.fr])
AM_INIT_AUTOMAKE([-Wall -Werror])
AM_SILENT_RULES([yes]) # non-verbose make; override with ./configure --disable-silent-rules
AM_PROG_AR
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile e-antic/test/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_LIBTOOL
AC_PROG_CXX

# AC_ARG_WITH([flint],
#   AS_HELP_STRING([--with-flint=DIR],
#      [provide location of FLINT (with arb, antic) installation prefix]),
#   [ LDFLAGS="-L$with_flint/lib $LDFLAGS -Wl,-rpath,$with_flint/lib"
#     CPPFLAGS="-I $with_flint/include $CPPFLAGS"
#   ])

AC_CHECK_LIB(arb, arb_init, [LIBS="-larb"],
 [AC_MSG_ERROR([libarb not found])])

AC_CHECK_LIB(flint, fmpz_init, [LIBS="-lflint"],
 [AC_MSG_ERROR([libflint not found])])

AC_CHECK_LIB(antic, nf_init, [LIBS="-lantic"],
 [AC_MSG_ERROR([libantic not found])])

AC_OUTPUT
