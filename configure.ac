#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_INIT(e-antic, 0.1)
AM_INIT_AUTOMAKE(foreign)
AM_SILENT_RULES([yes]) # non-verbose make; override with ./configure --disable-silent-rules
AC_CONFIG_FILES([Makefile e-antic/test/Makefile])
AC_CONFIG_MACRO_DIR([m4])
AC_PROG_LIBTOOL
AC_PROG_CXX

## Test for OpenMP parallelization.
AC_ARG_ENABLE([openmp],
    [AS_HELP_STRING([--enable-openmp@<:@=ARG@:>@],
      [enable parallelization with OpenMP @<:@default=check@:>@])],
    [:],
    [enable_openmp=check])
AS_IF([test "x$enable_openmp" != xno],
  [AX_OPENMP([have_openmp=yes],
	     [have_openmp=no])])

AS_IF([test x$have_openmp = xyes],
  [enable_openmp=yes],
  [AS_IF([test "x$enable_openmp" = xyes],
    [AC_MSG_ERROR([openmp not available])],
    [AS_IF([test "x$enable_openmp" != xno],
        [AC_MSG_NOTICE([Parallelization with OpenMP will not be available because no suitable compiler was found])])
     enable_openmp=no])])

AC_SUBST(OPENMP_CXXFLAGS)
AM_CONDITIONAL(ENABLE_OPENMP, [test x$enable_openmp = xyes])

AC_ARG_WITH([flint],
   AS_HELP_STRING([--with-flint=DIR],
      [provide location of FLINT (with arb, antic) installation prefix]),
   [ LDFLAGS="-L$with_flint/lib $LDFLAGS -Wl,-rpath,$with_flint/lib"
     CPPFLAGS="-I $with_flint/include $CPPFLAGS"
   ])

AC_OUTPUT
