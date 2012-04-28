#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([PolyR12], 1.0, [shiozaki.toru@gmail.com])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cc])
AM_CONFIG_HEADER([config.h])
AC_F77_LIBRARY_LDFLAGS

# Checks for programs.
AC_PROG_CXX
#AC_PROG_LIBTOOL
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_FC

AC_CONFIG_MACRO_DIR([m4])

# Checks for libraries.

# Checks for header files.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/wfn/Makefile
                 src/rysint/Makefile
                 src/osint/Makefile
                 src/slater/Makefile
                 src/scf/Makefile
                 src/df/Makefile
                 src/smith/Makefile
                 src/mp2/Makefile
                 src/pmp2/Makefile
                 src/util/Makefile
                 src/fci/Makefile
                 src/casscf/Makefile
                 src/pscf/Makefile])
AC_OUTPUT
