#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([BAGEL], 1.0, [shiozaki@northwestern.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cc])
AM_CONFIG_HEADER([config.h])
AC_F77_LIBRARY_LDFLAGS

# Checks for programs.
AC_PROG_CXX
AC_PROG_FC

AC_CHECK_LIB(boost_regex-gcc-mt, main, , [
    AC_CHECK_LIB(boost_regex-mt, main, , [
        AC_CHECK_LIB(boost_regex, main, , [
            AC_MSG_ERROR("Linking against boost::regex library failed.")
        ])
    ])
])
AC_CHECK_LIB(boost_thread-gcc-mt, main, , [
    AC_CHECK_LIB(boost_thread-mt, main, , [
        AC_CHECK_LIB(boost_thread, main, , [
            AC_MSG_ERROR("Linking against boost::thread library failed.")
        ])
    ])
])
AC_CHECK_LIB(boost_filesystem-gcc-mt, main, , [
    AC_CHECK_LIB(boost_filesystem-mt, main, , [
        AC_CHECK_LIB(boost_filesystem, main, , [
            AC_MSG_ERROR("Linking against boost::filesystem library failed.")
        ])
    ])
])
AC_CHECK_LIB(boost_system-gcc-mt, main, , [
    AC_CHECK_LIB(boost_system-mt, main, , [
        AC_CHECK_LIB(boost_system, main, , [
            AC_MSG_ERROR("Linking against boost::system library failed.")
        ])
    ])
])
AC_CHECK_LIB(boost_unit_test_framework-gcc-mt, main, , [
    AC_CHECK_LIB(boost_unit_test_framework-mt, main, , [
        AC_CHECK_LIB(boost_unit_test_framework, main, , [
            AC_MSG_ERROR("Linking against boost::unit_test_framework library failed.")
        ])
    ])
])

LT_INIT
AC_CONFIG_MACRO_DIR([m4])


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/parallel/Makefile
                 src/dimer/Makefile
                 src/opt/Makefile
                 src/wfn/Makefile
                 src/grad/Makefile
                 src/rysint/Makefile
                 src/osint/Makefile
                 src/slater/Makefile
                 src/scf/Makefile
                 src/df/Makefile
                 src/prop/Makefile
                 src/smith/Makefile
                 src/mp2/Makefile
                 src/util/Makefile
                 src/io/Makefile
                 src/fci/Makefile
                 src/casscf/Makefile
                 src/pscf/Makefile
                 src/pmp2/Makefile
                 src/rel/Makefile])
AC_OUTPUT
