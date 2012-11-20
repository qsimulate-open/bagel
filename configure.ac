#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.57)
AC_INIT([BAGEL], 1.0, [shiozaki@northwestern.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/main.cc])
AM_CONFIG_HEADER([config.h])

AC_PROG_CC([gcc])
AC_PROG_CXX([g++])
AC_PROG_F77([gfortran])
AC_PROG_FC([gfortran])

AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--enable-mkl],[enable MKL extensions.])], [use_mkl=yes], [use_mkl=no])
if test x${use_mkl} = xyes; then
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])])
    AC_CHECK_LIB(mkl_intel_lp64, main)
    AC_CHECK_LIB(mkl_intel_thread, main)
    AC_CHECK_LIB(mkl_core, main)
    AC_CHECK_LIB(iomp5, main)
    AC_CHECK_LIB(pthread, main)
    AC_CHECK_LIB(m, main)
fi

AC_ARG_ENABLE(mpi, [AS_HELP_STRING([--enable-mpi],[enable MPI. Please make sure mpicxx and mpif90 are in $PATH [default=no]])], [use_mpi=yes], [use_mpi=no])
if test x${use_mpi} = xyes; then
    AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
    AC_CHECK_PROGS(MPICXX, mpic++ mpiCC mpicxx mpCC hcp mpxlC mpxlC_r cmpic++, [AC_MSG_ERROR([mpic++ not found or not working])])

    AC_ARG_VAR(MPIF77,[MPI Fortran 77 compiler command])
    AC_CHECK_PROGS(MPIF77, mpif77 hf77 mpxlf mpf77 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, [AC_MSG_ERROR([mpif77 not found or not working])])

    AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
    AC_CHECK_PROGS(MPIFC, mpif90 hf90 mpxlf mpf90 mpif90 mpf90 mpxlf90 mpxlf95 mpxlf_r cmpifc cmpif90c, [AC_MSG_ERROR([mpif90 not found or not working])])

    CXX="${MPICXX}"
    F77="${MPIF77}"
    FC="${MPIFC}"
    LIBS="${MPILIBS} ${LIBS}"
    CXXFLAGS="${CXXFLAGS}"
    FFLAGS="${FFLAGS}"
    FCFLAGS="${FCFLAGS}"
fi

AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
if test "x${with_include}" != xno; then
    CFLAGS="${with_include} ${CFLAGS}"
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    FFLAGS="${with_include} ${FFLAGS}"
    FCFLAGS="${with_include} ${FCFLAGS}"
fi

if test x${use_mpi} = xyes; then
    AC_CHECK_HEADERS([mpi.h], [], [AC_MSG_ERROR([mpi.h not found or not working])])
    AC_CHECK_LIB(mpi, MPI_Init, [MPILIBS="-lmpi"])
    AC_CHECK_LIB(mpich, MPI_Init, [MPILIBS="-lmpich"])
fi

AC_F77_LIBRARY_LDFLAGS

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

LT_INIT([disable-static])
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
