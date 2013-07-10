#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([BAGEL], 1.0, [shiozaki@northwestern.edu])
AM_INIT_AUTOMAKE
AM_PATH_PYTHON
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([config.h])
# since pyconfig.h sometimes needs to be compiled by the C compiler
AC_LANG([C])

m4_include([m4/acx_mpi.m4])

AC_PROG_CC([gcc])
AC_PROG_CXX([g++])

AC_ARG_WITH(mpi, [AS_HELP_STRING([--with-mpi],[Use MPI. Please make sure mpicxx and mpif90 are in $PATH])], [with_mpi=$withval], [with_mpi=no])
if test x${with_mpi} = xopenmpi ; then
    echo "..openmpi requested.."
elif test x${with_mpi} = xintel ; then
    echo "..Intel mpi requested.."
elif test x${with_mpi} = xmvapich ; then
    echo "..mvapich requested.."
elif test x${with_mpi} != xno ; then
    AC_MSG_ERROR("Please specify MPI implementation (openmpi, intel, mvapich)")
fi
if test x${with_mpi} != xno; then
    ACX_MPI
    AC_PROG_CC([$MPICC])
    AC_PROG_CXX([$MPICXX])
    LIBS="${MPILIBS} ${LIBS}"
    CXXFLAGS="${CXXFLAGS}"
fi
AX_CXX_COMPILE_STDCXX_11([noext])

LT_INIT([disable-static])
AC_CONFIG_MACRO_DIR([m4])

AC_ARG_WITH([include], [AS_HELP_STRING([--with-include],[include flags])], [with_include=$withval], [with_include=no])
if test "x${with_include}" != xno; then
    CFLAGS="${with_include} ${CFLAGS}"
    CXXFLAGS="${with_include} ${CXXFLAGS}"
    CPPFLAGS="${with_include} ${CPPFLAGS}"
fi

AC_ARG_WITH([slater], [AS_HELP_STRING([--with-slater],[build F12 methods])], [slater=yes], [slater=no])
AC_ARG_WITH([libxc], [AS_HELP_STRING([--with-libxc],[build DFT with libxc])], [libxc=yes], [libxc=no])

#lapack and blas
AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--enable-mkl],[enable MKL extensions.])], [use_mkl=yes], [use_mkl=no])
AC_ARG_ENABLE(scalapack, [AS_HELP_STRING([--enable-scalapack],[enable scalapack - used for parallel diagonalization])], [use_scalapack=yes], [use_scalapack=no])

if test "x${use_scalapack}" = xyes; then
    if test x${with_mpi} != xno; then
        AH_TEMPLATE([HAVE_SCALAPACK], [a scalapack library will be linked.])
        if test x${use_mkl} = xno; then
            AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])
        fi
    else
        AC_MSG_ERROR("scalapack should be used with mpi")
    fi
fi

AX_BOOST_BASE([1.53], [], [AC_MSG_ERROR("BOOST version >= 1.53.0 required.")])

AC_CHECK_LIB(boost_regex-gcc-mt, main, , [
    AC_CHECK_LIB(boost_regex-mt, main, , [
        AC_CHECK_LIB(boost_regex, main, , [
            AC_MSG_ERROR("Linking against boost::regex library failed.")
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
AC_CHECK_LIB(python, main, , [
    AC_CHECK_LIB(python${PYTHON_VERSION}, main, , [
        AC_MSG_ERROR("Linking against python library failed.")
    ])
])
AC_CHECK_LIB(boost_python-gcc-mt, main, , [
    AC_CHECK_LIB(boost_python-mt, main, , [
        AC_CHECK_LIB(boost_python, main, , [
            AC_MSG_ERROR("Linking against boost::python library failed.")
        ])
    ])
])
AC_CHECK_HEADERS([pyconfig.h], [], [AC_MSG_ERROR(["Include path to pyconfig.h required"])], [])

# slatermem.h and srootlist.h need to be compiled with C++11
AC_LANG([C++])
if test x${slater} = xyes; then
  AC_CHECK_HEADERS([slatermem.h], [], [AC_MSG_ERROR(["slatermem.h not found"])], [])
  AC_CHECK_HEADERS([srootlist.h], [], [AC_MSG_ERROR(["srootlist.h not found"])], [])
  AC_CHECK_LIB(slater, main, , AC_MSG_ERROR(["libslater not found"]))
fi

if test x${libxc} = xyes; then
  AC_CHECK_HEADERS([xc.h], [], [AC_MSG_ERROR([xc.h not found or not working])], [])
  AC_CHECK_LIB(xc, main, , AC_MSG_ERROR([libxc not found]))
fi

if test x${use_mkl} = xyes; then
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
    AC_CHECK_LIB([mkl_intel_lp64], main, LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm $LIBS",
        AC_MSG_ERROR([Cannot find Intel Math Kernel Library]), [-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])

    if test "x${use_scalapack}" = xyes; then
        if test x${with_mpi} = xopenmpi; then
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
        else
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
        fi
    fi
fi
if test x${with_mpi} != xno; then
    AC_CHECK_HEADERS([mpi.h], [], [AC_MSG_ERROR([mpi.h not found or not working])], [])
fi

AC_PREFIX_DEFAULT(["/usr/local/bagel"])


AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/util/Makefile
                 src/parallel/Makefile
                 src/dimer/Makefile
                 src/meh/Makefile
                 src/wfn/Makefile
                 src/grad/Makefile
                 src/rysint/Makefile
                 src/osint/Makefile
                 src/scf/Makefile
                 src/df/Makefile
                 src/prop/Makefile
                 src/smith/Makefile
                 src/mp2/Makefile
                 src/math/Makefile
                 src/input/Makefile
                 src/io/Makefile
                 src/fci/Makefile
                 src/zfci/Makefile
                 src/casscf/Makefile
                 src/pscf/Makefile
                 src/pmp2/Makefile
                 src/rel/Makefile
                 src/ks/Makefile
                 src/opt/Makefile
                 src/transp/Makefile])
AC_OUTPUT
