#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.69)
AC_INIT([BAGEL], 1.2.0, [shiozaki@northwestern.edu], [], [http://www.nubakery.org])
AM_INIT_AUTOMAKE
AM_PATH_PYTHON
AC_CONFIG_SRCDIR([src/main.cc])
AC_CONFIG_HEADERS([bagel_config.h])
# since pyconfig.h sometimes needs to be compiled by the C compiler
AC_LANG([C])

AC_PREFIX_DEFAULT(["/usr/local/bagel"])

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
AC_ARG_ENABLE(mkl, [AS_HELP_STRING([--enable-mkl],[enable MKL extensions.])], [use_mkl=$enable_mkl], [use_mkl=no])
AC_ARG_ENABLE(scalapack, [AS_HELP_STRING([--disable-scalapack],[disable scalapack.])], [use_scalapack=$enable_scalapack], [use_scalapack=yes])

if test "x${use_scalapack}" = xyes; then
    if test "x${with_mpi}" != xno; then
        AH_TEMPLATE([HAVE_SCALAPACK], [a scalapack library will be linked.])
        if test "x${use_mkl}" = xno; then
            AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])
        fi
    fi
fi

AX_BOOST_BASE([1.53], [], [AC_MSG_ERROR("BOOST version >= 1.53.0 required.")])
AX_BOOST_REGEX
AX_BOOST_UNIT_TEST_FRAMEWORK
AX_BOOST_SERIALIZATION
CPPFLAGS="${CPPFLAGS} ${BOOST_CPPFLAGS}"
LDFLAGS="${BOOST_LDFLAGS} ${LDFLAGS}"
LIBS="${BOOST_REGEX_LIB} ${BOOST_UNIT_TEST_FRAMEWORK_LIB} ${BOOST_SERIALIZATION_LIB} ${LIBS}"

AC_LANG_PUSH([C++])

# srootlist.h need to be compiled with C++11
if test "x${slater}" = xyes; then
  AC_CHECK_HEADERS([srootlist.h], [], [AC_MSG_ERROR(["srootlist.h not found"])], [])
  AC_CHECK_LIB(slater, main, , AC_MSG_ERROR(["libslater not found"]))
fi

if test "x${libxc}" = xyes; then
  AC_CHECK_HEADERS([xc.h], [], [AC_MSG_ERROR([xc.h not found or not working])], [])
  AC_CHECK_LIB(xc, main, , AC_MSG_ERROR([libxc not found]))
fi

#toggle for smith
AC_ARG_ENABLE(smith, [AS_HELP_STRING([--disable-smith],[disable SMITH code.])], [compile_smith=$enable_smith], [compile_smith=yes])
if test "x${compile_smith}" != xno; then
  if test "x${with_mpi}" = xno; then
    AC_MSG_ERROR("SMITH requires MPI")
  fi
  AC_DEFINE([COMPILE_SMITH], [], [Compile SMITH generated code])
fi

AC_LANG_POP()

if test "x${use_mkl}" = xyes; then
    AC_CHECK_HEADERS([mkl.h], [], [AC_MSG_ERROR([mkl.h not found or not working])], [])
    AC_CHECK_LIB([mkl_intel_lp64], main, LIBS="-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm $LIBS",
        AC_MSG_ERROR([Cannot find Intel Math Kernel Library]), [-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])

    if test "x${use_scalapack}" = xyes; then
        if test "x${with_mpi}" = xopenmpi; then
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
        elif test "x${with_mpi}" != xno; then
            AC_CHECK_LIB([mkl_scalapack_lp64], main, AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 $LIBS",
                [AC_CHECK_LIB(scalapack, main,  [AC_DEFINE([HAVE_SCALAPACK]) LIBS="-lscalapack $LIBS"], [AC_MSG_ERROR("Linking against scalapack library failed.")])],
                [-lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lm])
        fi
    fi
fi

# ZGEMM3M is provided by MKL, OpenBLAS and perhaps others.
AC_CHECK_FUNC([zgemm3m_], [AC_DEFINE([HAVE_ZGEMM3M], [1], [Define if zgemm3m_ is present in BLAS .])])

if test "x${with_mpi}" != xno; then
    AC_CHECK_HEADERS([mpi.h], [], [AC_MSG_ERROR([mpi.h not found or not working])], [])
fi

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/util/Makefile
                 src/util/input/Makefile
                 src/util/io/Makefile
                 src/util/parallel/Makefile
                 src/util/math/Makefile
                 src/asd/multisite/Makefile
                 src/asd/dimer/Makefile
                 src/asd/dmrg/Makefile
                 src/asd/orbital/Makefile
                 src/asd/Makefile
                 src/wfn/Makefile
                 src/grad/Makefile
                 src/integral/Makefile
                 src/scf/Makefile
                 src/df/Makefile
                 src/prop/Makefile
                 src/smith/Makefile
                 src/molecule/Makefile
                 src/pt2/Makefile
                 src/ci/Makefile
                 src/ci/ciutil/Makefile
                 src/ci/fci/Makefile
                 src/ci/zfci/Makefile
                 src/ci/ras/Makefile
                 src/response/Makefile
                 src/multi/Makefile
                 src/multi/casscf/Makefile
                 src/multi/zcasscf/Makefile
                 src/mat1e/Makefile
                 src/opt/Makefile
                 src/basis/Makefile])
AC_OUTPUT
