AC_DEFUN([ACX_MPI], [
  AC_ARG_VAR(MPICC,[MPI C compiler command])
  AC_CHECK_PROGS(MPICC, mpicc, openmpicc, $CC)
  acx_mpi_save_CC="$CC"
  CC="$MPICC"
  AC_SUBST(MPICC)

  AC_ARG_VAR(MPICXX,[MPI C++ compiler command])
  AC_CHECK_PROGS(MPICXX, mpicxx mpic++, openmpic++, $CXX)
  acx_mpi_save_CXX="$CXX"
  CXX="$MPICXX"
  AC_SUBST(MPICXX)

  AC_ARG_VAR(MPIF77,[MPI Fortran compiler command])
  AC_CHECK_PROGS(MPIF77, mpif77, openmpif77, $F77)
  acx_mpi_save_F77="$F77"
  F77="$MPIF77"
  AC_SUBST(MPIF77)

  AC_ARG_VAR(MPIFC,[MPI Fortran compiler command])
  AC_CHECK_PROGS(MPIFC, mpif90, openmpif90, $FC)
  acx_mpi_save_FC="$FC"
  FC="$MPIFC"
  AC_SUBST(MPIFC)
])

