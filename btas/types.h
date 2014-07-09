#ifndef __BTAS_TYPES_H
#define __BTAS_TYPES_H 1

//
//  BLAS types
//

#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

#ifdef _HAS_CBLAS
#ifdef _HAS_INTEL_MKL
#include <mkl_cblas.h>
#else
#include <cblas.h>
#endif // _HAS_INTEL_MKL
#else

/// major order directive
enum CBLAS_ORDER { CblasRowMajor, CblasColMajor };

/// transposition directive in gemv and gemm
enum CBLAS_TRANSPOSE { CblasNoTrans, CblasTrans, CblasConjTrans };

/// upper / lower triangular directive (not used)
enum CBLAS_UPLO { CblasUpper, CblasLower };

/// diagonal matrix directive (not used)
enum CBLAS_DIAG { CblasNonUnit, CblasUnit };

/// transposition directive for symmetric matrix (not used)
enum CBLAS_SIDE { CblasLeft, CblasRight };

#endif // _HAS_CBLAS

#ifdef __cplusplus
}
#endif // __cplusplus

#ifdef _HAS_CBLAS
#include <complex>
#endif // _HAS_CBLAS

namespace btas {

//
//  Other aliases for convenience
//

/// default size type
typedef unsigned long size_type;

/// null deleter
struct nulldeleter { void operator() (void const*) { } };

} // namespace btas

#endif // __BTAS_TYPES_H
