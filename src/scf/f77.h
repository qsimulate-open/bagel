#ifndef __src_scf_f77_h
#define __src_scf_f77_h

extern "C" {
 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void dscal_(const int*, const double*, double*, const int*);
 const double ddot_(const int*, const double*, const int*, const double*, const int*); 
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);

 void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*); 

 void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, 
             const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, 
             const double* beta, double* c, const int* ldc);

 void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, 
             double* b, const int* ldb, double* work, const int* lwork, int* info);
};

#endif
