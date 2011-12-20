#ifndef __src_fci_f77_h
#define __src_fci_f77_h

extern "C" {
 void mytranspose_(const double*, const int*, const int*, double*);
 const double ddot_(const int*, const double*, const int*, const double*, const int*); 
 void dscal_(const int*, const double*, double*, const int*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, 
             const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, 
             const double* beta, double* c, const int* ldc);
};

#endif
