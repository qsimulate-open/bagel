//
// Author : Toru Shiozaki
// Date   : Dec 2011
// Blas & Lapack interface 
//

#ifndef __src_util_f77_h
#define __src_util_f77_h

#include <complex>
#include <memory>

extern "C" {

 // transposition
 void mytranspose_(const double*, const int*, const int*, double*);
 void mytranspose1_(const double*, const int*, const int*, double*);
 void mytranspose4_(const double*, const int*, const int*, double*);
 void mytranspose_complex_(const std::complex<double>*, const int*, const int*, std::complex<double>*);

 void dcopy_(const int*, const double*, const int*, double*, const int*);
 void daxpy_(const int*, const double*, const double*, const int*, double*, const int*);
 void dsyev_(const char*, const char*, const int*, double*, const int*, double*, double*, const int*, int*); 
 void dgesv_(const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, double* b, const int* ldb, int* info);
 void dscal_(const int*, const double*, double*, const int*);
 double ddot_(const int*, const double*, const int*, const double*, const int*); 
 void dgemv_(const char*, const int*, const int*, const double*, const double*, const int*, const double*, const int*,
             const double*, double*, const int*);
 void dgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, 
             const double* alpha, const double* a, const int* lda, const double* b, const int* ldb, 
             const double* beta, double* c, const int* ldc);
 void dsysv_(const char* uplo, const int* n, const int* nrhs, double* a, const int* lda, int* ipiv, 
             double* b, const int* ldb, double* work, const int* lwork, int* info);

 void zcopy_(const int*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
 void zscal_(const int*, const std::complex<double>*, std::complex<double>*, const int*);
#ifndef ZDOT_RETURN 
 void zdotu_(std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*); 
 void zdotc_(std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*); 
#else
 std::complex<double> zdotu_(const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*); 
 std::complex<double> zdotc_(const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*); 
#endif
 void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
 void zheev_(const char*, const char*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*, const int*, double*, int*); 
 void zgeev_(const char*, const char*, const int*, std::complex<double>*, const int*, std::complex<double>*, 
             std::complex<double>*, const int*, std::complex<double>*, const int*, std::complex<double>*, const int*, double*, int*);
 void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k, 
             const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb, 
             const std::complex<double>* beta, std::complex<double>* c, const int* ldc);
 void zhesv_(const char* uplo, const int* n, const int* nrhs, std::complex<double>* a, const int* lda, int* ipiv, 
             std::complex<double>* b, const int* ldb, std::complex<double>* work, const int* lwork, int* info);
 void zgesv_(const int* n, const int* nrhs, std::complex<double>* a, const int* lda, int* ipiv, 
             std::complex<double>* b, const int* ldb, int* info);
 void zgesvd_(const char*, const char*, const int*, const int*, std::complex<double>*, const int*, const double*,
     std::complex<double>*, const int*, std::complex<double>*, const int*,  std::complex<double>*, const int*, double*, int*);
 void zgesdd_(const char*, const int*, const int*, std::complex<double>*, const int*, const double*,
     std::complex<double>*, const int*, std::complex<double>*, const int*,  std::complex<double>*, const int*, double*, int*, int*);

};


static void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k, 
                   const double alpha, const double* a, const int lda, const double* b, const int ldb, 
                   const double beta, double* c, const int ldc) { dgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }; 
static void dgemm_(const char* transa, const char* transb, const int m, const int n, const int k, 
                   const double alpha, const std::unique_ptr<double []>& a, const int lda, const std::unique_ptr<double []>& b, const int ldb, 
                   const double beta, std::unique_ptr<double []>& c, const int ldc)
                   { dgemm_(transa,transb,&m,&n,&k,&alpha,a.get(),&lda,b.get(),&ldb,&beta,c.get(),&ldc); }; 
static void dgemv_(const char* a, const int b, const int c, const double d, const double* e, const int f, const double* g, const int h,
                   const double i, double* j, const int k) { dgemv_(a,&b,&c,&d,e,&f,g,&h,&i,j,&k); };
static void dgemv_(const char* a, const int b, const int c, const double d, const std::unique_ptr<double []>& e, const int f,
                   const std::unique_ptr<double []>& g, const int h, const double i, std::unique_ptr<double []>& j, const int k)
                   { dgemv_(a,&b,&c,&d,e.get(),&f,g.get(),&h,&i,j.get(),&k); };
static void daxpy_(const int a, const double b, const double* c, const int d, double* e, const int f) {daxpy_(&a,&b,c,&d,e,&f); };
static void daxpy_(const int a, const double b, const std::unique_ptr<double []>& c, const int d, std::unique_ptr<double []>& e, const int f)
                   {daxpy_(&a,&b,c.get(),&d,e.get(),&f); };
static void dcopy_(const int a, const double* b, const int c, double* d, const int e) {dcopy_(&a, b, &c, d, &e); };
static void dcopy_(const int a, const std::unique_ptr<double []>& b, const int c, std::unique_ptr<double []>& d, const int e)
                  {dcopy_(&a, b.get(), &c, d.get(), &e); };
static void dscal_(const int a, const double b, double* c, const int d) {dscal_(&a, &b, c, &d); };
static void dscal_(const int a, const double b, std::unique_ptr<double []>& c, const int d) {dscal_(&a, &b, c.get(), &d); };
static double ddot_(const int a, const double* b, const int c, const double* d, const int e) { return ddot_(&a,b,&c,d,&e); };
static double ddot_(const int a, const std::unique_ptr<double []>& b, const int c, const std::unique_ptr<double []>& d, const int e)
                   { return ddot_(&a,b.get(),&c,d.get(),&e); };
static double dsyev_(const char* a, const char* b, const int c, double* d, const int e, double* f, double* g, const int h, int& i)
                    {dsyev_(a,b,&c,d,&e,f,g,&h,&i);};
static double dsyev_(const char* a, const char* b, const int c, std::unique_ptr<double []>& d, const int e,
                     std::unique_ptr<double []>& f, std::unique_ptr<double []>& g, const int h, int& i)
                    {dsyev_(a,b,&c,d.get(),&e,f.get(),g.get(),&h,&i);};

#endif
