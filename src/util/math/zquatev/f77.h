//
// ZQUATEV: Diagonalization of quaternionic matrices
// File   : f77.h
// Copyright (c) 2013, Toru Shiozaki (shiozaki@northwestern.edu)
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// The views and conclusions contained in the software and documentation are those
// of the authors and should not be interpreted as representing official policies,
// either expressed or implied, of the FreeBSD Project.
//

#ifndef __TS_ZQUATEV_F77_H
#define __TS_ZQUATEV_F77_H

#include <bagel_config.h>
#include <complex>

// blas
extern "C" {

 void zscal_(const int*, const std::complex<double>*, std::complex<double>*, const int*);
 void zdotc_(std::complex<double>*, const int*, const std::complex<double>*, const int*, const std::complex<double>*, const int*);
 void zaxpy_(const int*, const std::complex<double>*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
 void zgemv_(const char*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*,
             const std::complex<double>*, std::complex<double>*, const int*);
 void ztrmv_(const char*, const char*, const char*, const int*, const std::complex<double>*, const int*, std::complex<double>*, const int*);
#ifdef HAVE_ZGEMM3M
 void zgemm3m_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
               const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
               const std::complex<double>* beta, std::complex<double>* c, const int* ldc);
#else
 void zgemm_(const char* transa, const char* transb, const int* m, const int* n, const int* k,
             const std::complex<double>* alpha, const std::complex<double>* a, const int* lda, const std::complex<double>* b, const int* ldb,
             const std::complex<double>* beta, std::complex<double>* c, const int* ldc);
#endif
 void zrot_(const int*, std::complex<double>*, const int*, std::complex<double>*, const int*, const double*, const std::complex<double>*);
 void zgerc_(const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*,
             std::complex<double>*, const int*);
 void zgeru_(const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const std::complex<double>*, const int*,
             std::complex<double>*, const int*);

 // lapack
 void zheev_(const char*, const char*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*, const int*, double*, int*);
 void zhbev_(const char*, const char*, const int*, const int*, std::complex<double>*, const int*, double*, std::complex<double>*, const int*,
             std::complex<double>*, double*, int*);
 void zlartg_(const std::complex<double>*, const std::complex<double>*, double*, std::complex<double>*, std::complex<double>*);
 void zlarfg_(const int*, const std::complex<double>*, std::complex<double>*, const int*, std::complex<double>*);
}


namespace {

 void zgemv_(const char* a, const int b, const int c, const std::complex<double>& d, const std::complex<double>* e, const int f, const std::complex<double>* g, const int h,
             const std::complex<double>& i, std::complex<double>* j, const int k) { ::zgemv_(a,&b,&c,&d,e,&f,g,&h,&i,j,&k); }
 void ztrmv_(const char* a, const char* b, const char* c, const int d, const std::complex<double>* e, const int f, std::complex<double>* g, const int h)
            { ::ztrmv_(a,b,c,&d,e,&f,g,&h); }
#ifdef HAVE_ZGEMM3M
 void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
               const std::complex<double>& alpha, const std::complex<double>* a, const int lda, const std::complex<double>* b, const int ldb,
               const std::complex<double>& beta, std::complex<double>* c, const int ldc) { ::zgemm3m_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }
#else
 void zgemm3m_(const char* transa, const char* transb, const int m, const int n, const int k,
             const std::complex<double> alpha, const std::complex<double>* a, const int lda, const std::complex<double>* b, const int ldb,
             const std::complex<double> beta, std::complex<double>* c, const int ldc) { ::zgemm_(transa,transb,&m,&n,&k,&alpha,a,&lda,b,&ldb,&beta,c,&ldc); }
#endif

 void zaxpy_(const int a, const std::complex<double>& b, const std::complex<double>* c, const int d, std::complex<double>* e, const int f) { ::zaxpy_(&a,&b,c,&d,e,&f); }
 void zscal_(const int a, const std::complex<double>& b, std::complex<double>* c, const int d) { ::zscal_(&a, &b, c, &d); }
 std::complex<double> zdotc_(const int b, const std::complex<double>* c, const int d, const std::complex<double>* e, const int f) {
   std::complex<double> a;
   ::zdotc_(&a,&b,c,&d,e,&f);
   return a;
 }
 void zheev_(const char* a, const char* b, const int c, std::complex<double>* d, const int e, double* f, std::complex<double>* g, const int h, double* i, int& j)
             { ::zheev_(a,b,&c,d,&e,f,g,&h,i,&j); }
 void zhbev_(const char* a, const char* b, const int c, const int d, std::complex<double>* e, const int f, double* g, std::complex<double>* h, const int i,
             std::complex<double>* j, double* k, int& l) { ::zhbev_(a, b, &c, &d, e, &f, g, h, &i, j, k, &l); }
 void zrot_(const int a, std::complex<double>* b, const int c, std::complex<double>* d, const int e, const double f, const std::complex<double>& g) {
            ::zrot_(&a, b, &c, d, &e, &f, &g); }
 void zgerc_(const int a, const int b, const std::complex<double>& c, const std::complex<double>* d, const int e, const std::complex<double>* f, const int g,
             std::complex<double>* h, const int i) { ::zgerc_(&a, &b, &c, d, &e, f, &g, h, &i); }
 void zgeru_(const int a, const int b, const std::complex<double>& c, const std::complex<double>* d, const int e, const std::complex<double>* f, const int g,
             std::complex<double>* h, const int i) { ::zgeru_(&a, &b, &c, d, &e, f, &g, h, &i); }
 void zlartg_(const std::complex<double>& a, const std::complex<double>& b, double& c, std::complex<double>& d, std::complex<double>& e) { ::zlartg_(&a, &b, &c, &d, &e); }
 void zlarfg_(const int a, std::complex<double>& b, std::complex<double>* c, const int d, std::complex<double>& e) { ::zlarfg_(&a, &b, c, &d, &e); }

}

#endif
