//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: scalapack.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#ifndef __SRC_PARALLEL_SCALAPACK_H
#define __SRC_PARALLEL_SCALAPACK_H

#include <bagel_config.h>
#ifdef HAVE_SCALAPACK
#include <memory>
#include <cassert>
#include <cmath>
#include <complex>

extern "C" {
  // scalapack routines
  void sl_init_(int*, const int*, const int*);
  void blacs_gridinfo_(const int*, const int*, const int*, int*, int*);
  void blacs_get_(const int*, const int*, int*);
  void blacs_gridmap_(int*, const int*, const int*, const int*, const int*);
  void blacs_gridexit_(const int*);
  void blacs_exit_(int*);
  int blacs_pnum_(const int*, const int*, const int*);

  int numroc_(const int* globalsize, const int* blocksize, const int* myrow, const int* startproc, const int* nproc);
  void descinit_(int* desc, const int* dimr, const int* dimc, const int* nbr, const int* nbc, const int* nsr, const int* nsc, const int* context, const int* ld, int* info);
  void pdelset_(double* mat, const int* i, const int* j, const int* desc, const double* a);
  void pdelget_(const char*, const char*, double* val, const double* mat, const int* i, const int* j, const int* desc);

  void pdgemm_(const char*, const char*, const int*, const int*, const int*, const double*, const double*, const int*, const int*, const int*,
               const double*, const int*, const int*, const int*, const double*, double*, const int*, const int*, const int*);
  void pzgemm_(const char*, const char*, const int*, const int*, const int*, const std::complex<double>*, const std::complex<double>*, const int*, const int*, const int*,
               const std::complex<double>*, const int*, const int*, const int*, const std::complex<double>*, std::complex<double>*, const int*, const int*, const int*);
  void pdsyev_(const char*, const char*, const int*, double*, const int*, const int*, const int*, double*, double*, const int*, const int*, const int*, double*, const int*, const int*);
  void pzheev_(const char*, const char*, const int*, std::complex<double>*, const int*, const int*, const int*, double*, std::complex<double>*, const int*, const int*, const int*,
               std::complex<double>*, const int*, double*, const int*, int*);
  void pdsyevd_(const char*, const char*, const int*, double*, const int*, const int*, const int*, double*, double*, const int*, const int*, const int*,
                double*, const int*, int*, const int*, const int*);
  void pzheevd_(const char*, const char*, const int*, std::complex<double>*, const int*, const int*, const int*, double*, std::complex<double>*, const int*, const int*, const int*,
                std::complex<double>*, const int*, double*, const int*, int*, const int*, int*);
}

static void sl_init_(int& i, const int j, const int k) { sl_init_(&i, &j, &k); }
static void blacs_gridinfo_(const int a, const int b, const int c, int& d, int& e) { blacs_gridinfo_(&a, &b, &c, &d, &e); }
static void blacs_get_(const int i, const int j, int& k) { blacs_get_(&i, &j, &k); }
static void blacs_gridmap_(int& i, const int* j, const int k, const int l, const int m) { blacs_gridmap_(&i, j, &k, &l, &m); }
static void blacs_gridexit_(const int i) { blacs_gridexit_(&i); }
static void blacs_exit_(int i) { blacs_exit_(&i); }
static int blacs_pnum_(const int a, const int b, const int c) { return blacs_pnum_(&a, &b, &c); }

static int numroc_(const int a, const int b, const int c, const int d, const int e)
  { return numroc_(&a, &b, &c, &d, &e); }
static void descinit_(int* a, const int b, const int c, const int d, const int e, const int f, const int g, const int h, const int i, int& j)
  { descinit_(a, &b, &c, &d, &e, &f, &g, &h, &i, &j); }
static void pdelset_(double* a, const int b, const int c, const int* d, const double e) { pdelset_(a, &b, &c, d, &e); }
static void pdelget_(const char* a, const char* b, double& val, const double* mat, const int i, const int j, const int* desc) { pdelget_(a, b, &val, mat, &i, &j, desc); }

static void pdgemm_(const char* transa, const char* transb, const int l, const int m, const int n, const double alpha, const double* a, const int ia, const int ja, const int* desca,
                    const double* b, const int ib, const int jb, const int* descb, const double beta, double* c, const int ic, const int jc, const int* descc) {
  pdgemm_(transa, transb, &l, &m, &n, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
}

static void pdgemm_(const char* transa, const char* transb, const int l, const int m, const int n, const double alpha, const double* a, const int* desca,
                    const double* b, const int* descb, const double beta, double* c, const int* descc) {
  pdgemm_(transa, transb, l, m, n, alpha, a, 1, 1, desca, b, 1, 1, descb, beta, c, 1, 1, descc);
}

static void pzgemm_(const char* transa, const char* transb, const int l, const int m, const int n, const std::complex<double> alpha,
                    const std::complex<double>* a, const int ia, const int ja, const int* desca,
                    const std::complex<double>* b, const int ib, const int jb, const int* descb, const std::complex<double> beta,
                    std::complex<double>* c, const int ic, const int jc, const int* descc) {
  pzgemm_(transa, transb, &l, &m, &n, &alpha, a, &ia, &ja, desca, b, &ib, &jb, descb, &beta, c, &ic, &jc, descc);
}

static void pzgemm_(const char* transa, const char* transb, const int l, const int m, const int n, const std::complex<double> alpha, const std::complex<double>* a, const int* desca,
                    const std::complex<double>* b, const int* descb, const std::complex<double> beta, std::complex<double>* c, const int* descc) {
  pzgemm_(transa, transb, l, m, n, alpha, a, 1, 1, desca, b, 1, 1, descb, beta, c, 1, 1, descc);
}

static void pdsyev_(const char* a, const char* b, const int dim, double* mat, const int* descm, double* eig, double* coeff, const int* descc, double* work, const int lwork, int& info) {
  const int one = 1;
  pdsyev_(a, b, &dim, mat, &one, &one, descm, eig, coeff, &one, &one, descc, work, &lwork, &info);
}

static void pzheev_(const char* a, const char* b, const int dim, std::complex<double>* mat, const int* descm, double* eig, std::complex<double>* coeff, const int* descc,
                    std::complex<double>* work, const int lwork, double* rwork, const int lrwork, int& info) {
  const int one = 1;
  pzheev_(a, b, &dim, mat, &one, &one, descm, eig, coeff, &one, &one, descc, work, &lwork, rwork, &lrwork, &info);
}

static void pdsyevd_(const char* a, const char* b, const int dim, double* mat, const int* descm, double* eig, double* coeff, const int* descc, double* work, const int lwork,
                     int* iwork, const int liwork, int& info) {
  const int one = 1;
  pdsyevd_(a, b, &dim, mat, &one, &one, descm, eig, coeff, &one, &one, descc, work, &lwork, iwork, &liwork, &info);
}

static void pzheevd_(const char* a, const char* b, const int dim, std::complex<double>* mat, const int* descm, double* eig, std::complex<double>* coeff, const int* descc,
                     std::complex<double>* work, const int lwork, double* rwork, const int lrwork, int* iwork, const int liwork, int& info) {
  const int one = 1;
  pzheevd_(a, b, &dim, mat, &one, &one, descm, eig, coeff, &one, &one, descc, work, &lwork, rwork, &lrwork, iwork, &liwork, &info);
}

static std::pair<int, int> numgrid(int numproc) {
  int sq = static_cast<int>(std::sqrt(static_cast<double>(numproc)))+1;
  for (int i = sq; i != 0; --i)
    if (numproc%i == 0) return {i, numproc/i};
  assert(false);
  return {0,0};
}

namespace bagel {
  const static int blocksize__ = 32;
}

#endif

#endif
