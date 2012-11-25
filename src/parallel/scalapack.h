//
// BAGEL - Parallel electron correlation program.
// Filename: scalapack.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#ifndef __SRC_PARALLEL_SCALAPACK_H
#define __SRC_PARALLEL_SCALAPACK_H

#ifdef HAVE_SCALAPACK

extern "C" {
  // scalapck routines
  void sl_init_(int*, const int*, const int*);
  void blacs_gridinfo_(const int*, const int*, const int*, int*, int*); 
  void blacs_gridexit_(const int*);
  void blacs_exit_(int*);

  int numroc_(const int* globalsize, const int* blocksize, const int* myrow, const int* startproc, const int* nproc); 
  void descinit_(int* desc, const int* dimr, const int* dimc, const int* nbr, const int* nbc, const int* nsr, const int* nsc, const int* context, const int* ld, int* info);
  void pdelset_(double* mat, const int* i, const int* j, const int* desc, const double* a); 

  void pdsyev_(const char*, const char*, const int*, double*, const int*, const int*, const int*, double*, double*, const int*, const int*, const int*, double*, const int*, const int*); 
}

static int numroc_(const int a, const int b, const int c, const int d, const int e)
  { return numroc_(&a, &b, &c, &d, &e); } 
static void descinit_(int* a, const int b, const int c, const int d, const int e, const int f, const int g, const int h, const int i, int& j)
  { descinit_(a, &b, &c, &d, &e, &f, &g, &h, &i, &j); }
static void pdelset_(double* a, const int b, const int c, const int* d, const double e) { pdelset_(a, &b, &c, d, &e); }

static void pdsyev_(const char* a, const char* b, const int dim, double* mat, const int* descm, double* eig, double* coeff, const int* descc, double* work, const int lwork, int& info) {
  const int one = 1;
  pdsyev_(a, b, &dim, mat, &one, &one, descm, eig, coeff, &one, &one, descc, work, &lwork, &info);
}

std::pair<int, int> numgrid(int numproc) {
  int sq = static_cast<int>(sqrt(static_cast<double>(numproc)))+1;
  for (int i = sq; i != 0; --i)
    if (numproc%i == 0) return make_pair(i, numproc/i); 
  assert(false);
  return make_pair(0,0);
}

#endif

#endif
