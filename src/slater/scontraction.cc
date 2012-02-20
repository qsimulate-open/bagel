//
// Newint - Parallel electron correlation program.
// Filename: scontraction.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/util/f77.h>
#include <src/slater/slaterbatch.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;

void SlaterBatch::perform_contraction_batch(const int nsize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                         const vector<vector<double> >& coeff0, const vector<int>& upper0, const vector<int>& lower0, const int cdim0, 
                                         const vector<vector<double> >& coeff1, const vector<int>& upper1, const vector<int>& lower1, const int cdim1) {

  // transformation of index1
  double* work = new double[cdim1 * pdim0];
  double* current_cont = cont;
  int poffset = 0;
  for (int n = 0; n != nsize; ++n) {
    double* work_pointer = work;
    for (int i = 0; i != pdim0; ++i, poffset += pdim1) { 
      const double* current_prim = &prim[poffset];
      for (int j = 0; j != cdim1; ++j, ++work_pointer) {
        const int begin = lower1[j];
        const int end   = upper1[j];
        *work_pointer = 0.0;
        for (int k = begin; k != end; ++k) {
          *work_pointer += coeff1[j][k] * current_prim[k]; 
        }
      }
    }
    // transformation of index2
    work_pointer = work;
    for (int i = 0; i != cdim0; ++i, current_cont += cdim1) {
      const int begin = lower0[i];
      const int end   = upper0[i];
      if (begin == end - 1) {
        memcpy(current_cont, &work_pointer[begin * cdim1], cdim1 * sizeof(double));
      } else {
        for (int j = begin; j != end; ++j) {
          const int unit = 1;
          daxpy_(&cdim1, &coeff0[i][j], &work_pointer[j * cdim1], &unit, current_cont, &unit); 
        }
      }
    }
  } 
  delete[] work;
}


void SlaterBatch::perform_contraction_new_outer(const int nsize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                             const vector<vector<double> >& coeff0, const vector<int>& upper0, const vector<int>& lower0, const int cdim0, 
                                             const vector<vector<double> >& coeff1, const vector<int>& upper1, const vector<int>& lower1, const int cdim1) {
  const int unit = 1;
  const int zeroint = 0;
  const double zero = 0.0;
  const int worksize = nsize * pdim1; 
  double* work = new double[worksize];
  double* current_cont = cont;

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = lower0[i];
    const int end0   = upper0[i];
    dcopy_(&worksize, &zero, &zeroint, work, &unit); 
    for (int j = begin0; j != end0; ++j) 
      daxpy_(&worksize, &coeff0[i][j], &prim[j * worksize], &unit, work, &unit); 

    for (int k = 0; k != cdim1; ++k, current_cont += nsize) {
      const int begin1 = lower1[k];
      const int end1   = upper1[k];
      fill(current_cont, current_cont + nsize, 0.0); 
      for (int j = begin1; j != end1; ++j)
        daxpy_(&nsize, &coeff1[k][j], &work[j * nsize], &unit, current_cont, &unit); 
    }
  }

  delete[] work;
}


void SlaterBatch::perform_contraction_new_inner(const int nsize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                             const vector<vector<double> >& coeff0, const vector<int>& upper0, const vector<int>& lower0, const int cdim0, 
                                             const vector<vector<double> >& coeff1, const vector<int>& upper1, const vector<int>& lower1, const int cdim1) {
  // transformation of index1
  const int unit = 1;
  const double zero = 0.0;
  const int zeroint = 0;
  const int ac = asize_ * csize_;
  const int worksize = pdim1 * ac;
  double* work = new double[worksize];
  double* current_cont = cont;

  for (int n = 0; n != nsize; ++n) { // loop of cdim * cdim
    const double* current_prim = &prim[ac * pdim1 * pdim0 * n];

    for (int i = 0; i != cdim0; ++i) {

      const int begin0 = lower0[i];
      const int end0   = upper0[i];
      dcopy_(&worksize, &zero, &zeroint, work, &unit); 
      for (int j = begin0; j != end0; ++j) 
        daxpy_(&worksize, &coeff0[i][j], &current_prim[j * worksize], &unit, work, &unit); 

      for (int k = 0; k != cdim1; ++k, current_cont += ac) {
        const int begin1 = lower1[k];
        const int end1   = upper1[k];
        fill(current_cont, current_cont + ac, 0.0); 
        for (int j = begin1; j != end1; ++j) {
          daxpy_(&ac, &coeff1[k][j], &work[j * ac], &unit, current_cont, &unit); 
        }
      }
    }
  } 
  delete[] work;
}


