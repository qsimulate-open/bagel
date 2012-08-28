//
// BAGEL - Parallel electron correlation program.
// Filename: contraction.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/util/f77.h>
#include <src/rysint/rysint.h>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;

void RysInt::perform_contraction_new_outer(const int nsize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                       const vector<vector<double> >& coeff0, const vector<int>& upper0, const vector<int>& lower0, const int cdim0, 
                       const vector<vector<double> >& coeff1, const vector<int>& upper1, const vector<int>& lower1, const int cdim1) {
  const int worksize = nsize * pdim1; 
  double* const work = stack_->get(worksize);
  double* current_cont = cont;

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = lower0[i];
    const int end0   = upper0[i];
    fill(work, work + worksize, 0.0);
    for (int j = begin0; j != end0; ++j) 
      daxpy_(worksize, coeff0[i][j], &prim[j * worksize], 1, work, 1); 

    for (int k = 0; k != cdim1; ++k, current_cont += nsize) {
      const int begin1 = lower1[k];
      const int end1   = upper1[k];
      fill(current_cont, current_cont + nsize, 0.0); 
      for (int j = begin1; j != end1; ++j)
        daxpy_(nsize, coeff1[k][j], &work[j * nsize], 1, current_cont, 1); 
    }
  }

  stack_->release(worksize, work);
}


void RysInt::perform_contraction_new_inner(const int nsize, const int ac, const double* prim, const int pdim0, const int pdim1, double* cont, 
                       const vector<vector<double> >& coeff0, const vector<int>& upper0, const vector<int>& lower0, const int cdim0, 
                       const vector<vector<double> >& coeff1, const vector<int>& upper1, const vector<int>& lower1, const int cdim1) {
  const int worksize = pdim1 * ac;
  double* const work = stack_->get(worksize);
  double* current_cont = cont;

  for (int n = 0; n != nsize; ++n) { // loop of cdim * cdim
    const double* current_prim = &prim[ac * pdim1 * pdim0 * n];

    for (int i = 0; i != cdim0; ++i) {

      const int begin0 = lower0[i];
      const int end0   = upper0[i];
      fill(work, work + worksize,  0.0);
      for (int j = begin0; j != end0; ++j) 
        daxpy_(worksize, coeff0[i][j], &current_prim[j * worksize], 1, work, 1); 

      for (int k = 0; k != cdim1; ++k, current_cont += ac) {
        const int begin1 = lower1[k];
        const int end1   = upper1[k];
        fill(current_cont, current_cont + ac, 0.0); 
        for (int j = begin1; j != end1; ++j) {
          daxpy_(ac, coeff1[k][j], &work[j * ac], 1, current_cont, 1); 
        }
      }
    }
  } 
  stack_->release(worksize, work);
}


void RysInt::perform_contraction(const int asize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                 const vector<vector<double> >& coeff0, const vector<pair<int, int> >& ranges0, const int cdim0, 
                                 const vector<vector<double> >& coeff1, const vector<pair<int, int> >& ranges1, const int cdim1) {
  // transformation of index1
  const int worksize = pdim1 * asize;
  double* const work = stack_->get(worksize);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    std::fill(work, work + worksize, 0.0);
    for (int j = begin0; j != end0; ++j) 
      daxpy_(worksize, coeff0[i][j], &prim[j * worksize], 1, work, 1); 

    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      fill(cont, cont + asize, 0.0);
      for (int j = begin1; j != end1; ++j) {
        daxpy_(asize, coeff1[k][j], &work[j * asize], 1, cont, 1); 
      }
    }
  }
  stack_->release(worksize, work);
}


