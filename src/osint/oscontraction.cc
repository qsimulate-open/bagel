//
// Newint - Parallel electron correlation program.
// Filename: oscontraction.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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
#include <src/osint/osint.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <src/stackmem.h>

using namespace std;

extern StackMem* stack;

void OSInt::perform_contraction(const int asize, const double* prim, const int pdim0, const int pdim1, double* cont, 
                                       const vector<vector<double> >& coeff0, const vector<pair<int, int> >& ranges0, const int cdim0, 
                                       const vector<vector<double> >& coeff1, const vector<pair<int, int> >& ranges1, const int cdim1) {
  // transformation of index1
  const int worksize = pdim1 * asize;
  double* const work = stack->get(worksize);
  fill(cont, cont+asize*pdim0*pdim1, 0.0);

  for (int i = 0; i != cdim0; ++i) {
    const int begin0 = ranges0[i].first;
    const int end0   = ranges0[i].second;
    fill(work, work + worksize, 0.0);
    for (int j = begin0; j != end0; ++j) 
      daxpy_(worksize, coeff0[i][j], prim+j*worksize, 1, work, 1); 

    for (int k = 0; k != cdim1; ++k, cont += asize) {
      const int begin1 = ranges1[k].first;
      const int end1   = ranges1[k].second;
      for (int j = begin1; j != end1; ++j) {
        daxpy_(asize, coeff1[k][j], work+j*asize, 1, cont, 1); 
      }
    }
  }
  stack->release(worksize, work);
}


