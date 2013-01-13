//
// BAGEL - Parallel electron correlation program.
// Filename: dfinttask.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Matthew Kelley <matthew.kelley@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can reblocksribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The BAGEL package is blocksributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the BAGEL package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//

#include <src/df/dfinttask.h>

using namespace std;
using namespace bagel;

DFIntTask_base::DFIntTask_base(array<shared_ptr<const Shell>,4>& a, vector<int>& b, vector<shared_ptr<DFBlock> >& df) : shell_(a), rank_(b.size()), dfblocks_(df) {
  int j = 0;
  for (auto i = b.begin(); i != b.end(); ++i, ++j) offset_[j] = *i;
}

void DFIntTask_base::compute() {

  shared_ptr<Integral> p = compute_batch(shell_);

  // all slot in
 for (int i = 0; i != nblocks(); ++i) {
    assert(dfblocks_[i]->b1size() == dfblocks_[i]->b2size());
    const size_t nbin = dfblocks_[i]->b1size();
    const size_t naux = dfblocks_[i]->asize();
    const double* ppt = p->data(i);
    if (rank_ == 3) {
      double* const data = dfblocks_[i]->get();
      for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
        for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1) {
          for (int j2 = offset_[2]; j2 != offset_[2] + shell_[1]->nbasis(); ++j2, ++ppt) {
            data[j2+naux*(j1+nbin*j0)] = data[j2+naux*(j0+nbin*j1)] = *ppt;
          }
        }
      }
    
#if 0
  } else if (rank_ == 2) {
    double* const data = df_->data2_.get();
    for (int j0 = offset_[0]; j0 != offset_[0] + shell_[2]->nbasis(); ++j0) {
      for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {
        data[j1+j0*naux] = data[j0+j1*naux] = *ppt;
      }
    }
#endif
    } else {
      assert(false);
    }
  }
}

