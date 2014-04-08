//
// BAGEL - Parallel electron correlation program.
// Filename: dfinttask_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/df/dfblock.h>
#include <src/df/df_london.h>
#include <src/df/dfinttask_london.h>
#include <src/molecule/shell.h>

using namespace std;
using namespace bagel;


void DFIntTask_London::compute() {
  std::shared_ptr<ComplexERIBatch> p = compute_batch(shell_);

  // all slot in
  for (int i = 0; i != N; ++i) {
    assert(dfblocks_[i]->b1size() == dfblocks_[i]->b2size());
    const size_t nbin = dfblocks_[i]->b1size();
    const size_t naux = dfblocks_[i]->asize();
    const std::complex<double>* ppt = p->data(i);
    std::complex<double>* const data = dfblocks_[i]->get();
    for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
      for (int j1 = offset_[1]; j1 != offset_[1] + shell_[2]->nbasis(); ++j1, ppt += shell_[1]->nbasis()) {
        for (int n=0; n!=shell_[1]->nbasis(); n++) {
          data[offset_[2]+naux*(j1+nbin*j0)+n] = ppt[n];
          if (N == 1) data[offset_[2]+naux*(j0+nbin*j1)+n] = conj(ppt[n]); // This shortcut only works with a real auxiliary basis
        }
      }
    }
  }
}


void DFIntTask_OLD_London::compute() {

  std::pair<const double*, std::shared_ptr<RysIntegral<double, Int_t::Standard>>> p = df_->compute_batch(shell_);
  const double* ppt = p.first;

  const size_t naux = df_->naux();
  // all slot in
  if (rank_ == 2) {
    double* const data = data_->data();
    for (int j0 = offset_[0]; j0 != offset_[0] + shell_[3]->nbasis(); ++j0) {
      for (int j1 = offset_[1]; j1 != offset_[1] + shell_[0]->nbasis(); ++j1, ++ppt) {
        data[j1+j0*naux] = data[j0+j1*naux] = *ppt; // Equivalent for a real auxiliary basis
        //data[j0+j1*naux] = *ppt;
        //data[j1+j0*naux] = std::conj(*ppt);
      }
    }
  } else {
    assert(false);
  }
}

