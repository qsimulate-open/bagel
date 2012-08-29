//
// BAGEL - Parallel electron correlation program.
// Filename: tildex.cc
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


#include <src/scf/tildex.h>
#include <src/util/f77.h>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <cmath>
#include <vector>

using namespace std;
using namespace bagel;

//#define USE_CANONICAL

TildeX::TildeX(const std::shared_ptr<Overlap> olp, const double thresh) : Matrix1e(olp->geom())  {

  // Use canonical orthogonalization (Szabo pp.144)
  nbasis_ = geom_->nbasis();
  ndim_ = nbasis_;
  const int size = nbasis_ * nbasis_;

  unique_ptr<double[]> eig(new double[nbasis_]);

  dcopy_(size, olp->data(), 1, data(), 1);

  {
    int info;
    const int lwork = 5 * nbasis_;
    unique_ptr<double[]> work(new double[lwork]);
    dsyev_("V", "L", ndim_, data(), ndim_, eig.get(), work.get(), lwork, info);
    if(info) throw runtime_error("dsyev in tildex failed.");
  }
  const double largest = fabs(eig[ndim_ - 1]);

  // counting how many orbital must be deleted owing to the linear dependency
  int cnt = 0;
  for (int i = 0; i != ndim_; ++i) {
    if (fabs(eig[i]) < largest * thresh) ++cnt;
    else break;
  }
  if (cnt != 0)
    cout << "  Caution: ignored " << cnt << " orbital" << (cnt == 1 ? "" : "s") << " in the orthogonalization." << endl << endl;

  for (int i = cnt; i != ndim_; ++i) {
#ifdef USE_CANONICAL
    const double scale = 1.0 / ::sqrt(eig[i]);
#else
    const double scale = 1.0 / ::sqrt(::sqrt(eig[i]));
#endif
    const int offset = i * ndim_;
    for (int j = 0; j != ndim_; ++j) {
      data_[j + offset] *= scale;
    }
  }
#ifdef USE_CANONICAL
  mdim_ = ndim_ - cnt;
  if (cnt != 0) {
    for (int i = 0; i != mdim_; ++i) {
      dcopy_(ndim_, data()+(i+cnt)*ndim_, 1, data()+i*ndim_, 1);
    }
  }
#else
  {
    mdim_ = ndim_;
    unique_ptr<double[]> tmp(new double[size]);
    dgemm_("N", "T", ndim_, ndim_, ndim_-cnt, 1.0, data()+cnt*ndim_, ndim_, data()+cnt*ndim_, ndim_, 0.0, tmp.get(), ndim_);
    dcopy_(size, tmp, 1, data_, 1);
  }
#endif


}


TildeX::~TildeX() {

}


