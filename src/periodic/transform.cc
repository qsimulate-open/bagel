//
// BAGEL - Parallel electron correlation program.
// Filename: transform.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/transform.h>

using namespace std;
using namespace bagel;

Transform::Transform(const int blocksize, shared_ptr<Data> data, shared_ptr<KData> kdata, const vector<array<double, 3>> gvec, const vector<array<double, 3>> kvec)
  : nbasis_(blocksize), data_(data), kdata_(kdata), gvector_(gvec), kvector_(kvec) {

  num_gvector_ = gvec.size();
  num_kvector_ = kvec.size();

  assert(data_->nblock()  == nbasis_ * nbasis_ * num_gvector_);
  assert(kdata_->nblock() == nbasis_ * nbasis_ * num_kvector_);
  assert(data_->blocksize() == nbasis_ && kdata_->blocksize() == nbasis_);
}

// Slow discrete FT for now
void Transform::ft() {

  int k = 0;
  for (auto& kvec : kvector_) {
    auto kblock = make_shared<ZMatrix>(nbasis_, nbasis_);
    kblock->zero();
    int g = 0;
    for (auto& gvec : gvector_) {
      shared_ptr<Matrix> gblock = (*data_)[g];
      complex<double> factor(0.0, gvec[0]* kvec[0] + gvec[1] * kvec[1] + gvec[2] * kvec[2]);
      factor = std::exp(factor);
      auto tmp = make_shared<ZMatrix>(*gblock, factor);
      *kblock += *tmp;
      ++g;
    }
    (*kdata_)[k] = kblock;
    ++k;
  }
}

void Transform::ift() {

  data_->zero();

  int g = 0;
  for (auto& gvec : gvector_) {
    auto gblock = make_shared<Matrix>(nbasis_, nbasis_);
    gblock->zero();
    int k = 0;
    for (auto& kvec : kvector_) {
      shared_ptr<ZMatrix> kblock = (*kdata_)[k];
      complex<double> factor(0.0, -gvec[0]* kvec[0] - gvec[1] * kvec[1] - gvec[2] * kvec[2]);
      factor = std::exp(factor);
      auto tmp1 = make_shared<ZMatrix>(*(kblock->get_real_part()), factor);
      auto tmp2 = make_shared<ZMatrix>(*(kblock->get_imag_part()), factor);
      *gblock += *(tmp1->get_real_part()) + *(tmp2->get_real_part());
      ++k;
    }
    (*data_)[g] = gblock;
    ++g;
  }
}
