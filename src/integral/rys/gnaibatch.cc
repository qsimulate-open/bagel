//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gnaibatch.cc
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


#include <src/integral/rys/gnaibatch.h>

using namespace std;
using namespace bagel;


GNAIBatch::GNAIBatch(const array<shared_ptr<const Shell>,2>& _info, const shared_ptr<const Molecule> mol, const tuple<int,int> i, shared_ptr<StackMem> stack)
  :  CoulombBatch_base(_info, mol, 1, 0, stack), iatom_(i) {
  if (swap01_) {
    swap(get<0>(iatom_), get<1>(iatom_));
  }
  set_exponents();
  const double integral_thresh = PRIM_SCREEN_THRESH;

  this->allocate_arrays(primsize_*natom_);
  compute_ssss(integral_thresh);
  root_weight(primsize_*natom_);
}


void GNAIBatch::set_exponents() {
  exponents_ = unique_ptr<double[]>(new double[primsize_*2]);
  assert(primsize_ == basisinfo_[0]->exponents().size() * basisinfo_[1]->exponents().size());
  double* tmp = exponents_.get();
  for (auto i0 = basisinfo_[0]->exponents().begin(); i0 != basisinfo_[0]->exponents().end(); ++i0) {
    for (auto i1 = basisinfo_[1]->exponents().begin(); i1 != basisinfo_[1]->exponents().end(); ++i1, tmp+=2) {
      tmp[0] = *i0;
      tmp[1] = *i1;
    }
  }
}


shared_ptr<GradFile> GNAIBatch::compute_gradient(shared_ptr<const Matrix> d, const int iatom0, const int iatom1, const int natom) const {
  auto out = make_shared<GradFile>(natom);
  for (int l = 0; l != natom; ++l)
    for (int k = 0; k != 3; ++k)
      out->element(k, l) += blas::dot_product(d->data(), d->size(), data_+size_block_*(k+3*l));
  return out;
}


void GNAIBatch::root_weight(const int ps) {
  if (amax_ + cmax_ == 0) {
    for (int j = 0; j != screening_size_; ++j) {
      int i = screening_[j];
      if (std::abs(T_[i]) < T_thresh__) {
        weights_[i] = 1.0;
      } else {
        const double sqrtt = sqrt(T_[i]);
        const double erfsqt = inline_erf(sqrtt);
        weights_[i] = erfsqt * sqrt(pi__) * 0.5 / sqrtt;
      }
    }
  } else {
    eriroot__.root(rank_, T_, roots_, weights_, ps);
  }
}
