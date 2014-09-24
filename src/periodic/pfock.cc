//
// BAGEL - Parallel electron correlation program.
// Filename: pfock.cc
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


#include <src/periodic/pfock.h>

using namespace std;
using namespace bagel;

PFock::PFock(shared_ptr<const Lattice> l, shared_ptr<const ZMatrix> h, shared_ptr<const ZMatrix> c)
  : ZMatrix(l->primitive_cell()->nbasis(), l->primitive_cell()->nbasis()), lattice_(l), previous_(h), pcoeff_(c) {

  nblock_ = l->num_kpoints();
  blocksize_ = l->primitive_cell()->nbasis();

}

void PFock::ft() {

  shared_ptr<ZMatrix> tmp = make_shared<ZMatrix>(*this);
  zero();

  int k = 0;
  for (auto& kvec : lattice_->lattice_rvectors()) {
    shared_ptr<ZMatrix> kblock = make_shared<ZMatrix>(blocksize_, blocksize_);
    kblock->zero();
    int g = 0;
    for (auto& gvec : lattice_->lattice_vectors()) {
      complex<double> factor(0.0, gvec[0] * kvec[0] + gvec[1] * kvec[1] + gvec[2] * kvec[2]);
      factor = std::exp(factor);
      const int offset = g * blocksize_;
      shared_ptr<ZMatrix> gblock = tmp->get_submatrix(offset, offset, blocksize_, blocksize_);
      gblock->scale(factor);
      *kblock += *gblock;
      ++g;
    }
    add_block(complex<double>(1.0, 0.0), blocksize_ * k, blocksize_ * k,  blocksize_, blocksize_, kblock);
    ++k;
  }

}

void PFock::ift() {

}
