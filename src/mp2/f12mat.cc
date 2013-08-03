//
// BAGEL - Parallel electron correlation program.
// Filename: f12mat.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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


#include <src/mp2/f12mat.h>

using namespace std;

void F12Mat::symmetrize(const bool braket) {
  // first electron 1-2 symmetry
  for (int i = 0; i != nocc_; ++i) {
    for (int j = 0; j != nocc_; ++j) {
      for (int k = 0; k != nocc_; ++k) {
        for (int l = 0; l != nocc_; ++l) {
          if ((i<<8) + j > (k<<8) + l) continue;
          const double ijkl = (data(l,k,j,i) + data(j,i,l,k)) * 0.5;
          data(l,k,j,i) = ijkl;
          data(j,i,l,k) = ijkl;
        }
      }
    }
  }
  // bra-ket symmetry
  if (braket) {
    for (int i = 0; i != nocc_; ++i) {
      for (int j = 0; j != nocc_; ++j) {
        for (int k = 0; k != nocc_; ++k) {
          for (int l = 0; l != nocc_; ++l) {
            if ((i<<8) + j > (k<<8) + l) continue;
            const double ijkl = (data(l,k,j,i) + data(j,i,l,k)) * 0.5;
            data(l,k,j,i) = ijkl;
            data(k,l,i,j) = ijkl;
          }
        }
      }
    }
  }

}

std::shared_ptr<F12Mat> F12Ten::contract(const std::shared_ptr<const F12Ten> o) const {
  if (dim0_ != o->dim0_ || dim1_ != o->dim0_) throw std::logic_error("bug in F12Ten::contract");
  std::unique_ptr<double[]> buf(new double[nocc_*nocc_*nocc_*nocc_]);
  dgemm_("T", "N", nocc_*nocc_, nocc_*nocc_, dim0_*dim1_, 1.0, data_, dim0_*dim1_, o->data_, dim0_*dim1_, 0.0, buf, nocc_*nocc_);
  return std::make_shared<F12Mat>(nocc_, std::move(buf));
};
