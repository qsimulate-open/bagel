//
// BAGEL - Parallel electron correlation program.
// Filename: matrix_interface.h
// Copyright (C) 2012 Toru Shiozaki
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

// The contents of this file are used to interface MixedBasis with the complex London integral codes

#ifndef __SRC_INTEGRAL_COMPOS_MATRIX_INTERFACE_H
#define __SRC_INTEGRAL_COMPOS_MATRIX_INTERFACE_H

#include <cstdlib>

namespace bagel {

template<typename Batch>
class RealLondon {
  protected:

    double* data_;
    std::shared_ptr<Batch> batch_;
    bool ReIm_;
    int nblocks_;
    int size_block_;

  public:
    RealLondon(const std::array<std::shared_ptr<const Shell>,2>& input, const bool ReIm, const int nblocks) {
      batch_ = std::make_shared<Batch>(input);
      ReIm_ = ReIm;
      nblocks_ = nblocks;

      // We could make this a friend class and pull size_block_ directly from the integral codes, but I don't want to modify them if we aren't keeping this
      const int prim01 = input[0]->num_primitive() * input[1]->num_primitive();
      const int ang0 = input[0]->angular_number();
      const int ang1 = input[1]->angular_number();
      const int amin = std::min(ang0, ang1);
      const int amax = ang0 + ang1;
      int asize = 0;
      for (int i=amin; i!=amax+1; i++) asize += (i+1)*(i+2)/2;
      const int asize_int = (ang0+1)*(ang0+2)*(ang1+1)*(ang1+2)/4;
      size_block_ = prim01 * std::max(asize, asize_int);
    }

    double* data() {return data_;}

    void compute() {
      batch_->compute();
      const std::complex<double>* cdata = batch_->data();
      data_ = (double*) std::malloc(nblocks_*size_block_*sizeof(double));
      if (!ReIm_) {
        for (int i=0; i!=nblocks_; i++) {
          for (int j=0; j!=size_block_; j++) {
            data_[i*nblocks_ + j] = cdata[i*nblocks_ + j].real();
          }
        }
      } else {
        for (int i=0; i!=nblocks_; i++) {
          for (int j=0; j!=size_block_; j++) {
            data_[i*nblocks_ + j] = cdata[i*nblocks_ + j].imag();
          }
        }
      }
    }
};

std::complex<Matrix> inverse (std::complex<Matrix> S) {
  Matrix S_Ri = S.real();
  S_Ri.inverse_symmetric();
  Matrix Si_R = S.real() + (S.imag() * S_Ri * S.imag());
  Si_R.inverse_symmetric();
  Matrix Si_I = S_Ri * S.imag() * Si_R;
  Si_I.scale(-1.0);
  std::complex<Matrix> out (Si_R, Si_I);
  return out;
}

}

#endif
