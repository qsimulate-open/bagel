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
    int block_;
    int size_block_;
    int size_final_;

  public:
    RealLondon(const std::array<std::shared_ptr<const Shell>,2>& input, const bool ReIm, const int nblocks, const int block) {
      batch_ = std::make_shared<Batch>(input);
      ReIm_ = ReIm;
      nblocks_ = nblocks;
      block_ = block;
      if (block_ > nblocks_) throw std::runtime_error ("You're asking for a block of data that does not exist - ??");

      // We could make this a friend class and pull size_block_ directly from the integral codes, but I don't want to modify them if we aren't keeping this
      const int prim01 = input[0]->num_primitive() * input[1]->num_primitive();
      const int ang0 = input[0]->angular_number();
      const int ang1 = input[1]->angular_number();
      const bool sph = input[0]->spherical();
      assert (sph == input[1]->spherical());
      const int amin = std::min(ang0, ang1);
      const int amax = ang0 + ang1;
      int asize = 0;
      for (int i=amin; i!=amax+1; i++) asize += (i+1)*(i+2)/2;
      const int asize_car = (ang0+1)*(ang0+2)*(ang1+1)*(ang1+2)/4;
      const int asize_sph = (2*ang0+1)*(2*ang1+1);
      size_block_ = prim01 * std::max(asize, asize_car);
      if (sph) { size_final_ = prim01 * asize_sph; }
      else size_final_ = prim01 * asize_car;
    }

    double* data() {return data_;}

    void compute() {
      batch_->compute();
      const std::complex<double>* cdata = batch_->data();
      data_ = (double*) std::malloc(size_final_*sizeof(double));
      const int i = block_-1;
      if (!ReIm_) {
        for (int j=0; j!=size_final_; j++) {
          data_[j] = cdata[i*size_block_ + j].real();
        }
      } else {
        for (int j=0; j!=size_final_; j++) {
          data_[j] = cdata[i*size_block_ + j].imag();
        }
      }
    }
};

std::complex<Matrix> transpose (std::complex<Matrix> S) {
  std::shared_ptr<Matrix> RT = S.real().transpose();
  std::shared_ptr<Matrix> IT = S.imag().transpose();
  std::complex<Matrix> out (*RT, *IT);
  return out;
}

std::complex<Matrix> multiply (std::complex<Matrix> A, std::complex<Matrix> B) {
  Matrix real = A.real() * B.real() - A.imag() * B.imag();
  Matrix imag = A.real() * B.imag() + A.imag() * B.real();
  std::complex<Matrix> out (real, imag);
  return out;
}

std::complex<Matrix> multiply (std::complex<Matrix> A, std::complex<Matrix> B, std::complex<Matrix> C) {
  std::complex<Matrix> D = multiply (A, B);
  std::complex<Matrix> out = multiply (D, C);
  return out;
}

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
