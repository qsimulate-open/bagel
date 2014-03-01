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

  public:
    RealLondon(const std::array<std::shared_ptr<const Shell>,2>& input, const bool ReIm, const int nblocks, const int block) {
      batch_ = std::make_shared<Batch>(input);
      ReIm_ = ReIm;
      nblocks_ = nblocks;
      block_ = block;
      if (block_ > nblocks_ || block_ < 1) throw std::runtime_error ("You're asking for a block of data that does not exist - ??");
      size_block_ = batch_->size_block();
    }

    double* data() {return data_;}

    void compute() {
      batch_->compute();
      const std::complex<double>* cdata = batch_->data();
      data_ = (double*) std::malloc(size_block_*sizeof(double));
      const int i = block_-1;
      if (!ReIm_) {
        for (int j=0; j!=size_block_; j++) {
          data_[j] = cdata[i*size_block_ + j].real();
        }
      } else {
        for (int j=0; j!=size_block_; j++) {
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

std::complex<Matrix> scale (std::complex<double> c, std::complex<Matrix> A) {
  double a = c.real();
  double b = c.imag();
  Matrix Ra = A.real();
  Matrix Ia = A.imag();
  Matrix Rb = Ra;
  Matrix Ib = Ia;
  Ra.scale(a);
  Rb.scale(b);
  Ia.scale(a);
  Ib.scale(b);
  std::complex<Matrix> out ((Ra-Ib),(Rb+Ia));
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
