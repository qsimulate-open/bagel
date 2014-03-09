//
// Filename: matrix_interface.h
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
//
// The contents of this file are used to interface MixedBasis with the complex London integral codes.
// It could conceivably be useful later on, but is intended as a crude patch to facilitate debugging.
//

#ifndef __SRC_INTEGRAL_COMPRYS_TEST_CODES_MATRIX_INTERFACE_H
#define __SRC_INTEGRAL_COMPRYS_TEST_CODES_MATRIX_INTERFACE_H

#include <cstdlib>

namespace ryan {

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
    RealLondon(const std::array<std::shared_ptr<const bagel::Shell>,2>& input, const bool ReIm, const int nblocks, const int block) {
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

std::complex<bagel::Matrix> transpose (std::complex<bagel::Matrix> S) {
  std::shared_ptr<bagel::Matrix> RT = S.real().transpose();
  std::shared_ptr<bagel::Matrix> IT = S.imag().transpose();
  std::complex<bagel::Matrix> out (*RT, *IT);
  return out;
}

std::complex<bagel::Matrix> scale (std::complex<double> c, std::complex<bagel::Matrix> A) {
  // Note: The accuracy of this function has not been checked!
  double a = c.real();
  double b = c.imag();
  bagel::Matrix Ra = A.real();
  bagel::Matrix Ia = A.imag();
  bagel::Matrix Rb = Ra;
  bagel::Matrix Ib = Ia;
  Ra.scale(a);
  Rb.scale(b);
  Ia.scale(a);
  Ib.scale(b);
  std::complex<bagel::Matrix> out ((Ra-Ib),(Rb+Ia));
  return out;
}

std::complex<bagel::Matrix> multiply (std::complex<bagel::Matrix> A, std::complex<bagel::Matrix> B) {
  bagel::Matrix real = A.real() * B.real() - A.imag() * B.imag();
  bagel::Matrix imag = A.real() * B.imag() + A.imag() * B.real();
  std::complex<bagel::Matrix> out (real, imag);
  return out;
}

std::complex<bagel::Matrix> multiply (std::complex<bagel::Matrix> A, std::complex<bagel::Matrix> B, std::complex<bagel::Matrix> C) {
  std::complex<bagel::Matrix> D = multiply (A, B);
  std::complex<bagel::Matrix> out = multiply (D, C);
  return out;
}

std::complex<bagel::Matrix> inverse (std::complex<bagel::Matrix> S) {
  bagel::Matrix S_Ri = S.real();
  S_Ri.inverse();
  bagel::Matrix Si_R = S.real() + (S.imag() * S_Ri * S.imag());
  Si_R.inverse();
  bagel::Matrix Si_I = S_Ri * S.imag() * Si_R;
  Si_I.scale(-1.0);
  std::complex<bagel::Matrix> out (Si_R, Si_I);
  return out;
}

}

#endif
