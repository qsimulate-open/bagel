//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: csymmatrix.h
// Copyright (C) 2013 Toru Shiozaki
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

#ifndef __BAGEL_MATH_CSYMMATRIX_H
#define __BAGEL_MATH_CSYMMATRIX_H

#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>

namespace bagel {

// compressed symmetric matrix
class CSymMatrix {
  protected:
    bool localized_;
    int nocc_;
    size_t size_;
    std::unique_ptr<double[]> data_;

  private:
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int version) {
      boost::serialization::split_member(ar, *this, version);
    }
    template<class Archive>
    void save(Archive& ar, const unsigned int) const {
      ar << localized_ << nocc_ << size_ << make_array(data(), size());
    }
    template<class Archive>
    void load(Archive& ar, const unsigned int) {
      ar >> localized_ >> nocc_ >> size_;
      data_ = std::unique_ptr<double[]>(new double[size()]);
      ar >> make_array(data(), size());
    }

  public:
    CSymMatrix() { }
    CSymMatrix(int n, bool l = true) : localized_(l), nocc_(n), size_(n*(n+1)/2), data_(new double[size_]) { std::fill_n(data_.get(), size_, 0.0); }

    CSymMatrix(std::shared_ptr<const Matrix> in) : nocc_(in->ndim()), size_(nocc_*(nocc_+1)/2), data_(new double[size_]) {
      assert(in->ndim() == in->mdim() && (*in - *in->transpose()).rms() < 1.0e-8);
      for (int i = 0; i != nocc_; ++i)
        for (int j = 0; j <= i; ++j)
          element(j,i) = in->element(j,i);
    }

    // get info
    int nocc() const { return nocc_; }
    int size() const { return size_; }

    // sequential access
    double* data() { return data_.get(); }
    const double* data() const { return data_.get(); }
    double& data(int i) { assert(i < size_); return data_[i]; }
    const double& data(int i) const { assert(i < size_); return data_[i]; }

    // element access
    double& element(int i, int j) { assert(i <= j); return data_[i+((j*(j+1))>>1)]; }
    const double& element(int i, int j) const { assert(i <= j); return data_[i+((j*(j+1))>>1)]; }

    // returns a full matrix
    std::shared_ptr<Matrix> matrix() const {
      auto out = std::make_shared<Matrix>(nocc_, nocc_, localized_);
      for (int i = 0; i != nocc_; ++i) {
        for (int j = 0; j != i; ++j)
          out->element(i,j) = out->element(j,i) = element(j,i);
        out->element(i,i) = element(i,i);
      }
      return out;
    }

    std::shared_ptr<CSymMatrix> clone() const { return std::make_shared<CSymMatrix>(nocc_, localized_); }
    std::shared_ptr<CSymMatrix> copy() const {
      auto out = clone();
      std::copy_n(data(), size_, out->data());
      return out;
    }
};

}

#include <src/util/archive.h>
BOOST_CLASS_EXPORT_KEY(bagel::CSymMatrix)

#endif
