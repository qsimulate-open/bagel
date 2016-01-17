//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diis.h
// Copyright (C) 2009 Toru Shiozaki
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


#ifndef __SRC_MATH_DIIS_H
#define __SRC_MATH_DIIS_H

#include <type_traits>
#include <src/util/math/matrix.h>
#include <src/util/math/zmatrix.h>
#include <src/util/math/matop.h>
#include <src/util/serialization.h>

// std::shared_ptr<T> is assumed to be a shared_pointer of some class
// which have daxpy and ddot functions.
// T must have clone() function that returns shared_ptr<T>
// Mat should be either Matrix or ZMatrix. Needs Mat::solve function

namespace bagel {

template <class T, typename Mat = Matrix>
class DIIS {
  protected:
    int ndiis_;

    std::list<std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>>> data_;

    std::shared_ptr<Mat> matrix_;
    std::shared_ptr<Mat> coeff_;

  private:
    // serialization
    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive& ar, const unsigned int) {
      ar & ndiis_ & data_ & matrix_ & coeff_;
    }

  public:
    DIIS() { }
    DIIS(const int ndiis) : ndiis_(ndiis), matrix_(std::make_shared<Mat>(ndiis+1, ndiis+1, true)), coeff_(std::make_shared<Mat>(ndiis+1, 1)) { }

    std::shared_ptr<T> extrapolate(const std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>> input) {
      std::shared_ptr<const T> v = input.first;
      std::shared_ptr<const T> e = input.second;
      data_.push_back(input);

      if (data_.size() > ndiis_) {
        data_.pop_front();
        matrix_->copy_block(0, 0, ndiis_-1, ndiis_-1, matrix_->get_submatrix(1, 1, ndiis_-1, ndiis_-1));
      }
      const int cnum = data_.size();
      auto data_iter = data_.begin();

      // left hand side
      for (int i = 0; i != cnum - 1; ++i, ++data_iter) {
        matrix_->element(cnum-1, i) = e->dot_product(*(data_iter->second));
        matrix_->element(i, cnum-1) = detail::conj(matrix_->element(cnum-1, i));
      }
      matrix_->element(cnum-1, cnum-1)= e->dot_product(e);
      for (int i = 0; i != cnum; ++i)
        matrix_->element(cnum, i) = matrix_->element(i, cnum) = -1.0;
      matrix_->element(cnum, cnum) = 0.0;

      // right hand side
      for (int i = 0; i != cnum; ++i)
        coeff_->element(i,0) = 0.0;
      coeff_->element(cnum,0) = -1.0;

      // solve the linear equation
      coeff_ = coeff_->solve(matrix_, cnum+1);

      // return a linear combination
      std::shared_ptr<T> out = input.first->clone();
      data_iter = data_.begin();
      for (int i = 0; i != cnum; ++i, ++data_iter)
        out->ax_plus_y(coeff_->element(i,0), *(data_iter->first));
      return out;
    }

};

}

#endif

