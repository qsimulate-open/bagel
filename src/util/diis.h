//
// BAGEL - Parallel electron correlation program.
// Filename: diis.h
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and\/or modify
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


#ifndef __src_util_diis_h
#define __src_util_diis_h

#include <list>
#include <src/util/f77.h>
#include <memory>
#include <type_traits>
#include <src/util/matrix.h>
#include <src/util/zmatrix.h>

// std::shared_ptr<T> is assumed to be a shared_pointer of some class
// which have daxpy and ddot functions.
// T must have clone() function that returns shared_ptr<T>
// Mat should be either Matrix or ZMatrix. Needs Mat::solve function

namespace bagel {

template <class T, typename Mat = Matrix>
class DIIS {
  using Container_type_ = std::list<std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>>>;
  using iterator = typename Container_type_::iterator;

  protected:
    const int ndiis_;

    Container_type_ data_;

    std::shared_ptr<Mat> matrix_;
    std::shared_ptr<Mat> coeff_;

  public:
    DIIS(const int ndiis) : ndiis_(ndiis), matrix_(new Mat(ndiis+1, ndiis+1, true)), coeff_(new Mat(ndiis+1, 1)) { }

    ~DIIS() { }

    std::shared_ptr<T> extrapolate(const std::pair<std::shared_ptr<const T>, std::shared_ptr<const T>> input) {
      std::shared_ptr<const T> v = input.first;
      std::shared_ptr<const T> e = input.second;
      data_.push_back(input);

      if (data_.size() > ndiis_) {
        data_.pop_front();
        for (int i = 1; i != ndiis_; ++i) {
          for (int j = 1; j != ndiis_; ++j)
            matrix_->element(j-1, i-1) = matrix_->element(j, i);
        }
      }
      const int cnum = data_.size();
      iterator data_iter = data_.begin();

      // left hand side
      for (int i = 0; i != cnum - 1; ++i, ++data_iter) {
        matrix_->element(cnum-1, i) = matrix_->element(i, cnum-1) = e->ddot(*(data_iter->second));
        if (std::is_same<ZMatrix, Mat>::value)
          matrix_->element(i, cnum-1) = std::conj(matrix_->element(i, cnum-1)); 
      }
      matrix_->element(cnum-1, cnum-1)= e->ddot(e);
      for (int i = 0; i != cnum; ++i)
        matrix_->element(cnum, i) = matrix_->element(i, cnum) = -1.0;
      matrix_->element(cnum, cnum) = 0.0;

      // right hand side
      for (int i = 0; i != cnum; ++i)
        coeff_->element(i,0) = 0.0;
      coeff_->element(cnum,0) = -1.0;

      // solve the linear equation
      coeff_ = coeff_->solve(matrix_, cnum+1);

      // take a linear combination
      std::shared_ptr<T> out = input.first->clone();
      data_iter = data_.begin();
      for (int i = 0; i != cnum; ++i, ++data_iter)
        out->daxpy(coeff_->element(i,0), *(data_iter->first));

      return out;
    }

};

}

#endif

