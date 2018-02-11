//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diagmatrix.h
// Copyright (C) 2017 Nils Strand
//
// Author: Nils Strand <nilsstrand2022@u.northwestern.edu>
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


#ifndef __SRC_MATH_DIAGMATRIX_H
#define __SRC_MATH_DIAGMATRIX_H

#include <src/util/math/matrix.h>
#include <src/util/math/vectorb.h>

namespace bagel {

class DiagMatrix {
  private:
    VectorB data_;

  public:
    DiagMatrix() { }
    DiagMatrix(const size_t n) : data_(n) { }
    DiagMatrix(const VectorB& o) : data_(o) { }
    DiagMatrix(const DiagMatrix& o) : data_(o.data_) { }
    DiagMatrix(DiagMatrix&& o) : data_(std::move(o.data_)) { }
    DiagMatrix(const Matrix& m) : data_(m.mdim()) {
      assert(m.ndim() == m.mdim());
      for (int i = 0; i != ndim(); ++i)
        data_(i) = m(i, i);
    }

    double& operator()(const int i) { return data_(i); }
    const double& operator()(const int i) const { return data_(i); }

    DiagMatrix& operator+=(const DiagMatrix& o) { data_ += o.data_; return *this; }
    DiagMatrix& operator-=(const DiagMatrix& o) { data_ -= o.data_; return *this; }
    DiagMatrix operator+(const DiagMatrix& o) const { DiagMatrix out(*this); out += o; return out; }
    DiagMatrix operator-(const DiagMatrix& o) const { DiagMatrix out(*this); out -= o; return out; }

    size_t ndim() const { return data_.size(); }
    size_t mdim() const { return data_.size(); }

    const VectorB& diag() const { return data_; }
};

extern Matrix operator*(const DiagMatrix&, const Matrix&);
extern Matrix operator*(const Matrix&, const DiagMatrix&);
extern DiagMatrix operator*(const DiagMatrix&, const DiagMatrix&);

}

#endif
