//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diagmatrix.h
// Copyright (C) 2017 Toru Shiozaki
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
    DiagMatrix(VectorB&& o) : data_(std::move(o)) { }
    DiagMatrix(const DiagMatrix& o) : data_(o.data_) { }
    DiagMatrix(DiagMatrix&& o) : data_(std::move(o.data_)) { }

    double& operator()(const int i) { return data_(i); }
    const double& operator()(const int i) const { return data_(i); }

    DiagMatrix& operator=(const DiagMatrix& o) { data_ = o.data_; return *this; }
    DiagMatrix& operator=(DiagMatrix&& o) { data_ = std::move(o.data_); return *this; }
    DiagMatrix& operator+=(const DiagMatrix& o) { data_ += o.data_; return *this; }
    DiagMatrix& operator-=(const DiagMatrix& o) { data_ -= o.data_; return *this; }
    DiagMatrix operator+(const DiagMatrix& o) const { DiagMatrix out(*this); out += o; return out; }
    DiagMatrix operator-(const DiagMatrix& o) const { DiagMatrix out(*this); out -= o; return out; }

    size_t ndim() const { return data_.size(); }
    size_t mdim() const { return data_.size(); }

    const VectorB& diag() const { return data_; }
};

namespace {

Matrix operator*(const DiagMatrix& v, const Matrix& m) {
  assert(m.ndim() == v.ndim());
  Matrix out(v.ndim(), v.mdim());
  for (int i = 0; i != m.mdim(); ++i)
    for (int j = 0; j != m.ndim(); ++j)
      out(j, i) = m(j, i) * v(j);
  return out;
}

Matrix operator*(const Matrix& m, const DiagMatrix& v) {
  assert(m.mdim() == v.mdim());
  Matrix out = m;
  for (int i = 0; i != m.mdim(); ++i)
    blas::scale_n(v(i), out.element_ptr(0, i), m.ndim());
  return out;
}

DiagMatrix operator*(const DiagMatrix& v1, const DiagMatrix& v2) {
  assert(v1.ndim() == v2.ndim());
  DiagMatrix out(v1.ndim());
  for (int i = 0; i != v1.ndim(); ++i)
    out(i) = v1(i) * v2(i);
  return out;
}

Matrix& operator+=(Matrix& out, const DiagMatrix& v) {
  assert(out.mdim() == v.mdim());
  for (int i = 0; i != out.mdim(); ++i)
    out(i, i) += v(i);
  return out;
}

Matrix& operator-=(Matrix& out, const DiagMatrix& v) {
  assert(out.mdim() == v.mdim());
  for (int i = 0; i != out.mdim(); ++i)
    out(i, i) -= v(i);
  return out;
}

Matrix operator+(const Matrix& m, const DiagMatrix& v) {
  Matrix out = m;
  out += v;
  return out;
}

Matrix operator-(const Matrix& m, const DiagMatrix& v) {
  Matrix out = m;
  out -= v;
  return out;
}

Matrix operator+(const DiagMatrix& v, const Matrix& m) { return m + v; }
Matrix operator-(const DiagMatrix& v, const Matrix& m) { return m - v; }


}

}

#endif
