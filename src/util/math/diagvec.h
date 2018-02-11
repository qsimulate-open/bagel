//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: diagvec.h
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


#ifndef __SRC_MATH_DIAGVEC_H
#define __SRC_MATH_DIAGVEC_H

#include <src/util/math/matrix.h>
#include <src/util/math/vectorb.h>

namespace bagel {

class DiagVec {
  private:
    VectorB data_;

  public:
    DiagVec() { }
    DiagVec(const size_t n) : data_(n) { }
    DiagVec(const VectorB& o) : data_(o) { }
    DiagVec(const Matrix& m) : data_(m.mdim()) {
      for (int i = 0; i != size(); ++i)
        data_(i) = m(i, i);
    }

    size_t size() const { return data_.size(); }
    const VectorB& diag() const { return data_; }
    double& operator()(const int i) { return data_(i); }
    const double& operator()(const int i) const { return data_(i); }

};

extern Matrix operator*(const DiagVec&, const Matrix&);
extern Matrix operator*(const Matrix&, const DiagVec&);
extern DiagVec operator*(const DiagVec&, const DiagVec&);

}

#endif
