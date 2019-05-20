//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: xyzfile.h
// Copyright (C) 2012 Toru Shiozaki
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


#ifndef __SRC_MATH_XYZFILE_H
#define __SRC_MATH_XYZFILE_H

#include <src/util/math/matrix.h>
#include <src/util/math/matop.h>

namespace bagel {

class XYZFile : public Matrix {
  protected:

  public:
    XYZFile(const size_t natom, const double a = 0.0) : Matrix(3, natom, true) { fill(a); }
    XYZFile(const XYZFile& o) : Matrix(o) { assert(o.ndim() == 3 && localized_); }
    XYZFile(const Matrix& o) : Matrix(o) { localize(); assert(o.ndim() == 3); }

    XYZFile& operator=(const XYZFile& o) = default;
    XYZFile& operator=(XYZFile&& o) = default;

    std::shared_ptr<XYZFile> clone() const;

    void print(const std::string in = "", const int dum = 0) const override;
    void print_export(const std::string in = "", const int dum = 0) const;
    double maximum(int atomno) const;

    // this function assumes that double[] has data_.size()*data_size() elements.
    std::shared_ptr<XYZFile> transform(const std::shared_ptr<const Matrix>, const bool transpose) const;

};

// to make the code readable
using GradFile = XYZFile;

}

#endif
