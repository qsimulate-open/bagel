//
// BAGEL - Parallel electron correlation program.
// Filename: gradfile.h
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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


#ifndef __SRC_GRAD_GRADFILE_H
#define __SRC_GRAD_GRADFILE_H

// a class for using the BFGS solver, which requires
// "clone, ddot and daxpy, along with overloaded operators and a copy constructor"

#include <cassert>
#include <iomanip>
#include <src/math/matrix.h>
#include <src/wfn/geometry.h>
#include <src/util/f77.h>

namespace bagel {

class GradFile : public Matrix {
  protected:

  public:
    GradFile(const size_t natom, const double a = 0.0) : Matrix(3, natom, true) { fill(a); }
    GradFile(const GradFile& o) : Matrix(o) { assert(o.ndim() == 3 && localized_); }
    GradFile(const Matrix& o) : Matrix(o) { localize(); assert(o.ndim() == 3); }

    std::shared_ptr<GradFile> clone() const;

    void print() const;

    // this function assumes that double[] has data_.size()*data_size() elements.
    std::shared_ptr<GradFile> transform(const std::shared_ptr<const Matrix>, const bool transpose) const;

};

}

#endif
