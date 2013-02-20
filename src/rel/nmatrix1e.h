//
// BAGEL - Parallel electron correlation program.
// Filename: nmatrix1e.h
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


#ifndef __src_rel_nmatrix1e_h
#define __src_rel_nmatrix1e_h

#include <cassert>
#include <src/wfn/shell.h>
#include <src/wfn/geometry.h>
#include <src/util/matrix.h>
#include <string>
#include <algorithm>
#include <memory>

namespace bagel {

// specialized matrix for 1e integrals
class NMatrix1e {
  protected:
    std::shared_ptr<const Geometry> geom_;

    virtual void init() = 0;

    std::vector<std::shared_ptr<Matrix>> matrix_data_;

  public:
    const std::shared_ptr<Matrix>& operator[](const int i) const { return matrix_data_[i]; }
    const std::shared_ptr<Matrix> data(const int i) const { return matrix_data_[i]; }

    NMatrix1e(const std::shared_ptr<const Geometry> geom) : geom_(geom) { }

    virtual void print() const = 0;

    const std::shared_ptr<const Geometry> geom() const { return geom_; }

    std::vector<std::shared_ptr<Matrix>> data() const { return matrix_data_; }

};

}

#endif
