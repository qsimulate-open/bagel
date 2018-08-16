//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: grid.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_KS_GRID_H
#define __SRC_KS_GRID_H

#include <array>
#include <memory>
#include <src/molecule/molecule.h>
#include <src/util/math/xyzfile.h>

namespace bagel {

class Grid {
  protected:
    const std::shared_ptr<const Molecule> mol_;
    const std::shared_ptr<const Matrix> data_; // x,y,z,weight

    // basis functions and derivaties on this grid
    std::shared_ptr<Matrix> basis_;
    std::shared_ptr<Matrix> gradx_;
    std::shared_ptr<Matrix> grady_;
    std::shared_ptr<Matrix> gradz_;

  public:
    Grid(std::shared_ptr<const Molecule> g, std::shared_ptr<const Matrix>& o)
      : mol_(g), data_(o) { assert(data_->ndim() == 4); }

    std::shared_ptr<const Matrix> basis() const { return basis_; }
    std::shared_ptr<const Matrix> gradx() const { return gradx_; }
    std::shared_ptr<const Matrix> grady() const { return grady_; }
    std::shared_ptr<const Matrix> gradz() const { return gradz_; }
    const double& weight(const size_t i) const { return data_->element(3,i); }
    size_t size() const { return data_->mdim(); }
    std::shared_ptr<const Matrix> data() const { return data_; }

    std::array<std::shared_ptr<Matrix>,6> compute_grad2() const;

    void init();

};

}

#endif
