//
// BAGEL - Parallel electron correlation program.
// Filename: grid.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
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

#ifndef __SRC_KS_GRID_H
#define __SRC_KS_GRID_H

#include <array>
#include <memory>
#include <src/wfn/geometry.h>
#include <src/grad/gradfile.h>

namespace bagel {

class Grid {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    const std::shared_ptr<const Matrix> data_; // x,y,z,weight

    // basis functions and derivaties on this grid
    std::shared_ptr<Matrix> basis_;
    std::shared_ptr<Matrix> gradx_;
    std::shared_ptr<Matrix> grady_;
    std::shared_ptr<Matrix> gradz_;

  public:
    Grid(std::shared_ptr<const Geometry> g, std::shared_ptr<const Matrix>& o)
      : geom_(g), data_(o) { assert(data_->ndim() == 4); }

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
