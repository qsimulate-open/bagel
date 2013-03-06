//
// BAGEL - Parallel electron correlation program.
// Filename: dftgrid.h
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: NU theory
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

#ifndef __SRC_KS_DFTGRID_H
#define __SRC_KS_DFTGRID_H

#include <array>
#include <memory>
#include <src/wfn/geometry.h>

namespace bagel {

class DFTGridPoint {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::array<double,4> data_; // x,y,z,weight

    // basis functions and derivaties on this grid
    std::unique_ptr<double[]> basis_;
    std::unique_ptr<double[]> gradx_;
    std::unique_ptr<double[]> grady_;
    std::unique_ptr<double[]> gradz_;

    void init();

  public:
    DFTGridPoint(std::shared_ptr<const Geometry> g, const std::array<double,4>& o) : geom_(g), data_(o) { init(); }

    double* pointer(const int i = 0) { return &data_[i]; }
};


class DFTGrid_base {
  protected:
    std::vector<std::shared_ptr<const DFTGridPoint> > grid_;

  public:
    DFTGrid_base() { }
};


// Becke-Chebyshev-Lebedev
class BLGrid : public DFTGrid_base {
  protected:
    const std::shared_ptr<const Geometry> geom_;

  public:
    BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom);

};

}

#endif
