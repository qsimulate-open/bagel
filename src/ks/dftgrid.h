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
    std::array<double,4> data; // x,y,z,weight

  public:
    DFTGridPoint() { }
    DFTGridPoint(const std::array<double,4>& o) : data(o) { }
    DFTGridPoint(const DFTGridPoint& o) : data(o.data) { }

    DFTGridPoint& operator=(const std::array<double,4>& o) { data = o; return *this; }
    DFTGridPoint& operator=(const DFTGridPoint& o) { data = o.data; return *this; }

    double* pointer(const int i = 0) { return &data[i]; }
};


class DFTGrid_base {
  protected:
    std::vector<DFTGridPoint> grid_;

  public:
    DFTGrid_base() { }
};


// Becke-Chebyshev-Lebedev
class BLGrid : public DFTGrid_base {
  protected:

  public:
    BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom);

};

}

#endif
