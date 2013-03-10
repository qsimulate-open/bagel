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
    const double* const basis() const { return basis_.get(); } 
    const double* const gradx() const { return gradx_.get(); } 
    const double* const grady() const { return grady_.get(); } 
    const double* const gradz() const { return gradz_.get(); } 
    const double& weight() const { return data_[3]; }
};


class DFTGrid_base {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::vector<std::shared_ptr<const DFTGridPoint>> grid_;

    // TODO to be controlled by the input deck
    constexpr static double grid_thresh_ = 0.0;

    double fuzzy_cell(std::shared_ptr<const Atom> a, std::array<double,3>&& x) const;
    void add_grid(const int nrad, const int nang, const std::unique_ptr<double[]>& r_ch, const std::unique_ptr<double[]>& w_ch,
                  const std::unique_ptr<double[]>& x, const std::unique_ptr<double[]>& y, const std::unique_ptr<double[]>& z, const std::unique_ptr<double[]>& w);

  public:
    DFTGrid_base(std::shared_ptr<const Geometry> geom) : geom_(geom) { }

    double integrate(std::shared_ptr<const Matrix> mat, const int power = 2);
};


// Becke-Chebyshev-Lebedev
class BLGrid : public DFTGrid_base {
  protected:

  public:
    BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom);

};

}

#endif
