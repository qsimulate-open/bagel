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
    std::shared_ptr<const Matrix> data_; // x,y,z,weight

    const size_t ngrid_;

    // basis functions and derivaties on this grid
    std::shared_ptr<Matrix> basis_;
    std::shared_ptr<Matrix> gradx_;
    std::shared_ptr<Matrix> grady_;
    std::shared_ptr<Matrix> gradz_;

    void init();

  public:
    DFTGridPoint(std::shared_ptr<const Geometry> g, std::shared_ptr<const Matrix>& o)
      : geom_(g), data_(o), ngrid_(o->mdim()) { assert(data_->ndim() == 4); init(); }

    std::shared_ptr<const Matrix> basis() const { return basis_; } 
    std::shared_ptr<const Matrix> gradx() const { return gradx_; } 
    std::shared_ptr<const Matrix> grady() const { return grady_; } 
    std::shared_ptr<const Matrix> gradz() const { return gradz_; } 
    const double& weight(const size_t i) const { return data_->element(3,i); }
    size_t size() const { return ngrid_; }
    std::shared_ptr<const Matrix> data() const { return data_; }
};


class DFTGrid_base {
  protected:
    const std::shared_ptr<const Geometry> geom_;
    std::shared_ptr<const DFTGridPoint> grid_;

    // TODO to be controlled by the input deck
    constexpr static double grid_thresh_ = 1.0e-8;

    double fuzzy_cell(std::shared_ptr<const Atom> a, std::array<double,3>&& x) const;
    void add_grid(const int nrad, const int nang, const std::unique_ptr<double[]>& r_ch, const std::unique_ptr<double[]>& w_ch,
                  const std::unique_ptr<double[]>& x, const std::unique_ptr<double[]>& y, const std::unique_ptr<double[]>& z, const std::unique_ptr<double[]>& w);

  public:
    DFTGrid_base(std::shared_ptr<const Geometry> geom) : geom_(geom) { }

    std::tuple<std::shared_ptr<const Matrix>,double> compute_xc(const std::string name, std::shared_ptr<const Matrix> mat) const;
};


// Becke-Chebyshev-Lebedev
class BLGrid : public DFTGrid_base {
  public:
    BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom);
};

// Treutler-Ahlrichs-Chebyshev-Lebedev 
class TALGrid : public DFTGrid_base {
  public:
    TALGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom);
};

// Pruned Grid
class DefaultGrid : public DFTGrid_base {
  public:
    DefaultGrid(std::shared_ptr<const Geometry> geom);
};

}

#endif
