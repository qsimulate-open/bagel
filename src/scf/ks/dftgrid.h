//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: dftgrid.h
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

#ifndef __SRC_KS_DFTGRID_H
#define __SRC_KS_DFTGRID_H

#include <src/scf/ks/xcfunc.h>
#include <src/scf/ks/grid.h>

namespace bagel {

class DFTGrid_base {
  protected:
    const std::shared_ptr<const Molecule> mol_;
    std::shared_ptr<Grid> grid_;

    // TODO to be controlled by the input deck
    constexpr static double grid_thresh_ = 1.0e-10;

    void add_grid(const int nrad, const int nang, const std::unique_ptr<double[]>& r_ch, const std::unique_ptr<double[]>& w_ch,
                  const std::unique_ptr<double[]>& x, const std::unique_ptr<double[]>& y, const std::unique_ptr<double[]>& z, const std::unique_ptr<double[]>& w);
    void remove_redgrid();

    std::vector<std::shared_ptr<const Matrix>> compute_rho_sigma(std::shared_ptr<const XCFunc> func, std::shared_ptr<const Matrix> mat,
                                                    std::unique_ptr<double[]>& rho, std::unique_ptr<double[]>& sigma,
                                                    std::unique_ptr<double[]>& rhox, std::unique_ptr<double[]>& rhoy, std::unique_ptr<double[]>& rhoz) const;
  public:
    DFTGrid_base(std::shared_ptr<const Molecule> mol) : mol_(mol) { }

    std::tuple<std::shared_ptr<const Matrix>,double> compute_xc(std::shared_ptr<const XCFunc> func, std::shared_ptr<const Matrix> mat) const;
    std::shared_ptr<const GradFile> compute_xcgrad(std::shared_ptr<const XCFunc> func, std::shared_ptr<const Matrix> mat) const;
    double fuzzy_cell(std::shared_ptr<const Atom> a, std::array<double,3>&& x) const;

    std::shared_ptr<const Grid> grid() const { return grid_; }
};


// Becke-Chebyshev-Lebedev
class BLGrid : public DFTGrid_base {
  public:
    BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Molecule> mol);
};

// Treutler-Ahlrichs-Chebyshev-Lebedev
class TALGrid : public DFTGrid_base {
  public:
    TALGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Molecule> mol);
};

// Pruned Grid
class DefaultGrid : public DFTGrid_base {
  public:
    DefaultGrid(std::shared_ptr<const Molecule> mol);
};

}

#endif
