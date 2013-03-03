//
// BAGEL - Parallel electron correlation program.
// Filename: dftgrid.cc
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

#include <src/ks/dftgrid.h>
#include <src/ks/lebedevlist.h>
#include <src/util/constants.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

const static LebedevList lebedev;

// grid without 'pruning'. Becke's original mapping
BLGrid::BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom) {
  // first allocate Grids
  const size_t gridsize = nrad*nang*geom->natom();
  grid_ = vector<DFTGridPoint>(gridsize);

  // construct Lebedev grid
  unique_ptr<double[]> x(new double[nang]);
  unique_ptr<double[]> y(new double[nang]);
  unique_ptr<double[]> z(new double[nang]);
  unique_ptr<double[]> w(new double[nang]);
  lebedev.root(nang, x.get(), y.get(), z.get(), w.get());

  // construct Chebyshev grid 
  unique_ptr<double[]> r_ch(new double[nrad]);
  unique_ptr<double[]> w_ch(new double[nrad]);
  for (int i = 0; i != nrad; ++i) {
    const double t = cos(i*pi__/(nrad+1)); 
    r_ch[i] = (1.0+t)/(1.0-t); 
    w_ch[i] = 2.0 / pow(1.0-t, 2.0)              // due to mapping from [0,infty) to [-1, 1]
            * pi__/(nrad+1)*sin(i*pi__/(nrad+1)) // Gauss-Chebyshev weight
            * r_ch[i]*r_ch[i];                   // due to r^2 in the spherical coordinate integration
  }

  auto iter = grid_.begin(); 

  for (auto& i : geom->atoms()) {
    // TODO the value of R_BS is atom-specific
    const double rbs = 1.0;
    for (int i = 0; i != nrad; ++i) {
      for (int j = 0; j != nang; ++j) {
        const double xg = x[j] * r_ch[i] * rbs;
        const double yg = y[j] * r_ch[i] * rbs;
        const double zg = z[j] * r_ch[i] * rbs;
        const double weight = w[j] * w_ch[i] * pow(rbs,3.0);
        // TODO multiply fuzzy cell factors

        // set to data 
        *iter++ = array<double,4>{{xg, yg, zg, weight}};
      }
    }
  }
  assert(iter == grid_.end());

  // to make sure about alignment
  static_assert(sizeof(DFTGridPoint) == 4*sizeof(double), "something is wrong with the alignment of the DFTGrid class");
  mpi__->allreduce(grid_[0].pointer(), grid_.size()*4); 
}
