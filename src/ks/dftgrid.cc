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

#include <numeric>
#include <src/ks/dftgrid.h>
#include <src/ks/lebedevlist.h>
#include <src/util/constants.h>
#include <src/parallel/mpi_interface.h>

using namespace std;
using namespace bagel;

const static LebedevList lebedev;


void DFTGridPoint::init() {
  basis_ = unique_ptr<double[]>(new double[geom_->nbasis()]);
  gradx_ = unique_ptr<double[]>(new double[geom_->nbasis()]);
  grady_ = unique_ptr<double[]>(new double[geom_->nbasis()]);
  gradz_ = unique_ptr<double[]>(new double[geom_->nbasis()]);

  const double weight = data_[3]; 

  int pos = 0;
  for (auto& i : geom_->atoms()) {
    // xyz coordinate relative to the atom i
    const double x = data_[0] - i->position(0);
    const double y = data_[1] - i->position(1);
    const double z = data_[2] - i->position(2);
    for (auto& j : i->shells()) {
      // angular number
      j->compute_grid_value(basis_.get()+pos, gradx_.get()+pos, grady_.get()+pos, gradz_.get()+pos, x, y, z);
      pos += j->nbasis();
    }
  }
  assert(pos == geom_->nbasis());
}


double DFTGrid_base::integrate(std::shared_ptr<const Matrix> mat, const int power) {
  assert(power == 2);
  double sum = 0.0;
  for (int m = 0; m != mat->mdim(); ++m) {
    for (auto& i : grid_) {
      const double rho = ddot_(geom_->nbasis(), i->basis(), 1, mat->element_ptr(0,m), 1);
      const double gx = ddot_(geom_->nbasis(), i->gradx(), 1, mat->element_ptr(0,m), 1);
      const double gy = ddot_(geom_->nbasis(), i->grady(), 1, mat->element_ptr(0,m), 1);
      const double gz = ddot_(geom_->nbasis(), i->gradz(), 1, mat->element_ptr(0,m), 1);
      // this should give the number of electrons
//    sum += rho * rho * i->weight();
      // this should give the kinetic energy * 2
      sum += (gx*gx + gy*gy + gz*gz) * i->weight() * (-0.5);
    }
  }
  return sum;
}


double DFTGrid_base::fuzzy_cell(std::shared_ptr<const Atom> atom, array<double,3>&& xyz) const {
  double fuzzy = -1.0;
  double total = 0.0;
  for (auto& b : geom_->atoms()) {
    const double rbs1 = b->radius();
    double tmp = 1.0;
    for (auto& c : geom_->atoms()) {
      if (b != c) {
        const double rbs2 = c->radius();
        const double xi = sqrt(rbs1/rbs2); // sqrt. see JCP 102, 346 (1995)
        const double uij = (xi-1.0)/(xi+1.0);
        const double aij = uij / (uij*uij-1.0);

        const double distbc = b->distance(c); 
        const double distbg = b->distance(xyz);
        const double distcg = c->distance(xyz);
        const double muij = (distbg - distcg) / distbc;

        // see Becke's appendix
        double nuij = muij + aij*(1.0-muij*muij); // eq. a2
        for (int i = 0; i != 3; ++i)
          nuij = (1.5-0.5*nuij*nuij)*nuij; // eq. 19

        tmp *= 0.5*(1.0-nuij); // eq. 21
      }
    }
    if (b == atom) fuzzy = tmp; 
    total += tmp;
  }
  assert(fuzzy >= 0);
  return fuzzy / total; // Eq. 22
}


void DFTGrid_base::add_grid(const int nrad, const int nang, const unique_ptr<double[]>& r_ch, const unique_ptr<double[]>& w_ch,
                            const unique_ptr<double[]>& x, const unique_ptr<double[]>& y, const unique_ptr<double[]>& z, const unique_ptr<double[]>& w) {
  for (auto& a : geom_->atoms()) {
    const double rbs = a->radius();
    for (int i = 0; i != nrad; ++i) {
      for (int j = 0; j != nang; ++j) {
        const double xg = x[j] * r_ch[i] * rbs + a->position(0);
        const double yg = y[j] * r_ch[i] * rbs + a->position(1);
        const double zg = z[j] * r_ch[i] * rbs + a->position(2);
        double weight = w[j] * w_ch[i] * pow(rbs,3) * 4.0*pi__ * fuzzy_cell(a, array<double,3>{{xg, yg, zg}});

        // set to data 
        if (weight > grid_thresh_)
          grid_.push_back(shared_ptr<const DFTGridPoint>(new DFTGridPoint(geom_, array<double,4>{{xg, yg, zg, weight}})));
      }
    }
  }
}


// grid without 'pruning'. Becke's original mapping
BLGrid::BLGrid(const size_t nrad, const size_t nang, std::shared_ptr<const Geometry> geom) : DFTGrid_base(geom) {
  // first allocate Grids
  const size_t gridsize = nrad*nang*geom->natom();
  grid_.reserve(gridsize);

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
    const double t = cos((i+1)*pi__/(nrad+1)); 
    r_ch[i] = (1.0+t)/(1.0-t); 
    w_ch[i] = 2.0 / pow(1.0-t, 2.0)                  // due to mapping from [0,infty) to [-1, 1]
            * pi__/(nrad+1)*sin((i+1)*pi__/(nrad+1)) // Gauss-Chebyshev weight
            * r_ch[i]*r_ch[i];                       // due to r^2 in the spherical coordinate integration
  }

  add_grid(nrad, nang, r_ch, w_ch, x, y, z, w);
}
