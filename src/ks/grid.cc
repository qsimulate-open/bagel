//
// BAGEL - Parallel electron correlation program.
// Filename: grid.cc
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

#include <src/ks/grid.h>
#include <src/parallel/resources.h>

namespace bagel {
class GridBasisTask {
  protected:
    double* const a;
    double* const b;
    double* const c;
    double* const d;
    const double x;
    const double y;
    const double z;
    std::shared_ptr<const Geometry> geom_;
  public:
    GridBasisTask(double* aa, double* bb, double* cc, double* dd, const double xx, const double yy, const double zz, std::shared_ptr<const Geometry> g)
    : a(aa), b(bb), c(cc), d(dd), x(xx), y(yy), z(zz), geom_(g) { }
    void compute() {
      int pos = 0;
      for (auto& i : geom_->atoms()) {
        // xyz coordinate relative to the atom i
        const double rx = x - i->position(0);
        const double ry = y - i->position(1);
        const double rz = z - i->position(2);
        for (auto& j : i->shells()) {
          j->compute_grid_value(a+pos, b+pos, c+pos, d+pos, rx, ry, rz);
          pos += j->nbasis();
        }
      }
    }
};
class GridDeriv2Task {
  protected:
    double* const a;
    double* const b;
    double* const c;
    double* const d;
    double* const e;
    double* const f;
    const double x;
    const double y;
    const double z;
    std::shared_ptr<const Geometry> geom_;
  public:
    GridDeriv2Task(double* aa, double* bb, double* cc, double* dd, double* ee, double* ff,
                   const double xx, const double yy, const double zz, std::shared_ptr<const Geometry> g)
    : a(aa), b(bb), c(cc), d(dd), e(ee), f(ff), x(xx), y(yy), z(zz), geom_(g) { }
    void compute() {
      int pos = 0;
      for (auto& i : geom_->atoms()) {
        // xyz coordinate relative to the atom i
        const double rx = x - i->position(0);
        const double ry = y - i->position(1);
        const double rz = z - i->position(2);
        for (auto& j : i->shells()) {
          j->compute_grid_value_deriv2(a+pos, b+pos, c+pos, d+pos, e+pos, f+pos, rx, ry, rz);
          pos += j->nbasis();
        }
      }
    }
};
}


using namespace std;
using namespace bagel;

void Grid::init() {
  const int ngrid = size();
  basis_ = make_shared<Matrix>(geom_->nbasis(), ngrid);
  gradx_ = make_shared<Matrix>(geom_->nbasis(), ngrid);
  grady_ = make_shared<Matrix>(geom_->nbasis(), ngrid);
  gradz_ = make_shared<Matrix>(geom_->nbasis(), ngrid);

  // TODO I guess this should be more efficient..
  vector<GridBasisTask> tasks;
  tasks.reserve(ngrid);
  for (size_t g = 0; g != ngrid; ++g) {
    tasks.push_back(GridBasisTask(basis_->element_ptr(0,g), gradx_->element_ptr(0,g), grady_->element_ptr(0,g), gradz_->element_ptr(0,g),
                                  data_->element(0,g), data_->element(1,g), data_->element(2,g), geom_));
  }
  TaskQueue<GridBasisTask> tq(tasks);
  tq.compute(resources__->max_num_threads());
}


array<shared_ptr<Matrix>,6> Grid::compute_grad2() const {
  array<shared_ptr<Matrix>,6> out;
  for (auto& i : out)
   i = make_shared<Matrix>(geom_->nbasis(), size());

  vector<GridDeriv2Task> tasks;
  tasks.reserve(size());
  for (size_t g = 0; g != size(); ++g) {
    tasks.push_back(GridDeriv2Task(out[0]->element_ptr(0,g), out[1]->element_ptr(0,g), out[2]->element_ptr(0,g),
                                   out[3]->element_ptr(0,g), out[4]->element_ptr(0,g), out[5]->element_ptr(0,g),
                                   data_->element(0,g), data_->element(1,g), data_->element(2,g), geom_));
  }
  TaskQueue<GridDeriv2Task> tq(tasks);
  tq.compute(resources__->max_num_threads());
  return out;
}
