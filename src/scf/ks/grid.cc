//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: grid.cc
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

#include <src/scf/ks/grid.h>
#include <src/util/taskqueue.h>
#include <src/util/parallel/resources.h>

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
    std::shared_ptr<const Molecule> mol_;
  public:
    GridBasisTask(double* aa, double* bb, double* cc, double* dd, const double xx, const double yy, const double zz, std::shared_ptr<const Molecule> g)
    : a(aa), b(bb), c(cc), d(dd), x(xx), y(yy), z(zz), mol_(g) { }
    void compute() {
      int pos = 0;
      for (auto& i : mol_->atoms()) {
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
    std::shared_ptr<const Molecule> mol_;
  public:
    GridDeriv2Task(double* aa, double* bb, double* cc, double* dd, double* ee, double* ff,
                   const double xx, const double yy, const double zz, std::shared_ptr<const Molecule> g)
    : a(aa), b(bb), c(cc), d(dd), e(ee), f(ff), x(xx), y(yy), z(zz), mol_(g) { }
    void compute() {
      int pos = 0;
      for (auto& i : mol_->atoms()) {
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
  basis_ = make_shared<Matrix>(mol_->nbasis(), ngrid);
  gradx_ = make_shared<Matrix>(mol_->nbasis(), ngrid);
  grady_ = make_shared<Matrix>(mol_->nbasis(), ngrid);
  gradz_ = make_shared<Matrix>(mol_->nbasis(), ngrid);

  // TODO I guess this should be more efficient..
  TaskQueue<GridBasisTask> tasks(ngrid);
  for (size_t g = 0; g != ngrid; ++g) {
    tasks.emplace_back(basis_->element_ptr(0,g), gradx_->element_ptr(0,g), grady_->element_ptr(0,g), gradz_->element_ptr(0,g),
                                  data_->element(0,g), data_->element(1,g), data_->element(2,g), mol_);
  }
  tasks.compute();
}


array<shared_ptr<Matrix>,6> Grid::compute_grad2() const {
  array<shared_ptr<Matrix>,6> out;
  for (auto& i : out)
   i = make_shared<Matrix>(mol_->nbasis(), size());

  TaskQueue<GridDeriv2Task> tasks(size());
  for (size_t g = 0; g != size(); ++g) {
    tasks.emplace_back(out[0]->element_ptr(0,g), out[1]->element_ptr(0,g), out[2]->element_ptr(0,g),
                                   out[3]->element_ptr(0,g), out[4]->element_ptr(0,g), out[5]->element_ptr(0,g),
                                   data_->element(0,g), data_->element(1,g), data_->element(2,g), mol_);
  }
  tasks.compute();
  return out;
}
