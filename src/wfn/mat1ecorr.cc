//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: mat1ecorr.cc
// Copyright (C) 2017 Toru Shiozaki
//
// Author: Jae Woo Park <jwpk1201@northwestern.edu>
// Maintainer: Shiozaki group
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


#include <src/mat1e/hcore.h>
#include <src/wfn/mat1ecorr.h>

using namespace std;
using namespace bagel;


BOOST_CLASS_EXPORT_IMPLEMENT(Mat1eCorr)

vector<shared_ptr<Matrix>> Mat1eCorr::dkh_grad(shared_ptr<const Molecule> current) const {
  int natom = current->natom();
  vector<shared_ptr<Matrix>> dkhgrad;

  for (int i = 0; i != natom; ++i) {
    for (int j = 0; j != 3; ++j) {
      shared_ptr<Matrix> h_plus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = mat1e_dx();
        auto geom_plus = make_shared<Molecule>(*current, displ, false);
        auto hd_plus = make_shared<Hcore>(geom_plus, /* nodkh = */false);
        auto ho_plus = make_shared<Hcore>(geom_plus, /* nodkh = */true);

        h_plus = make_shared<Matrix>(*hd_plus - *ho_plus);
      }

      shared_ptr<Matrix> h_minus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = -mat1e_dx();
        auto geom_minus = make_shared<Molecule>(*current, displ, false);
        auto hd_minus = make_shared<Hcore>(geom_minus, /*nodkh = */false);
        auto ho_minus = make_shared<Hcore>(geom_minus, /*nodkh = */true);

        h_minus = make_shared<Matrix>(*hd_minus - *ho_minus);
      }

      dkhgrad.push_back(make_shared<Matrix>(*h_plus - *h_minus));
      dkhgrad[j+i*3]->scale(1.0 / (2.0 * mat1e_dx()));
    }
  }

  return dkhgrad;
}


shared_ptr<Matrix> Mat1eCorr::compute_grad_dkh(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) const {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3,natom);
  vector<shared_ptr<Matrix>> dkhg = dkh_grad(current);

  for (int i = 0; i != natom; ++i)
    for (int j = 0; j != 3; ++j)
      out->element(j,i) += dkhg[j+i*3]->dot_product(den);

  return out;
}


shared_ptr<Matrix> Mat1eCorr::compute_grad(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) const {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3, natom);

  if (dkh())
    out = compute_grad_dkh(current, den);

  return out;
}
