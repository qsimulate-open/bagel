//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: hcoreinfo.cc
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
#include <src/mat1e/dkhcore.h>
#include <src/wfn/hcoreinfo.h>

using namespace std;
using namespace bagel;


HcoreInfo::HcoreInfo(shared_ptr<const PTree> idata) : type_(HcoreType::standard), gradtype_(true) {
  // DKH
  const bool dkh = idata->get<bool>("dkh", false);
  if (dkh)
    type_ = HcoreType::dkh;
  mat1e_dx_ = idata->get<double>("mat1e_dx", 0.001);

  // ECP
  const string basisfile = idata->get<string>("basis", "");
  const bool ecp = basisfile.find("ecp") != string::npos;
  if (ecp) {
    if (dkh)
      throw runtime_error("DKH and ECP cannot be used simultaneously");
    type_ = HcoreType::ecp;
  }
}

bool HcoreInfo::gradtype() const {
  assert(dkh());
  return gradtype_;
}

vector<shared_ptr<Matrix>> HcoreInfo::dkh_grad(shared_ptr<const Molecule> current) {
  int natom = current->natom();
  vector<shared_ptr<Matrix>> dkhgrad;

  for (int i = 0; i != natom; ++i) {
    for (int j = 0; j != 3; ++j) {
      shared_ptr<Matrix> h_plus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = mat1e_dx();
        auto geom_plus = make_shared<Molecule>(*current, displ, false);
        shared_ptr<Matrix> hd_plus = compute_dkh(geom_plus);
        auto ho_plus = make_shared<Hcore>(geom_plus);

        h_plus = make_shared<Matrix>(*hd_plus - *ho_plus);
      }

      shared_ptr<Matrix> h_minus;
      {
        auto displ = make_shared<XYZFile>(natom);
        displ->element(j,i) = -mat1e_dx();
        auto geom_minus = make_shared<Molecule>(*current, displ, false);
        shared_ptr<Matrix> hd_minus = compute_dkh(geom_minus);
        auto ho_minus = make_shared<Hcore>(geom_minus);

        h_minus = make_shared<Matrix>(*hd_minus - *ho_minus);
      }

      dkhgrad.push_back(make_shared<Matrix>(*h_plus - *h_minus));
      dkhgrad[j+i*3]->scale(1.0 / (2.0 * mat1e_dx()));
    }
  }

  return dkhgrad;
}


shared_ptr<Matrix> HcoreInfo::compute_grad_dkh(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3,natom);
  vector<shared_ptr<Matrix>> dkhg = dkh_grad(current);

  for (int i = 0; i != natom; ++i)
    for (int j = 0; j != 3; ++j)
      out->element(j,i) += dkhg[j+i*3]->dot_product(den);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute_grad(shared_ptr<const Molecule> current, shared_ptr<const Matrix> den) {
  int natom = current->natom();
  auto out = make_shared<Matrix>(3, natom);

  if (dkh())
    out = compute_grad_dkh(current, den);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute_dkh(shared_ptr<const Molecule> current) const {
  auto out = make_shared<DKHcore>(current);

  return out;
}


shared_ptr<Matrix> HcoreInfo::compute(shared_ptr<const Molecule> current) const {
  shared_ptr<Matrix> out;

  if (dkh())
    out = compute_dkh(current);

  return out;
}


void HcoreInfo::print() const {
  if (dkh())
    cout << "      - Using DKHcore" << endl;
}
