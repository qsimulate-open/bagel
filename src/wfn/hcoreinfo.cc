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


HcoreInfo::HcoreInfo(shared_ptr<const PTree> idata) : type_(HcoreType::standard) {
  // DKH
  const bool dkh = idata->get<bool>("dkh", false);
  if (dkh)
    type_ = HcoreType::dkh;

  // ECP
  const string basisfile = idata->get<string>("basis", "");
  const bool ecp = basisfile.find("ecp") != string::npos;
  if (ecp) {
    if (dkh)
      throw runtime_error("DKH and ECP cannot be used simultaneously");
    type_ = HcoreType::ecp;
  }
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
