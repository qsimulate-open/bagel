//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: ecp.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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

#include <iostream>
#include <src/molecule/ecp.h>

using namespace std;
using namespace bagel;

ECP::ECP() : ecp_ncore_(0), ecp_maxl_(0), shells_ecp_(1, make_shared<const Shell_ECP>()), ishell_maxl_(-1) {}

ECP::ECP(const int ncore, const int maxl, vector<shared_ptr<const Shell_ECP>> shells_ecp)
  : ecp_ncore_(ncore), ecp_maxl_(maxl), shells_ecp_(shells_ecp) {
  nshell_ = shells_ecp_.size();
  get_shell_maxl_ecp();
}

void ECP::get_shell_maxl_ecp() {
  ishell_maxl_ = -1;
  for (auto ish = shells_ecp_.begin(); ish != shells_ecp_.end(); ++ish)
    if ((*ish)->angular_number() == ecp_maxl_) {
      ishell_maxl_ = distance(shells_ecp_.begin(), ish);
      break;
    }

  if (ishell_maxl_ != -1)
    for (int i = 0; i != 3; ++i)
      nr_[i] = count(shells_ecp_[ishell_maxl_]->ecp_r_power().begin(), shells_ecp_[ishell_maxl_]->ecp_r_power().end(), abs(i-2));

}

shared_ptr<const Shell_ECP> ECP::shell_maxl_ecp() const {

  shared_ptr<const Shell_ECP> shell_maxl;
  if (ishell_maxl_ < 0) {
    shell_maxl = make_shared<const Shell_ECP>();
  } else {
    shell_maxl = shells_ecp_[ishell_maxl_];
  }
  return shell_maxl;

}

void ECP::print() const {
  cout << "+++ ECP Parameters +++" << endl;
  cout << "Number of core electrons = " << ecp_ncore_ << endl;
  cout << "Max angular number       = " << ecp_maxl_ << endl;
  cout << "Number of ECP shells     = " << nshell_ << endl;
  for (auto& i : shells_ecp_) cout << i->show() << endl;
}
