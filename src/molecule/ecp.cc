//
// BAGEL - Parallel electron correlation program.
// Filename: ecp.cc
// Copyright (C) 2009 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// The BAGEL package is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 3, or (at your option)
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

#include <src/molecule/ecp.h>
using namespace std;
using namespace bagel;


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

  if (ishell_maxl_ > 0)
    for (int i = 0; i != 3; ++i)
      nr_[i] = count(shells_ecp_[ishell_maxl_]->ecp_r_power().begin(), shells_ecp_[ishell_maxl_]->ecp_r_power().end(), abs(i-2));

}

shared_ptr<const Shell_ECP> ECP::shell_maxl_ecp() const {

  shared_ptr<const Shell_ECP> shell_maxl;
  if (ishell_maxl_ < 0) {
    shell_maxl = make_shared<const Shell_ECP>(shells_ecp_[0]->position(), ecp_maxl_, vector<double>(1, 0.0), vector<double>(1, 0.0), vector<int>(1, 2));
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
