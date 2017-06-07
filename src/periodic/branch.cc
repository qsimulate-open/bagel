//
// BAGEL - Parallel electron correlation program.
// Filename: branch.cc
// Copyright (C) 2016 Toru Shiozaki
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


#include <src/periodic/branch.h>

using namespace bagel;
using namespace std;

Branch::Branch(const int ws, const array<double, 3>& c, const vector<shared_ptr<const ShellPair>>& sp) : iws_(ws), centre_(c), sp_(sp) {

  extent_ = 0.0;
  for (auto& i : sp_) {
    const double rr = sqrt(pow(i->centre(0)-c[0], 2.0) + pow(i->centre(1)-c[1], 2.0) + pow(i->centre(2)-c[2], 2.0)) + i->extent();
    extent_ = max(extent_, rr);
  }
}


bool Branch::is_neigh(shared_ptr<const Branch> b) const {

  double rr = 0;
  for (int i = 0; i != 3; ++i)
    rr += pow(centre_[i] - b->centre()[i], 2);

  const bool out = (sqrt(rr) <= 1.2*(extent_ + b->extent()));
  return out;

}


void Branch::get_neigh(const vector<shared_ptr<Branch>>& branch) {

  int nsp = 0;
  for (auto& b : branch)
    nsp += b->sp().size();

  neigh_.reserve(nsp);
  non_neigh_.reserve(nsp);
  for (auto& br : branch) {
    if (is_neigh(br)) {
      neigh_.insert(neigh_.end(), br->sp().begin(), br->sp().end());    
    } else {
      non_neigh_.insert(neigh_.end(), br->sp().begin(), br->sp().end());    
    }
  }
}
