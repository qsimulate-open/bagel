//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: multisite.cc
// Copyright (C) 2014 Shane Parker
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/multisite/multisite.h>
#include <src/mat1e/overlap.h>
#include <src/util/muffle.h>

using namespace std;
using namespace bagel;

///> this constructor takes localized HF orbitals with a clear definition of active subspaces
MultiSite::MultiSite(shared_ptr<const PTree> input, shared_ptr<const Reference> ref, const int nsites) : input_(input), hf_ref_(ref), nsites_(nsites) {
  cout << string(60, '=') << endl;
  cout << string(15, ' ') << "Construct MultiSite" << endl;
  cout << string(60, '=') << endl;

  charge_ = input_->get<int>("charge", 0);
  nspin_ = input_->get<int>("spin", 0);

  active_electrons_ = input_->get_vector<int>("active_electrons");
  active_sizes_ = input_->get_vector<int>("active_sizes");

  const int nactele = accumulate(active_electrons_.begin(), active_electrons_.end(), 0);
  nclosed_ = (hf_ref_->geom()->nele() - charge_ - nactele) / 2;
  assert((hf_ref_->geom()->nele() - charge_ - nactele) % 2 == 0);

  nactive_ = accumulate(active_sizes_.begin(), active_sizes_.end(), 0);
}

void MultiSite::compute() {
  set_active();
}


void MultiSite::set_active() {
  auto active = input_->get_child("active");
  if (!active)
    throw logic_error("\"active\" keyword must be specified for multisite construction");
  if (active->size() != nsites_)
    throw logic_error("Must specify active spaces for all of the sites");

  vector<set<int>> active_orbitals;
  set<int> active_set;
  for (auto site : *active) {
    vector<int> active_orbs = site->get_vector<int>("");
    for_each(active_orbs.begin(), active_orbs.end(), [](int& x) { x--; });
    active_orbitals.emplace_back(active_orbs.begin(), active_orbs.end());
    active_set.insert(active_orbs.begin(), active_orbs.end());
  }
  
  int closed_position = 0;
  int active_position = nclosed_;
  int virt_position = nclosed_ + nactive_;
  const int multisitebasis = hf_ref_->geom()->nbasis();

  auto hf_coeff = hf_ref_->coeff();
  auto out_coeff = hf_coeff->clone();
  assert(hf_coeff->ndim() == multisitebasis);

  for (auto site : active_orbitals) {
    for (int aorb : site)
      copy_n(hf_coeff->element_ptr(0, aorb), multisitebasis, out_coeff->element_ptr(0, active_position++));
  }
  for (int i = 0; i != hf_coeff->mdim(); ++i) {
    if (active_set.count(i) == 0)
      copy_n(hf_coeff->element_ptr(0, i), multisitebasis, out_coeff->element_ptr(0, (closed_position < nclosed_ ? closed_position++ : virt_position++)));
  }
  assert(virt_position == hf_coeff->mdim());

  const int nvirt = multisitebasis - nclosed_ - nactive_;

  sref_ = make_shared<Reference>(hf_ref_->geom(), make_shared<Coeff>(move(*out_coeff)), nclosed_, nactive_, nvirt);
}
