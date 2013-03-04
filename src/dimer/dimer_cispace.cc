//
// BAGEL - Parallel electron correlation program.
// Filename: dimer_cispace.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <utility>

#include <src/dimer/dimer.h>
#include <src/dimer/dimer_cispace.h>

using namespace std;
using namespace bagel;

DimerCISpace::DimerCISpace(const shared_ptr<Dimer> dimer) : dimer_(dimer), norb_(dimer->nact()) {
  // This will probably need to be changed
  nelea_ = make_pair(dimer->nfilledactive().first/2, dimer->nfilledactive().second/2);
  neleb_ = make_pair(dimer->nfilledactive().first/2, dimer->nfilledactive().second/2);
}

// Adds extra Determinants into the mapping that will be needed for Hamiltonian computation
void DimerCISpace::complete() {

  for (auto& imap : cispaceA_) {
    int qa, qb;
    tie(qa, qb) = imap.first;

    auto idet = detspaceA_.find(make_pair(qa-1, qb));
    if (idet == detspaceA_.end()) add_det<0>(qa-1, qb);

    idet = detspaceA_.find(make_pair(qa, qb-1));
    if (idet == detspaceA_.end()) add_det<0>(qa, qb-1);

    idet = detspaceA_.find(make_pair(qa-1, qb-1));
    if (idet == detspaceA_.end()) add_det<0>(qa-1, qb-1);
  }

  for (auto& imap : cispaceB_) {
    int qa, qb;
    tie(qa, qb) = imap.first;

    auto idet = detspaceB_.find(make_pair(qa-1, qb));
    if (idet == detspaceB_.end()) add_det<1>(qa-1, qb);

    idet = detspaceB_.find(make_pair(qa, qb-1));
    if (idet == detspaceB_.end()) add_det<1>(qa, qb-1);

    idet = detspaceB_.find(make_pair(qa-1, qb-1));
    if (idet == detspaceB_.end()) add_det<1>(qa-1, qb-1);
  }
}
