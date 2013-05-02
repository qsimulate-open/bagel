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

#include <src/dimer/dimer_cispace.h>

using namespace std;
using namespace bagel;

// Completes the spin case and adds extra Determinants into the mapping that will be needed for Hamiltonian computation
void DimerCISpace::complete() {
  {
    vector<SpaceKey> references;
    for (auto& imap : cispaceA_) { references.push_back(imap.first); }
    // These spaces are assumed to be high-spin
    for (auto& ispace : references) {
      if (ispace.S > 0) {
        const int S = ispace.S;

        shared_ptr<Dvec> ref_state = ccvec<0>(ispace);

        int ref_qa, ref_qb;
        std::tie(ref_qa, ref_qb) = detkey<0>(ispace);
        const int mult = S + 1;

        for (int i = 1; i < mult; ++i) {
          const int nqa = ref_qa - i;
          const int nqb = ref_qb + i;

          shared_ptr<Determinants> det = add_det<0>(nqa, nqb);

          ref_state = ref_state->spin_lower(det);
          for (int istate = 0; istate < ref_state->ij(); ++istate) {
            const double norm = ref_state->data(istate)->norm();
            if ( norm < numerical_zero__ ) throw std::runtime_error("Spin lowering operator yielded no state.");
            ref_state->data(istate)->scale(1.0/norm);
          }
          insert<0>(ref_state, S);
        }
      }
    }
  }
  

  {
    vector<SpaceKey> references;
    for (auto& imap : cispaceB_) { references.push_back(imap.first); }
    // These spaces are assumed to be high-spin
    for (auto& ispace : references) {
      if (ispace.S > 0) {
        const int S = ispace.S;

        shared_ptr<Dvec> ref_state = ccvec<1>(ispace);

        int ref_qa, ref_qb;
        std::tie(ref_qa, ref_qb) = detkey<1>(ispace);
        const int mult = S + 1;

        for (int i = 1; i < mult; ++i) {
          const int nqa = ref_qa - i;
          const int nqb = ref_qb + i;

          shared_ptr<Determinants> det = add_det<1>(nqa, nqb);

          ref_state = ref_state->spin_lower(det);
          for (int istate = 0; istate < ref_state->ij(); ++istate) {
            const double norm = ref_state->data(istate)->norm();
            if ( norm < numerical_zero__ ) throw std::runtime_error("Spin lowering operator yielded no state.");
            ref_state->data(istate)->scale(1.0/norm);
          }
          insert<1>(ref_state, S);
        }
      }
    }
  }

  for (auto& imap : cispaceA_) {
    int qa, qb;
    tie(qa, qb) = detkey<0>(imap.first);

    add_det<0>(qa-1, qb);
    add_det<0>(qa, qb-1);
    add_det<0>(qa-1, qb-1);
  }

  for (auto& imap : cispaceB_) {
    int qa, qb;
    tie(qa, qb) = detkey<1>(imap.first);

    add_det<1>(qa-1, qb);
    add_det<1>(qa, qb-1);
    add_det<1>(qa-1, qb-1);
  }
}
