//
// BAGEL - Parallel electron correlation program.
// Filename: zharrison_rdm.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/zfci/zharrison.h>

using namespace std;
using namespace bagel;

void ZHarrison::compute_rdm12() {

  // first compute <00|00>
  for (int istate = 0; istate != nstate_; ++istate) {
    // loop over n-2 determinant spaces
    for (int nelea = 0; nelea != nele_-2; ++nelea) {
      const int neleb = nele_-2 - nelea; 
      if (nelea > norb_ || neleb > norb_) continue;

      // current intermediate determinant
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea, neleb);

      // map 
      unordered_map<bitset<2>, shared_ptr<ZDvec>> intermediates;

      const int ij = norb_*norb_;
      const int nri = int_det->lena() * int_det->lenb(); 

      if (nelea+1 <= norb_ || neleb+1 <= norb_) {
        auto cc = cc_->find(nelea+1, neleb+1)->data(istate);
        auto d = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_ab(cc, d);
        intermediates[bitset<2>("01")] = d; 
      }
      if (neleb+2 <= norb_) {
        // transpose the civec
        auto cc = cc_->find(nelea, neleb+2)->data(istate);
        auto tmp = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_aa(cc->transpose(), tmp); 
        // transpose back
        auto d = tmp->clone(); 
        auto jptr = tmp->dvec().begin();
        for (auto& i : d->dvec())
          *i = **jptr++;
        intermediates[bitset<2>("11")] = d; 
      }
      if (nelea+2 <= norb_) {
        // to be implemeneted
        auto cc = cc_->find(nelea+2, neleb)->data(istate);
        auto d = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_aa(cc, d); 
        intermediates[bitset<2>("00")] = d; 
      }

      for (auto& i : intermediates) {
        for (auto& j : intermediates) {
          if (i.first.to_ulong() <= j.first.to_ulong()) { 
            auto rdm2 = make_shared<ZRDM<2>>(norb_);
            zgemm3m_("c", "n", ij, ij, nri, 1.0, i.second->data(0)->data(), nri, j.second->data(0)->data(), nri, 0.0, rdm2->data(), ij);

            bitset<4> key((i.first.to_ulong() << 2) + j.first.to_ulong());

            if (rdm2_[istate].find(key) == rdm2_[istate].end()) {
              rdm2_[istate][key] = rdm2; 
            } else {
              *rdm2_[istate].at(key) += *rdm2;
            }
          }
        }
      }
    
    }
  }
// TODO one-body part has not been implemented
//    zgemv_("c", nri, ij, 1.0, d->data(0)->data(), nri, cc->data(), 1, 0.0, rdm1->data(), 1);

}
