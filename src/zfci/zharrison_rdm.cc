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
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

void ZHarrison::compute_rdm12() {

  // for one-body RDM
  auto space1 = make_shared<RelSpace>(norb_, nele_-1, 0);
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);

  for (int istate = 0; istate != nstate_; ++istate) {
    // loop over n-2 determinant spaces
    for (int nelea = 0; nelea != nele_-2; ++nelea) {
      const int neleb = nele_-2 - nelea; 
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      // current intermediate determinant
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea, neleb);

      // map 
      unordered_map<bitset<2>, shared_ptr<ZDvec>> intermediates;

      const int ij = norb_*norb_;
      const int nri = int_det->lena() * int_det->lenb(); 

      if (nelea+1 <= norb_ && neleb+1 <= norb_) {
        auto cc = cc_->find(nelea+1, neleb+1)->data(istate);
        auto d = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_ab(cc, d);
        intermediates[bitset<2>("01")] = d; 
      }
      if (neleb+2 <= norb_) {
        // transpose the civec
        auto cc = cc_->find(nelea, neleb+2)->data(istate);
        auto tmp = make_shared<ZDvec>(int_space_->finddet(neleb, nelea), ij);
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

    // one body RDM
    // loop over n-1 determinant spaces
    for (int nelea = 0; nelea != nele_-1; ++nelea) {
      const int neleb = nele_-1 - nelea;
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      shared_ptr<const Determinants> int_det = space1->finddet(nelea, neleb);
      unordered_map<bitset<1>, shared_ptr<ZDvec>> intermediates; 
      const int nri = int_det->lena() * int_det->lenb(); 
  
      if (nelea+1 <= norb_) {
        auto cc = cc_->find(nelea+1, neleb)->data(istate);
        auto d = make_shared<ZDvec>(int_det, norb_);
        const int lenb = int_det->lenb();
        for (auto& s : int_det->stringa()) {
          for (int i = 0; i != norb_; ++i) {
            if (s[i]) continue;
            const double sign = int_det->sign<0>(s, i);
            auto cs = s; cs.set(i);
            complex<double>* target = d->data(i)->data()+int_det->lexical<0>(s)*lenb;
            const complex<double>* source = cc->data()+cc->det()->lexical<0>(cs)*lenb;
            transform(source, source+lenb, target, target, [&sign](complex<double> a, complex<double> b) { return a*sign+b; }); 
          }
        }
        intermediates[bitset<1>("0")] = d;
      }
      if (neleb+1 <= norb_) {
        auto cc = cc_->find(nelea, neleb+1)->data(istate);
        auto d = make_shared<ZDvec>(int_det, norb_);
        const int lenbt = int_det->lenb();
        const int lenbs = cc->det()->lenb();
        for (size_t ia = 0; ia != int_det->lena(); ++ia) {
          const complex<double>* source = cc->data()+ia*lenbs;
          for (int i = 0; i != norb_; ++i) {
            complex<double>* target = d->data(i)->data()+ia*lenbt;
            for (auto& s : int_det->stringb()) {
              if (s[i]) continue;
              const double sign = int_det->sign<1>(s, i);
              auto cs = s; cs.set(i);
              target[int_det->lexical<1>(s)] += sign*source[cc->det()->lexical<1>(cs)];
            }
          }
        }
        intermediates[bitset<1>("1")] = d;
      }
      // almost the same code as above
      for (auto& i : intermediates) {
        for (auto& j : intermediates) {
          if (i.first.to_ulong() <= j.first.to_ulong()) { 
            auto rdm1 = make_shared<ZRDM<1>>(norb_);
            zgemm3m_("c", "n", norb_, norb_, nri, 1.0, i.second->data(0)->data(), nri, j.second->data(0)->data(), nri, 0.0, rdm1->data(), norb_);
            bitset<2> key((i.first.to_ulong() << 1) + j.first.to_ulong());
            if (rdm1_[istate].find(key) == rdm1_[istate].end()) {
              rdm1_[istate][key] = rdm1; 
            } else {
              *rdm1_[istate].at(key) += *rdm1;
            }
          }
        }
      }


    }
  }

  if (nstate_ > 1) {
    for (int istate = 0; istate != nstate_; ++istate) {
      for (auto& i : rdm1_[istate]) {
        if (rdm1_av_.find(i.first) == rdm1_av_.end()) {
          rdm1_av_[i.first] = i.second->copy();
        } else {
          *rdm1_av_.at(i.first) += *i.second;
        }
      }
      for (auto& i : rdm2_[istate]) {
        if (rdm2_av_.find(i.first) == rdm2_av_.end()) {
          rdm2_av_[i.first] = i.second->copy();
        } else {
          *rdm2_av_.at(i.first) += *i.second;
        }
      }
    }
    for (auto& i : rdm1_av_) i.second->scale(1.0/nstate_);
    for (auto& i : rdm2_av_) i.second->scale(1.0/nstate_);
  } else {
    rdm1_av_ = rdm1_.front(); 
    rdm2_av_ = rdm2_.front();
  }


#ifndef NDEBUG
  // Check the FCI energies computed by RDMs and integrals
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
  auto tmp0101 = jop_->mo2e(bitset<4>("0101"))->copy();
  SMITH::sort_indices<1,0,2,3,1,1,-1,1>(jop_->mo2e(bitset<4>("1001"))->data(), tmp0101->data(), norb_, norb_, norb_, norb_);

  cout << "    *  recalculated FCI energy (state averaged)" << endl;
  cout << "           " <<
          zdotc_(norb_*norb_, jop_->mo1e(bitset<2>("00"))->get_conjg()->data(), 1, rdm1_av(bitset<2>("00"))->data(), 1) 
        + zdotc_(norb_*norb_, jop_->mo1e(bitset<2>("01"))->get_conjg()->data(), 1, rdm1_av(bitset<2>("01"))->data(), 1)*2.0 
        + zdotc_(norb_*norb_, jop_->mo1e(bitset<2>("11"))->get_conjg()->data(), 1, rdm1_av(bitset<2>("11"))->data(), 1) 
        + nuc_core 
        + zdotc_(norb_*norb_*norb_*norb_, jop_->mo2e(bitset<4>("0000"))->get_conjg()->data(), 1, rdm2_av(bitset<4>("0000"))->data(), 1)*0.5   // 1
        + zdotc_(norb_*norb_*norb_*norb_, jop_->mo2e(bitset<4>("1111"))->get_conjg()->data(), 1, rdm2_av(bitset<4>("1111"))->data(), 1)*0.5   // 1
        + zdotc_(norb_*norb_*norb_*norb_, tmp0101->get_conjg()->data(), 1, rdm2_av(bitset<4>("0101"))->data(), 1)                             // 4
       << endl;
#endif

}
