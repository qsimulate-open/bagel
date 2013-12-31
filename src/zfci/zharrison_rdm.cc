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
  rdm1_.clear();
  rdm2_.clear();
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);

  for (int istate = 0; istate != nstate_; ++istate) {
    // loop over n-2 determinant spaces
    for (int nelea = 0; nelea <= nele_-2; ++nelea) {
      const int neleb = nele_-2 - nelea;
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      // current intermediate determinant
      shared_ptr<const Determinants> int_det = int_space_->finddet(nelea, neleb);

      // map
      map<string, shared_ptr<ZDvec>> intermediates;

      const int ij = norb_*norb_;
      const int nri = int_det->lena() * int_det->lenb();

      if (nelea+1 <= norb_ && neleb+1 <= norb_) {
        auto cc = cc_->find(nelea+1, neleb+1)->data(istate);
        auto d = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_ab(cc, d);
        intermediates["01"] = d;
      }
      if (neleb+2 <= norb_) {
        // transpose the civec
        auto cc = cc_->find(nelea, neleb+2)->data(istate);
        auto tmp = make_shared<ZDvec>(int_space_->finddet(neleb, nelea), ij);
        sigma_2e_annih_aa(cc->transpose(), tmp);
        // transpose back
        auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb), ij);
        auto jptr = tmp->dvec().begin();
        for (auto& i : d->dvec()) {
          *i = *(*jptr)->transpose();
          ++jptr;
        }
        intermediates["11"] = d;
      }
      if (nelea+2 <= norb_) {
        // to be implemeneted
        auto cc = cc_->find(nelea+2, neleb)->data(istate);
        auto d = make_shared<ZDvec>(int_det, ij);
        sigma_2e_annih_aa(cc, d);
        intermediates["00"] = d;
      }

      for (auto& i : intermediates) {
        for (auto& j : intermediates) {
          if (i.first >= j.first) {
            auto rdm2 = make_shared<ZRDM<2>>(norb_);
            zgemm3m_("c", "n", ij, ij, nri, 1.0, i.second->data(0)->data(), nri, j.second->data(0)->data(), nri, 0.0, rdm2->data(), ij);

            bitset<4> key(i.first + j.first);

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
    for (int nelea = 0; nelea <= nele_-1; ++nelea) {
      const int neleb = nele_-1 - nelea;
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      shared_ptr<const Determinants> int_det = space1->finddet(nelea, neleb);
      map<string, shared_ptr<ZDvec>> intermediates;
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
            blas::ax_plus_y_n(sign, source, lenb, target);
          }
        }
        intermediates["0"] = d;
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
        intermediates["1"] = d;
      }
      // almost the same code as above
      for (auto& i : intermediates) {
        for (auto& j : intermediates) {
          if (i.first >= j.first) {
            auto rdm1 = make_shared<ZRDM<1>>(norb_);
            zgemm3m_("c", "n", norb_, norb_, nri, 1.0, i.second->data(0)->data(), nri, j.second->data(0)->data(), nri, 0.0, rdm1->data(), norb_);
            bitset<2> key(i.first + j.first);
            if (rdm1_[istate].find(key) == rdm1_[istate].end()) {
              rdm1_[istate][key] = rdm1;
            } else {
              *rdm1_[istate].at(key) += *rdm1;
            }
          }
        }
      }
    }

    // for completeness we can compute all the blocks (2RDM), which are useful in CASSCF
    if (rdm1_[istate].find(bitset<2>("10")) == rdm1_[istate].end())
      rdm1_[istate][bitset<2>("10")] = rdm1_[istate].at(bitset<2>("00"))->clone();

    for (int i = 0; i != 4; ++i) {
      for (int j = 0; j != 4; ++j) {
        bitset<4> target((j << 2) + i);
        bitset<4> s2301((i << 2) + j);
        bitset<4> s1032; s1032.set(0, target[1]); s1032.set(1, target[0]); s1032.set(2, target[3]); s1032.set(3, target[2]);
        bitset<4> s3210; s3210.set(0, target[3]); s3210.set(1, target[2]); s3210.set(2, target[1]); s3210.set(3, target[0]);
        bitset<4> s0132; s0132.set(0, target[0]); s0132.set(1, target[1]); s0132.set(2, target[3]); s0132.set(3, target[2]);
        if (rdm2_[istate].find(target) != rdm2_[istate].end()) {
          continue;
        } else if (rdm2_[istate].find(s2301) != rdm2_[istate].end()) {
          rdm2_[istate][target] = make_shared<ZRDM<2>>(norb_);
          SMITH::sort_indices<2,3,0,1,0,1,1,1>(rdm2_[istate].at(s2301)->data(), rdm2_[istate].at(target)->data(), norb_, norb_, norb_, norb_);
          transform(rdm2_[istate].at(target)->data(), rdm2_[istate].at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate].at(target)->data(), [](complex<double> a){ return conj(a); });
        } else if (rdm2_[istate].find(s1032) != rdm2_[istate].end()) {
          rdm2_[istate][target] = make_shared<ZRDM<2>>(norb_);
          SMITH::sort_indices<1,0,3,2,0,1,1,1>(rdm2_[istate].at(s1032)->data(), rdm2_[istate].at(target)->data(), norb_, norb_, norb_, norb_);
        } else if (rdm2_[istate].find(s3210) != rdm2_[istate].end()) {
          rdm2_[istate][target] = make_shared<ZRDM<2>>(norb_);
          SMITH::sort_indices<3,2,1,0,0,1,1,1>(rdm2_[istate].at(s3210)->data(), rdm2_[istate].at(target)->data(), norb_, norb_, norb_, norb_);
          transform(rdm2_[istate].at(target)->data(), rdm2_[istate].at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate].at(target)->data(), [](complex<double> a){ return conj(a); });
        } else if (rdm2_[istate].find(s0132) != rdm2_[istate].end()) {
          rdm2_[istate][target] = make_shared<ZRDM<2>>(norb_);
          SMITH::sort_indices<1,0,2,3,0,1,1,1>(rdm2_[istate].at(s0132)->data(), rdm2_[istate].at(target)->data(), norb_, norb_, norb_, norb_);
          transform(rdm2_[istate].at(target)->data(), rdm2_[istate].at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate].at(target)->data(), [](complex<double> a){ return -a; });
        } else {
          // This is dangerous.. TODO
          rdm2_[istate][target] = make_shared<ZRDM<2>>(norb_);
#if 0
          cout << target << endl;
          for (auto& i : rdm2_[istate]) cout << i.first << endl;
          throw logic_error("debug .. ZHarrison::compute_rdm12()");
#endif
        }
      }
    }
  }

  rdm1_av_.clear();
  rdm2_av_.clear();
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

#if 0
  // Check the FCI energies computed by RDMs and integrals
  const double nuc_core = geom_->nuclear_repulsion() + jop_->core_energy();
  auto tmp0101 = jop_->mo2e(bitset<4>("0101"))->copy();
  SMITH::sort_indices<1,0,2,3,1,1,-1,1>(jop_->mo2e(bitset<4>("1001"))->data(), tmp0101->data(), norb_, norb_, norb_, norb_);

  const int n = norb_;
  auto trace1 = [this](const string st) {
    shared_ptr<const ZMatrix> a = jop_->mo1e(bitset<2>(st));
    shared_ptr<const ZRDM<1>> b = rdm1_av(bitset<2>(st));
    return inner_product(a->data(), a->data()+a->size(), b->data(), complex<double>(0.0),
                         std::plus<complex<double>>(), [](const complex<double>& i, const complex<double>& j){ return i*j; });
  };
  auto trace2 = [this](const string st) {
    shared_ptr<const ZMatrix> a = jop_->mo2e(bitset<4>(st));
    shared_ptr<const ZRDM<2>> b = rdm2_av(bitset<4>(st));
    return inner_product(a->data(), a->data()+a->size(), b->data(), complex<double>(0.0),
                         std::plus<complex<double>>(), [](const complex<double>& i, const complex<double>& j){ return i*j; });
  };

  const double recomp_energy =
          (trace1("00") + trace1("10")*2.0 + trace1("11")).real() + nuc_core
        + (trace2("0000")*0.5 + trace2("1111")*0.5
            + zdotc_(norb_*norb_*norb_*norb_, tmp0101->get_conjg()->data(), 1, rdm2_av(bitset<4>("0101"))->data(), 1)
            + trace2("1100") + trace2("0100")*2.0 + trace2("1101")*2.0).real();
  cout << "    *  recalculated FCI energy (state averaged)" << endl;
  cout << setw(29) << setprecision(8) << recomp_energy << endl;

  // checking against the original energies
  const double orig_energy = accumulate(energy_.begin(), energy_.end(), 0.0) / energy_.size();
  assert(fabs(orig_energy - recomp_energy) < 1.0e-8);

#if 0
  complex<double> recomp_energy2 = (trace1("00") + trace1("10")*2.0 + trace1("11")).real() + nuc_core;
  for (int i = 0; i != 16; ++i) {
    stringstream ss; ss << bitset<4>(i);
    recomp_energy2 += trace2(ss.str())*0.5;
  }
  cout << "    *  recalculated FCI energy (state averaged)" << endl;
  cout << setw(49) << setprecision(12) << recomp_energy2 << endl;
  assert(fabs(orig_energy - recomp_energy2.real()) < 1.0e-8);
#endif
#endif

}
