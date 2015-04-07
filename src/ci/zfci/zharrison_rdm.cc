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

#include <src/ci/zfci/zharrison.h>
#include <src/util/prim_op.h>
#include <src/mat1e/rel/reloverlap.h>

using namespace std;
using namespace bagel;
using namespace btas;

shared_ptr<Kramers<2,ZDvec>> ZHarrison::two_down_from_civec(const int nelea, const int neleb, const int istate) const {
  auto out = make_shared<Kramers<2,ZDvec>>();

  if (nelea+1 <= norb_ && neleb+1 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb+1)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb), norb_*norb_);
    sigma_2e_annih_ab(cc, d);
    out->emplace({0,1}, d);
  }
  if (neleb+2 <= norb_) {
    // transpose the civec
    shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb+2)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(neleb, nelea), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    // transpose back
    for (auto& i : d->dvec())
      *i = *i->transpose();
    out->emplace({1,1}, d);
  }
  if (nelea+2 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+2, neleb)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb), norb_*norb_);
    sigma_2e_annih_aa(cc, d);
    out->emplace({0,0}, d);
  }
  return out;
}


void ZHarrison::compute_rdm12() {

  // for one-body RDM
  auto space1 = make_shared<RelSpace>(norb_, nele_-1);
  rdm1_.clear();
  rdm2_.clear();
  for (int istate = 0; istate != nstate_; ++istate)  {
    rdm1_.push_back(make_shared<Kramers<2,ZRDM<1>>>());
    rdm2_.push_back(make_shared<Kramers<4,ZRDM<2>>>());
  }

  for (int istate = 0; istate != nstate_; ++istate) {
    // loop over n-2 determinant spaces
    for (int nelea = 0; nelea <= nele_-2; ++nelea) {
      const int neleb = nele_-2 - nelea;
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      // map
      shared_ptr<Kramers<2,ZDvec>> interm = two_down_from_civec(nelea, neleb, istate);

      for (auto& i : *interm) {
        for (auto& j : *interm) {
          if (i.first >= j.first) {
            auto rdm2 = make_shared<ZRDM<2>>(norb_);
            auto rdm2grouped = group(group(*rdm2, 2,4), 0,2);
            contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, rdm2grouped, {1,2}, true, false);
            rdm2_[istate]->add(merge(i.first, j.first), rdm2);
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
      Kramers<1,ZDvec> intermediates;

      if (nelea+1 <= norb_) {
        shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb)->data(istate);
        auto d = make_shared<ZDvec>(int_det, norb_);
        const int lenb = int_det->lenb();
        for (auto& s : int_det->string_bits_a()) {
          for (int i = 0; i != norb_; ++i) {
            if (s[i]) continue;
            const double sign = int_det->sign<0>(s, i);
            auto cs = s; cs.set(i);
            complex<double>* target = d->data(i)->data()+int_det->lexical<0>(s)*lenb;
            const complex<double>* source = cc->data()+cc->det()->lexical<0>(cs)*lenb;
            blas::ax_plus_y_n(sign, source, lenb, target);
          }
        }
        intermediates.emplace({0}, d);
      }
      if (neleb+1 <= norb_) {
        shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb+1)->data(istate);
        auto d = make_shared<ZDvec>(int_det, norb_);
        const int lenbt = int_det->lenb();
        const int lenbs = cc->det()->lenb();
        for (size_t ia = 0; ia != int_det->lena(); ++ia) {
          const complex<double>* source = cc->data()+ia*lenbs;
          for (int i = 0; i != norb_; ++i) {
            complex<double>* target = d->data(i)->data()+ia*lenbt;
            for (auto& s : int_det->string_bits_b()) {
              if (s[i]) continue;
              const double sign = int_det->sign<1>(s, i);
              auto cs = s; cs.set(i);
              target[int_det->lexical<1>(s)] += sign*source[cc->det()->lexical<1>(cs)];
            }
          }
        }
        intermediates.emplace({1}, d);
      }
      // almost the same code as above
      for (auto& i : intermediates) {
        for (auto& j : intermediates) {
          if (i.first >= j.first) {
            auto rdm1 = make_shared<ZRDM<1>>(norb_);
            contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, *rdm1, {1,2}, true, false);
            rdm1_[istate]->add(merge(i.first, j.first), rdm1);
          }
        }
      }
    }

    // for completeness we can compute all the blocks (2RDM), which are useful in CASSCF
    if (!rdm1_[istate]->exist({1,0}))
      rdm1_[istate]->at({1,0}) = rdm1_[istate]->at({0,0})->clone();

    for (int i = 0; i != 4; ++i) {
      for (int j = 0; j != 4; ++j) {
        bitset<4> target((j << 2) + i);

        if (!rdm2_[istate]->exist(target)) {
          bitset<4> s2301((i << 2) + j);
          bitset<4> s1032; s1032.set(0, target[1]); s1032.set(1, target[0]); s1032.set(2, target[3]); s1032.set(3, target[2]);
          bitset<4> s3210; s3210.set(0, target[3]); s3210.set(1, target[2]); s3210.set(2, target[1]); s3210.set(3, target[0]);
          bitset<4> s0132; s0132.set(0, target[0]); s0132.set(1, target[1]); s0132.set(2, target[3]); s0132.set(3, target[2]);

          // TODO I don't know why I cannot move emplace line outside if block (if I do, assert fails)
          if (rdm2_[istate]->exist(s2301)) {
            rdm2_[istate]->emplace(target, make_shared<ZRDM<2>>(norb_));
            sort_indices<2,3,0,1,0,1,1,1>(rdm2_[istate]->at(s2301)->data(), rdm2_[istate]->at(target)->data(), norb_, norb_, norb_, norb_);
            transform(rdm2_[istate]->at(target)->data(), rdm2_[istate]->at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate]->at(target)->data(), [](complex<double> a){ return conj(a); });
          } else if (rdm2_[istate]->exist(s1032)) {
            rdm2_[istate]->emplace(target, make_shared<ZRDM<2>>(norb_));
            sort_indices<1,0,3,2,0,1,1,1>(rdm2_[istate]->at(s1032)->data(), rdm2_[istate]->at(target)->data(), norb_, norb_, norb_, norb_);
          } else if (rdm2_[istate]->exist(s3210)) {
            rdm2_[istate]->emplace(target, make_shared<ZRDM<2>>(norb_));
            sort_indices<3,2,1,0,0,1,1,1>(rdm2_[istate]->at(s3210)->data(), rdm2_[istate]->at(target)->data(), norb_, norb_, norb_, norb_);
            transform(rdm2_[istate]->at(target)->data(), rdm2_[istate]->at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate]->at(target)->data(), [](complex<double> a){ return conj(a); });
          } else if (rdm2_[istate]->exist(s0132)) {
            rdm2_[istate]->emplace(target, make_shared<ZRDM<2>>(norb_));
            sort_indices<1,0,2,3,0,1,1,1>(rdm2_[istate]->at(s0132)->data(), rdm2_[istate]->at(target)->data(), norb_, norb_, norb_, norb_);
            transform(rdm2_[istate]->at(target)->data(), rdm2_[istate]->at(target)->data()+norb_*norb_*norb_*norb_, rdm2_[istate]->at(target)->data(), [](complex<double> a){ return -a; });
          } else {
            rdm2_[istate]->emplace(target, make_shared<ZRDM<2>>(norb_));
          }
        }
      }
    }
  }

  rdm1_av_ = make_shared<Kramers<2,ZRDM<1>>>();
  rdm2_av_ = make_shared<Kramers<4,ZRDM<2>>>();
  if (nstate_ > 1) {
    for (int istate = 0; istate != nstate_; ++istate) {
      for (auto& i : *rdm1_[istate])
        rdm1_av_->add(i.first, i.second->copy());
      for (auto& i : *rdm2_[istate])
        rdm2_av_->add(i.first, i.second->copy());
    }
    for (auto& i : *rdm1_av_) i.second->scale(1.0/nstate_);
    for (auto& i : *rdm2_av_) i.second->scale(1.0/nstate_);
  } else {
    rdm1_av_ = rdm1_.front();
    rdm2_av_ = rdm2_.front();
  }
}


shared_ptr<const ZMatrix> ZHarrison::rdm1_av() const {
  // RDM transform as D_rs = C*_ri D_ij (C*_rj)^+
  auto rdm1_tot = make_shared<ZMatrix>(norb_*2, norb_*2);
  rdm1_tot->copy_block(    0,     0, norb_, norb_, rdm1_av_kramers("00"));
  rdm1_tot->copy_block(norb_, norb_, norb_, norb_, rdm1_av_kramers("11"));
  rdm1_tot->copy_block(norb_,     0, norb_, norb_, rdm1_av_kramers("10"));
  rdm1_tot->copy_block(    0, norb_, norb_, norb_, rdm1_tot->get_submatrix(norb_, 0, norb_, norb_)->transpose_conjg());
  return rdm1_tot;
}


shared_ptr<const ZMatrix> ZHarrison::rdm2_av() const {
  // transformed 2RDM ; input format is i^+ k^+ j l ; output format is i^+ j k^+ l
  // TODO : slot in 2RDM by hand to avoid matrix multiplication ; should be slightly cheaper

  unordered_map<bitset<1>, shared_ptr<const ZMatrix>> trans;
  auto unit = make_shared<ZMatrix>(norb_*2,norb_*2);
  unit->unit();
  for (int i = 0; i != 2; ++i) {
    trans.emplace(i, unit->slice_copy(i*norb_,(i+1)*norb_));
  }

  // loop over each component
  auto ikjl = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  auto out  = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  for (auto& irdm : *rdm2_av_) {
    bitset<4> ib = irdm.first.tag();
    shared_ptr<const ZRDM<2>> rdm2 = irdm.second;
    // TODO to be updated once the Tensor branch comes out
    const int norb2 = norb_*norb_;
    const int norb3 = norb2*norb_;
    const int norb4 = norb3*norb_;
    unique_ptr<complex<double>[]> tmp1(new complex<double>[2*norb4]);
    unique_ptr<complex<double>[]> tmp2(new complex<double>[4*norb4]);
    unique_ptr<complex<double>[]> tmp3(new complex<double>[8*norb4]);
    zgemm3m_("N", "N", 2*norb_, norb3, norb_, 1.0, trans[ib[3] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, rdm2->data(), norb_, 0.0, tmp1.get(), 2*norb_);
    for (int i = 0; i != norb2; ++i)
      zgemm3m_("N", "T", 2*norb_, 2*norb_, norb_, 1.0, tmp1.get()+i*2*norb2, 2*norb_, trans[ib[2] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 0.0, tmp2.get()+i*4*norb2, 2*norb_);
    for (int i = 0; i != norb_; ++i)
      zgemm3m_("N", "C", 4*norb2, 2*norb_, norb_, 1.0, tmp2.get()+i*4*norb3, 4*norb2, trans[ib[1] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 0.0, tmp3.get()+i*8*norb3, 4*norb2);
    zgemm3m_("N", "C", 8*norb3, 2*norb_, norb_, 1.0, tmp3.get(),           8*norb3, trans[ib[0] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 1.0, ikjl->data()         , 8*norb3);
  }

  // sort indices : G(ik|jl) -> G(ij|kl)
  sort_indices<0,2,1,3,0,1,1,1>(ikjl->data(), out->data(), 2*norb_, 2*norb_, 2*norb_, 2*norb_);
  return out;
}


shared_ptr<const ZMatrix> ZHarrison::mo2e_full() const {
  // transformed two-electron energies ; input format is i^+ k^+ j l ; output format is i^+ j k^+ l

  unordered_map<bitset<1>, shared_ptr<const ZMatrix>> trans;
  auto unit = make_shared<ZMatrix>(norb_*2,norb_*2);
  unit->unit();
  for (int i = 0; i != 2; ++i) {
    auto co = make_shared<const ZMatrix>(unit->slice(i*norb_,(i+1)*norb_));
    bitset<1> b(i);
    trans.insert(make_pair(b, co));
  }

  // loop over each component
  auto ikjl = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  auto out  = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  for (auto& irdm : *jop()->mo2e()) {
    bitset<4> ib = irdm.first.tag();
    shared_ptr<const ZMatrix> rdm2 = irdm.second;
    // TODO to be updated once the Tensor branch comes out
    const int norb2 = norb_*norb_;
    const int norb3 = norb2*norb_;
    const int norb4 = norb3*norb_;
    unique_ptr<complex<double>[]> tmp1(new complex<double>[2*norb4]);
    unique_ptr<complex<double>[]> tmp2(new complex<double>[4*norb4]);
    unique_ptr<complex<double>[]> tmp3(new complex<double>[8*norb4]);
    zgemm3m_("N", "N", 2*norb_, norb3, norb_, 1.0, trans[ib[3] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, rdm2->data(), norb_, 0.0, tmp1.get(), 2*norb_);
    for (int i = 0; i != norb2; ++i)
      zgemm3m_("N", "T", 2*norb_, 2*norb_, norb_, 1.0, tmp1.get()+i*2*norb2, 2*norb_, trans[ib[2] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 0.0, tmp2.get()+i*4*norb2, 2*norb_);
    for (int i = 0; i != norb_; ++i)
      zgemm3m_("N", "C", 4*norb2, 2*norb_, norb_, 1.0, tmp2.get()+i*4*norb3, 4*norb2, trans[ib[1] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 0.0, tmp3.get()+i*8*norb3, 4*norb2);
    zgemm3m_("N", "C", 8*norb3, 2*norb_, norb_, 1.0, tmp3.get(),           8*norb3, trans[ib[0] ? bitset<1>(1) : bitset<1>(0)]->data(), 2*norb_, 1.0, ikjl->data()         , 8*norb3);
  }

  // sort indices : G(ik|jl) -> G(ij|kl)
  sort_indices<0,2,1,3,0,1,1,1>(ikjl->data(), out->data(), 2*norb_, 2*norb_, 2*norb_, 2*norb_);
  return out;
}
