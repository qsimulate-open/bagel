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

shared_ptr<Kramers<1,ZDvec>> ZHarrison::one_down_from_civec(const int nelea, const int neleb, const int istate, shared_ptr<const RelSpace> space1) const {
  auto out = make_shared<Kramers<1,ZDvec>>();

  shared_ptr<const Determinants> int_det = space1->finddet(nelea, neleb);
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
    out->emplace({0}, d);
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
    out->emplace({1}, d);
  }
  return out;
}


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


shared_ptr<Kramers<4,ZDvec>> ZHarrison::four_down_from_civec(const int nelea, const int neleb, const int istate, shared_ptr<const RelSpace> space4) const {
  auto out = make_shared<Kramers<4,ZDvec>>();

  // aaaa
  if (nelea+4 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+4, neleb)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea+2, neleb), norb_*norb_);
    sigma_2e_annih_aa(cc, d);
    auto tmp = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_aa(d->data(i), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j);
    }
    out->emplace({0,0,0,0}, e);
  }
  // aaab
  if (nelea+3 <= norb_ && neleb+1 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+3, neleb+1)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea+2, neleb), norb_*norb_);
    sigma_2e_annih_ab(cc, d);
    auto tmp = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_aa(d->data(i), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j);
    }
    out->emplace({0,0,0,1}, e);
  }
  // aabb
  if (nelea+2 <= norb_ && neleb+2 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+2, neleb+2)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(neleb, nelea+2), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    auto tmp = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_aa(d->data(i)->transpose(), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j);
    }
    out->emplace({0,0,1,1}, e);
  }
  // abbb
  if (nelea+1 <= norb_ && neleb+3 <= norb_) {
#if 0
    shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb+3)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(neleb+1, nelea+1), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    auto tmp = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_ab(d->data(i)->transpose(), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j);
    }
    out->emplace({0,1,1,1}, e);
#else
    shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb+3)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb+2), norb_*norb_);
    sigma_2e_annih_ab(cc, d);
    auto tmp = make_shared<ZDvec>(space4->finddet(neleb, nelea), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_aa(d->data(i)->transpose(), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j)->transpose();
    }
    out->emplace({1,1,0,1}, e);
#endif
  }
  // bbbb
  if (neleb+4 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb+4)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(neleb+2, nelea), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    auto tmp = make_shared<ZDvec>(space4->finddet(neleb, nelea), norb_*norb_);
    auto e = make_shared<ZDvec>(space4->finddet(nelea, neleb), norb_*norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_2e_annih_aa(d->data(i), tmp);
      for (int j = 0; j != norb_*norb_; ++j)
        *e->data(j+i*norb_*norb_) += *tmp->data(j)->transpose();
    }
    out->emplace({1,1,1,1}, e);
  }
  return out;
}


void ZHarrison::compute_rdm34(const int jst, const int ist) {
  auto space4 = make_shared<RelSpace>(norb_, nele_-4);
  auto rdm4 = make_shared<Kramers<8,ZRDM<4>>>();

  // loop over n-4 determinant spaces
  for (int nelea = 0; nelea <= nele_-4; ++nelea) {
    const int neleb = nele_-4 - nelea;
    if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

    shared_ptr<Kramers<4,ZDvec>> interm = four_down_from_civec(nelea, neleb, ist, space4);
    shared_ptr<Kramers<4,ZDvec>> interm2 = ist == jst ? interm : four_down_from_civec(nelea, neleb, jst, space4);
    for (auto& i : *interm2)
      for (auto& j : *interm) {
        auto tmp = make_shared<ZRDM<4>>(norb_);
        auto grouped = group(group(*tmp, 4,8), 0,4);
        contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, grouped, {1,2}, true, false);
        rdm4->add(merge(i.first, j.first), tmp);
      }
  }
  // TODO how can I automate it?
  map<array<int,4>,double> elem;
  elem.emplace(array<int,4>{{0,1,2,3}},  1.0); elem.emplace(array<int,4>{{0,1,3,2}}, -1.0); elem.emplace(array<int,4>{{0,2,1,3}}, -1.0);
  elem.emplace(array<int,4>{{0,2,3,1}},  1.0); elem.emplace(array<int,4>{{0,3,1,2}},  1.0); elem.emplace(array<int,4>{{0,3,2,1}}, -1.0);
  elem.emplace(array<int,4>{{1,0,2,3}}, -1.0); elem.emplace(array<int,4>{{1,0,3,2}},  1.0); elem.emplace(array<int,4>{{1,2,0,3}},  1.0);
  elem.emplace(array<int,4>{{1,2,3,0}}, -1.0); elem.emplace(array<int,4>{{1,3,0,2}}, -1.0); elem.emplace(array<int,4>{{1,3,2,0}},  1.0);
  elem.emplace(array<int,4>{{2,0,1,3}},  1.0); elem.emplace(array<int,4>{{2,0,3,1}}, -1.0); elem.emplace(array<int,4>{{2,1,0,3}}, -1.0);
  elem.emplace(array<int,4>{{2,1,3,0}},  1.0); elem.emplace(array<int,4>{{2,3,0,1}},  1.0); elem.emplace(array<int,4>{{2,3,1,0}}, -1.0);
  elem.emplace(array<int,4>{{3,0,1,2}}, -1.0); elem.emplace(array<int,4>{{3,0,2,1}},  1.0); elem.emplace(array<int,4>{{3,1,0,2}},  1.0);
  elem.emplace(array<int,4>{{3,1,2,0}}, -1.0); elem.emplace(array<int,4>{{3,2,0,1}}, -1.0); elem.emplace(array<int,4>{{3,2,1,0}},  1.0);

  for (auto& i : elem) {
    for (auto& j : elem) {
      array<int,8> perm;
      for (int k = 0; k != 4; ++k) {
        perm[k]   = j.first[k];
        perm[k+4] = i.first[k]+4;
      }
      rdm4->emplace_perm(perm, j.second*i.second);
    }
  }


  auto rdm3 = make_shared<Kramers<6,ZRDM<3>>>();
  for (int i = 0; i != 256; ++i) {
    const bitset<8> source(i);
    if (source[3] == source[7]) {
      shared_ptr<const ZRDM<4>> r4 = rdm4->get_data(KTag<8>(source));
      if (!r4) continue;
      bitset<6> target;
      target[0] = source[0];
      target[1] = source[1];
      target[2] = source[2];
      target[3] = source[4];
      target[4] = source[5];
      target[5] = source[6];

      auto r3 = make_shared<ZRDM<3>>(norb_);
      const size_t norb3 = norb_*norb_*norb_;
      for (int a = 0; a != norb_; ++a)
        for (int j = 0; j != norb3; ++j)
          for (int k = 0; k != norb3; ++k)
            *(r3->data()+k+norb3*j) += *(r4->data()+k+norb3*(a+norb_*(j+norb3*a)));
      rdm3->add(KTag<6>(target), r3);
    }
  }
  for (auto& i : *rdm3)
    i.second->scale(1.0/(nele_-3));
  map<array<int,3>,double> elem3;
  elem3.emplace(array<int,3>{{0,1,2}},  1.0); elem3.emplace(array<int,3>{{0,2,1}}, -1.0); elem3.emplace(array<int,3>{{1,0,2}}, -1.0);
  elem3.emplace(array<int,3>{{1,2,0}},  1.0); elem3.emplace(array<int,3>{{2,0,1}},  1.0); elem3.emplace(array<int,3>{{2,1,0}}, -1.0);
  for (auto& i : elem3) {
    for (auto& j : elem3) {
      array<int,6> perm;
      for (int k = 0; k != 3; ++k) {
        perm[k]   = j.first[k];
        perm[k+3] = i.first[k]+3;
      }
      rdm3->emplace_perm(perm, j.second*i.second);
    }
  }

  auto rdm2 = make_shared<Kramers<4,ZRDM<2>>>();
  for (int i = 0; i != 64; ++i) {
    const bitset<6> source(i);
    if (source[2] == source[5]) {
      shared_ptr<const ZRDM<3>> r3 = rdm3->get_data(KTag<6>(source));
      if (!r3) continue;
      bitset<4> target;
      target[0] = source[0];
      target[1] = source[1];
      target[2] = source[3];
      target[3] = source[4];

      auto r2 = make_shared<ZRDM<2>>(norb_);
      const size_t norb2 = norb_*norb_;
      for (int a = 0; a != norb_; ++a)
        for (int j = 0; j != norb2; ++j)
          for (int k = 0; k != norb2; ++k)
            *(r2->data()+k+norb2*j) += *(r3->data()+k+norb2*(a+norb_*(j+norb2*a)));
      rdm2->add(KTag<4>(target), r2);
    }
  }
  for (auto& i : *rdm2)
    i.second->scale(1.0/(nele_-2));
  rdm2->emplace_perm({{0,1,2,3}},  1.0);
  rdm2->emplace_perm({{0,1,3,2}}, -1.0);
  rdm2->emplace_perm({{1,0,2,3}}, -1.0);
  rdm2->emplace_perm({{1,0,3,2}},  1.0);

  auto rdm1 = make_shared<Kramers<2,ZRDM<1>>>();
  for (int i = 0; i != 16; ++i) {
    const bitset<4> source(i);
    if (source[1] == source[3]) {
      shared_ptr<const ZRDM<2>> r2 = rdm2->get_data(KTag<4>(source));
      if (!r2) continue;
      bitset<2> target;
      target[0] = source[0];
      target[1] = source[2];

      auto r1 = make_shared<ZRDM<1>>(norb_);
      for (int a = 0; a != norb_; ++a)
        for (int j = 0; j != norb_; ++j)
          for (int k = 0; k != norb_; ++k)
            *(r1->data()+k+norb_*j) += *(r2->data()+k+norb_*(a+norb_*(j+norb_*a)));
      rdm1->add(KTag<2>(target), r1);
    }
  }
  for (auto& i : *rdm1)
    i.second->scale(1.0/(nele_-1));

  for (auto& i : *rdm1) {
    cout << i.first.tag()[0] << i.first.tag()[1] << endl;
    i.second->print();
  }

}


void ZHarrison::compute_rdm12() {

  auto space1 = make_shared<RelSpace>(norb_, nele_-1);

  // for one-body RDM
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

      shared_ptr<Kramers<2,ZDvec>> interm = two_down_from_civec(nelea, neleb, istate);
      for (auto& i : *interm)
        for (auto& j : *interm)
          if (i.first >= j.first) {
            auto rdm2 = make_shared<ZRDM<2>>(norb_);
            auto rdm2grouped = group(group(*rdm2, 2,4), 0,2);
            contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, rdm2grouped, {1,2}, true, false);
            rdm2_[istate]->add(merge(i.first, j.first), rdm2);
          }
    }
    // utilize bra-ket symmetry
    for (int i = 0; i != 4; ++i) {
      for (int j = 0; j != 4; ++j) {
        const KTag<4> target((j << 2) + i);
        if (!rdm2_[istate]->exist(target)) {
          const KTag<4> s2301((i << 2) + j);
          if (rdm2_[istate]->exist(s2301)) {
            auto tmp = make_shared<ZRDM<2>>(norb_);
            sort_indices<2,3,0,1,0,1,1,1>(rdm2_[istate]->at(s2301)->data(), tmp->data(), norb_, norb_, norb_, norb_);
            transform(tmp->data(), tmp->data()+norb_*norb_*norb_*norb_, tmp->data(), [](complex<double> a){ return conj(a); });
            rdm2_[istate]->emplace(target, tmp);
          }
        }
      }
    }

    // one body RDM
    // loop over n-1 determinant spaces
    for (int nelea = 0; nelea <= nele_-1; ++nelea) {
      const int neleb = nele_-1 - nelea;
      if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

      shared_ptr<Kramers<1,ZDvec>> interm = one_down_from_civec(nelea, neleb, istate, space1);
      for (auto& i : *interm)
        for (auto& j : *interm)
          if (i.first >= j.first) {
            auto rdm1 = make_shared<ZRDM<1>>(norb_);
            contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, *rdm1, {1,2}, true, false);
            rdm1_[istate]->add(merge(i.first, j.first), rdm1);
          }
    }

    if (!rdm1_[istate]->exist({1,0}))
      rdm1_[istate]->at({1,0}) = rdm1_[istate]->at({0,0})->clone();

    // append permutation information
    rdm2_[istate]->emplace_perm({{0,1,3,2}},-1);
    rdm2_[istate]->emplace_perm({{1,0,3,2}}, 1);
    rdm2_[istate]->emplace_perm({{1,0,2,3}},-1);
  }

  rdm1_av_ = make_shared<Kramers<2,ZRDM<1>>>();
  rdm2_av_ = make_shared<Kramers<4,ZRDM<2>>>();
  if (nstate_ > 1) {
    for (int istate = 0; istate != nstate_; ++istate) {
      for (auto& i : *rdm1_[istate])
        rdm1_av_->add(i.first, i.second);
      for (auto& i : *rdm2_[istate])
        rdm2_av_->add(i.first, i.second);
    }
    for (auto& i : *rdm1_av_) i.second->scale(1.0/nstate_);
    for (auto& i : *rdm2_av_) i.second->scale(1.0/nstate_);
  } else {
    rdm1_av_ = rdm1_.front();
    rdm2_av_ = rdm2_.front();
  }
  rdm2_av_->emplace_perm({{0,1,3,2}},-1);
  rdm2_av_->emplace_perm({{1,0,3,2}}, 1);
  rdm2_av_->emplace_perm({{1,0,2,3}},-1);
}


shared_ptr<const ZMatrix> ZHarrison::rdm1_av() const {
  // RDM transform as D_rs = C*_ri D_ij (C*_rj)^+
  auto out = make_shared<ZMatrix>(norb_*2, norb_*2);
  out->copy_block(    0,     0, norb_, norb_, rdm1_av_kramers("00"));
  out->copy_block(norb_, norb_, norb_, norb_, rdm1_av_kramers("11"));
  out->copy_block(norb_,     0, norb_, norb_, rdm1_av_kramers("10"));
  out->copy_block(    0, norb_, norb_, norb_, out->get_submatrix(norb_, 0, norb_, norb_)->transpose_conjg());
  return out;
}


shared_ptr<const ZMatrix> ZHarrison::rdm2_av() const {
  // transformed 2RDM ; input format is i^+ k^+ j l ; output format is i^+ j k^+ l
  // loop over each component
  auto ikjl = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  for (int i = 0; i != 16; ++i) {
    array<int,4> off {{i/8, (i%8)/4, (i%4)/2, i%2}};
    shared_ptr<const ZRDM<2>> block = rdm2_av_->get_data(i);
    for (int a = 0, abcd = 0; a != norb_; ++a)
      for (int b = 0; b != norb_; ++b)
        for (int c = 0; c != norb_; ++c) {
          copy_n(block->data() + abcd, norb_,
                 ikjl->element_ptr(norb_*off[0]+2*norb_*(c+norb_*off[1]), b+norb_*off[2]+2*norb_*(a+norb_*off[3])));
          abcd += norb_;
        }
  }
  // sort indices : G(ik|jl) -> G(ij|kl)
  auto out  = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  sort_indices<0,2,1,3,0,1,1,1>(ikjl->data(), out->data(), 2*norb_, 2*norb_, 2*norb_, 2*norb_);
  return out;
}
