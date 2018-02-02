//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zharrison_rdm.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
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

#include <src/ci/zfci/zharrison.h>
#include <src/util/prim_op.h>
#include <src/mat1e/rel/reloverlap.h>

using namespace std;
using namespace bagel;
using namespace btas;

ZFCI_bare::ZFCI_bare(shared_ptr<const RelCIWfn> ci) {
  print_thresh_ = 1.0e-8;
  nele_ = ci->det()->first->nele();
  ncore_ = ci->ncore();
  norb_  = ci->nact();
  nstate_ = ci->nstates();
  energy_ = ci->energies();
  cc_ = ci->civectors() ? ci->civectors()->copy() : nullptr;
  space_ = ci->det()->first;
  int_space_ = ci->det()->second;
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);
}


void ZHarrison::sigma_1e_annih_a(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  shared_ptr<const Determinants> int_det = d->det();
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
}


void ZHarrison::sigma_1e_annih_b(shared_ptr<const ZCivec> cc, shared_ptr<ZDvec> d) const {
  shared_ptr<const Determinants> int_det = d->det();
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
}


shared_ptr<Kramers<1,ZDvec>> ZHarrison::one_down_from_civec(const int nelea, const int neleb, const int istate, shared_ptr<const RelSpace> space1) const {
  auto out = make_shared<Kramers<1,ZDvec>>();

  shared_ptr<const Determinants> int_det = space1->finddet(nelea, neleb);
  if (nelea+1 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb)->data(istate);
    auto d = make_shared<ZDvec>(int_det, norb_);
    sigma_1e_annih_a(cc, d);
    out->emplace({0}, d);
  }
  if (neleb+1 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb+1)->data(istate);
    auto d = make_shared<ZDvec>(int_det, norb_);
    sigma_1e_annih_b(cc, d);
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
    auto d2 = make_shared<ZDvec>(int_space_->finddet(nelea, neleb), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    // transpose back
    auto iter = d2->dvec().begin();
    for (auto& i : d->dvec())
      **iter++ = *i->transpose();
    out->emplace({1,1}, d2);
  }
  if (nelea+2 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+2, neleb)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb), norb_*norb_);
    sigma_2e_annih_aa(cc, d);
    out->emplace({0,0}, d);
  }
  return out;
}


shared_ptr<Kramers<3,ZDvec>> ZHarrison::three_down_from_civec(const int nelea, const int neleb, const int istate, shared_ptr<const RelSpace> space3) const {
  auto out = make_shared<Kramers<3,ZDvec>>();

  // aaa
  if (nelea+3 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+3, neleb)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea+1, neleb), norb_*norb_);
    sigma_2e_annih_aa(cc, d);
    auto tmp = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_);
    auto e = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_1e_annih_a(d->data(i), tmp);
      for (int j = 0; j != norb_; ++j)
        *e->data(j+i*norb_) += *tmp->data(j);
    }
    out->emplace({0,0,0}, e);
  }
  // aab
  if (nelea+2 <= norb_ && neleb+1 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+2, neleb+1)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea+1, neleb), norb_*norb_);
    sigma_2e_annih_ab(cc, d);
    auto tmp = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_);
    auto e = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_1e_annih_a(d->data(i), tmp);
      for (int j = 0; j != norb_; ++j)
        *e->data(j+i*norb_) += *tmp->data(j);
    }
    out->emplace({0,0,1}, e);
  }
  // bab
  if (nelea+1 <= norb_ && neleb+2 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea+1, neleb+2)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(nelea, neleb+1), norb_*norb_);
    sigma_2e_annih_ab(cc, d);
    auto tmp = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_);
    auto e = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_1e_annih_b(d->data(i), tmp);
      for (int j = 0; j != norb_; ++j)
        *e->data(j+i*norb_) += *tmp->data(j);
    }
    out->emplace({1,0,1}, e);
  }
  // bbb
  if (neleb+3 <= norb_) {
    shared_ptr<const ZCivec> cc = cc_->find(nelea, neleb+3)->data(istate);
    auto d = make_shared<ZDvec>(int_space_->finddet(neleb+1, nelea), norb_*norb_);
    sigma_2e_annih_aa(cc->transpose(), d);
    auto tmp = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_);
    auto e = make_shared<ZDvec>(space3->finddet(nelea, neleb), norb_*norb_*norb_);
    for (int i = 0; i != norb_*norb_; ++i) {
      tmp->zero();
      sigma_1e_annih_b(d->data(i)->transpose(), tmp);
      for (int j = 0; j != norb_; ++j)
        *e->data(j+i*norb_) += *tmp->data(j);
    }
    out->emplace({1,1,1}, e);
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


shared_ptr<Kramers<8,ZRDM<4>>> ZHarrison::rdm4(const int jst, const int ist) const {
  // loop over n-4 determinant spaces
  auto rdm4 = make_shared<Kramers<8,ZRDM<4>>>();
  if (nele_ < 4) return rdm4;

  auto space4 = make_shared<RelSpace>(norb_, nele_-4);
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

        auto sorted = tmp->clone();
        sort_indices<0,4,1,5,2,6,3,7,0,1,1,1>(tmp->data(), sorted->data(), norb_, norb_, norb_, norb_, norb_, norb_, norb_, norb_);
        auto ijtag = merge(i.first, j.first);
        rdm4->add(ijtag.perm({{0,4,1,5,2,6,3,7}}), sorted);
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
      vector<int> perm(8);
      for (int k = 0; k != 4; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      rdm4->emplace_perm(perm, j.second*i.second);
    }
  }
  return rdm4;
}


tuple<shared_ptr<Kramers<6,ZRDM<3>>>,shared_ptr<Kramers<6,ZRDM<3>>>> ZHarrison::rdm34f(const int jst, const int ist, shared_ptr<const ZMatrix> fockact) const {
  shared_ptr<Kramers<6,ZRDM<3>>> rdm3t = rdm3(jst, ist);
  shared_ptr<Kramers<8,ZRDM<4>>> rdm4t = rdm4(jst, ist);
  auto rdm4f = make_shared<Kramers<6,ZRDM<3>>>();
  const int n = fockact->ndim()/2;

  Kramers<2,ZMatrix> fock;
  fock.emplace(0, fockact->get_submatrix(0, 0, n, n));
  fock.emplace(1, fockact->get_submatrix(n, 0, n, n));
  fock.emplace(2, fockact->get_submatrix(0, n, n, n));
  fock.emplace(3, fockact->get_submatrix(n, n, n, n));
  for (int i = 0; i != 64; ++i) {
    for (int j = 0; j != 4; ++j) {
      auto work = make_shared<ZRDM<3>>(n);
      shared_ptr<const ZMatrix> cfock = fock.at(j);
      shared_ptr<const ZRDM<4>> crdm = rdm4t->get_data(i * 4 + j);
      if (!crdm) continue;

      auto wgr = btas::group(*work, 0,6);
      auto crdmgr = btas::group(btas::group(*crdm, 6,8),0,6);
      auto fgr = btas::group(*cfock, 0,2);
      btas::contract(1.0, crdmgr, {0,1}, fgr, {1}, 0.0, wgr, {0});
      rdm4f->add(i, work);
    }
  }

  return make_tuple(rdm3t, rdm4f);
}


shared_ptr<Kramers<6,ZRDM<3>>> ZHarrison::rdm3(const int jst, const int ist) const {
  // loop over n-3 determinant spaces
  auto rdm3 = make_shared<Kramers<6,ZRDM<3>>>();
  if (nele_ < 3) return rdm3;

  auto space3 = make_shared<RelSpace>(norb_, nele_-3);
  for (int nelea = 0; nelea <= nele_-3; ++nelea) {
    const int neleb = nele_-3 - nelea;
    if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

    shared_ptr<Kramers<3,ZDvec>> interm = three_down_from_civec(nelea, neleb, ist, space3);
    shared_ptr<Kramers<3,ZDvec>> interm2 = ist == jst ? interm : three_down_from_civec(nelea, neleb, jst, space3);
    for (auto& i : *interm2)
      for (auto& j : *interm) {
        auto tmp = make_shared<ZRDM<3>>(norb_);
        auto grouped = group(group(*tmp, 3,6), 0,3);
        contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, grouped, {1,2}, true, false);

        auto sorted = tmp->clone();
        sort_indices<0,3,1,4,2,5,0,1,1,1>(tmp->data(), sorted->data(), norb_, norb_, norb_, norb_, norb_, norb_);
        auto ijtag = merge(i.first, j.first);
        rdm3->add(ijtag.perm({{0,3,1,4,2,5}}), sorted);
      }
  }

  map<array<int,3>,double> elem3;
  elem3.emplace(array<int,3>{{0,1,2}},  1.0); elem3.emplace(array<int,3>{{0,2,1}}, -1.0); elem3.emplace(array<int,3>{{1,0,2}}, -1.0);
  elem3.emplace(array<int,3>{{1,2,0}},  1.0); elem3.emplace(array<int,3>{{2,0,1}},  1.0); elem3.emplace(array<int,3>{{2,1,0}}, -1.0);
  for (auto& i : elem3) {
    for (auto& j : elem3) {
      vector<int> perm(6);
      for (int k = 0; k != 3; ++k) {
        perm[k*2]   = j.first[k]*2;
        perm[k*2+1] = i.first[k]*2+1;
      }
      rdm3->emplace_perm(perm, j.second*i.second);
    }
  }
  return rdm3;
}


shared_ptr<Kramers<4,ZRDM<2>>> ZHarrison::rdm2(const int jst, const int ist) const {
  auto out = make_shared<Kramers<4,ZRDM<2>>>();
  if (nele_ < 2) return out;

  // loop over n-2 determinant spaces
  for (int nelea = 0; nelea <= nele_-2; ++nelea) {
    const int neleb = nele_-2 - nelea;
    if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

    shared_ptr<Kramers<2,ZDvec>> interm = two_down_from_civec(nelea, neleb, ist);
    shared_ptr<Kramers<2,ZDvec>> interm2 = ist == jst ? interm : two_down_from_civec(nelea, neleb, jst);
    // TODO bra-ket symmetry is not used for jst==ist. If this takes too long in CASSCF, please implement transpose below
    for (auto& i : *interm2)
      for (auto& j : *interm) {
        auto rdm2 = make_shared<ZRDM<2>>(norb_);
        auto rdm2grouped = group(group(*rdm2, 2,4), 0,2);
        contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, rdm2grouped, {1,2}, true, false);

        auto sorted = rdm2->clone();
        sort_indices<0,2,1,3,0,1,1,1>(rdm2->data(), sorted->data(), norb_, norb_, norb_, norb_);
        auto ijtag = merge(i.first, j.first);
        out->add(ijtag.perm({{0,2,1,3}}), sorted);
      }
  }
  out->emplace_perm({0,1,2,3},  1.0);
  out->emplace_perm({0,3,2,1}, -1.0);
  out->emplace_perm({2,1,0,3}, -1.0);
  out->emplace_perm({2,3,0,1},  1.0);
  return out;
}


shared_ptr<Kramers<2,ZRDM<1>>> ZHarrison::rdm1(const int jst, const int ist) const {
  auto out = make_shared<Kramers<2,ZRDM<1>>>();
  if (nele_ < 1) return out;

  auto space1 = make_shared<RelSpace>(norb_, nele_-1);
  // loop over n-1 determinant spaces
  for (int nelea = 0; nelea <= nele_-1; ++nelea) {
    const int neleb = nele_-1 - nelea;
    if (nelea > norb_ || neleb > norb_ || neleb < 0) continue;

    shared_ptr<Kramers<1,ZDvec>> interm  = one_down_from_civec(nelea, neleb, ist, space1);
    shared_ptr<Kramers<1,ZDvec>> interm2 = jst == ist ? interm : one_down_from_civec(nelea, neleb, jst, space1);
    // TODO bra-ket symmetry is not used for jst==ist. If this takes too long in CASSCF, please implement transpose below
    for (auto& i : *interm2)
      for (auto& j : *interm) {
        auto rdm1 = make_shared<ZRDM<1>>(norb_);
        contract(1.0, group(*i.second,0,2), {0,1}, group(*j.second,0,2), {0,2}, 0.0, *rdm1, {1,2}, true, false);
        out->add(merge(i.first, j.first), rdm1);
      }
  }
  return out;
}


void ZHarrison::compute_rdm12() {
  // for one-body RDM
  rdm1_.clear();
  rdm2_.clear();
  rdm1_.resize(nstate_);
  rdm2_.resize(nstate_);

  for (int istate = 0; istate != nstate_; ++istate) {
    // one body RDM
    rdm1_[istate] = rdm1(istate, istate);

    // if off-diagonals are zero, generate a blank RDM for completeness
    if (!rdm1_[istate]->exist({1,0}))
      rdm1_[istate]->add({1,0}, rdm1_[istate]->at({0,0})->clone());

    // two body RDM
    rdm2_[istate] = rdm2(istate, istate);

    // append permutation information
    rdm2_[istate]->emplace_perm({{0,3,2,1}},-1);
    rdm2_[istate]->emplace_perm({{2,3,0,1}}, 1);
    rdm2_[istate]->emplace_perm({{2,1,0,3}},-1);
  }

  if (nstate_ > 1) {
    rdm1_av_ = make_shared<Kramers<2,ZRDM<1>>>();
    rdm2_av_ = make_shared<Kramers<4,ZRDM<2>>>();
    for (int istate = 0; istate != nstate_; ++istate) {
      for (auto& i : *rdm1_[istate])
        rdm1_av_->add(i.first, i.second);
      for (auto& i : *rdm2_[istate])
        rdm2_av_->add(i.first, i.second);
    }
    for (auto& i : *rdm1_av_) i.second->scale(1.0/nstate_);
    for (auto& i : *rdm2_av_) i.second->scale(1.0/nstate_);
    rdm2_av_->emplace_perm({{0,3,2,1}},-1);
    rdm2_av_->emplace_perm({{2,3,0,1}}, 1);
    rdm2_av_->emplace_perm({{2,1,0,3}},-1);
  } else {
    rdm1_av_ = rdm1_.front();
    rdm2_av_ = rdm2_.front();
  }

  // set expanded
  rdm1_av_expanded_ = expand_kramers<1,complex<double>>(rdm1_av_, norb_);
  rdm2_av_expanded_ = expand_kramers<2,complex<double>>(rdm2_av_, norb_);
}


vector<shared_ptr<const ZMatrix>> ZHarrison::rdm1_matrix() const {
  // RDM transform as D_rs = C*_ri D_ij (C*_rj)^+
  vector<shared_ptr<const ZMatrix>> out = {};
  for (int i=0; i!=nstate_; ++i) {
    shared_ptr<const ZRDM<1>> tmp = expand_kramers<1,complex<double>>(rdm1_[i], norb_);
    auto mat = make_shared<ZMatrix>(norb_*2, norb_*2);
    copy_n(tmp->data(), tmp->size(), mat->data());
    out.push_back(mat);
  }
  assert(out.size() == nstate_);
  return out;
}


shared_ptr<ZMatrix> ZHarrison::rdm1_av() const {
  auto out = make_shared<ZMatrix>(norb_*2, norb_*2);
  copy_n(rdm1_av_expanded_->data(), out->size(), out->data());
  return out;
}


shared_ptr<ZMatrix> ZHarrison::rdm2_av() const {
  auto out  = make_shared<ZMatrix>(4*norb_*norb_, 4*norb_*norb_);
  copy_n(rdm2_av_expanded_->data(), out->size(), out->data());
  return out;
}
