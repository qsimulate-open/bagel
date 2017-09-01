//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd/asd_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/asd.h>
#include <src/asd/state_tensor.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace btas;

void ASD_base::compute_rdm12_dimer() {

  statetensor_ = make_shared<StateTensor>(adiabats_, subspaces_base());
  if (print_info_)
    statetensor_->print();

  for (int i = 0; i != nstates_; ++i) {
    shared_ptr<RDM<1>> rdm1;
    shared_ptr<RDM<2>> rdm2;
    tie(rdm1, rdm2) = compute_rdm12_dimer(i);
    rdm1_[i] = rdm1;
    rdm2_[i] = rdm2;
  }

}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_rdm12_dimer(const int istate) const {
  const int norbA = dimer_->active_refs().first->nact();
  const int norbB = dimer_->active_refs().second->nact();

  auto rdm1 = make_shared<RDM<1>>(norbA+norbB);
  auto rdm2 = make_shared<RDM<2>>(norbA+norbB);

  const auto subspaces = subspaces_base();

  // diagonal subspaces
  for (auto& subspace : subspaces) {
    shared_ptr<RDM<1>> r1;
    shared_ptr<RDM<2>> r2;
    tie(r1, r2) = compute_diagonal_block(subspace, istate);
    if (r1) *rdm1 += *r1;
    if (r2) *rdm2 += *r2;
  }

  // off diagonal subspaces
  for (auto iAB = subspaces.begin(); iAB != subspaces.end(); ++iAB) {
    for (auto jAB = subspaces.begin(); jAB != iAB; ++jAB) {
      shared_ptr<RDM<1>> r1;
      shared_ptr<RDM<2>> r2;
      tie(r1, r2) = couple_blocks(*jAB, *iAB, istate); //Lower-triangular (i<->j)
      if (r1) *rdm1 += *r1;
      if (r2) *rdm2 += *r2;
    }
  }

  // fill up redundnat part
  symmetrize_rdm12(rdm1, rdm2);

  return make_tuple(rdm1, rdm2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_diagonal_block(const DimerSubspace_base& subspace, const int istate ) const {

  array<MonomerKey,4> keys {{ subspace.monomerkey<0>(), subspace.monomerkey<1>(), subspace.monomerkey<0>(), subspace.monomerkey<1>() }};
  auto out = compute_inter_2e(keys, istate, /*subspace diagonal*/true);

  return out;
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::couple_blocks(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp, const int istate) const {

  Coupling term_type = coupling_type(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  const bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    swap(space1,space2);
  }

  tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> out;
  array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = make_tuple(nullptr,nullptr); break;
    case Coupling::diagonal :
      out = compute_inter_2e(keys, istate, /*subspace diagonal*/false); break;
    case Coupling::aET :
      out = compute_aET(keys, istate); break;
    case Coupling::bET :
      out = compute_bET(keys, istate); break;
    case Coupling::abFlip :
      out = compute_abFlip(keys, istate); break;
    case Coupling::abET :
      out = compute_abET(keys, istate); break;
    case Coupling::aaET :
      out = compute_aaET(keys, istate); break;
    case Coupling::bbET :
      out = compute_bbET(keys, istate); break;
    default :
      throw logic_error("Asking for a coupling type that has not been written.");
  }

  return out;
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_inter_2e(const array<MonomerKey,4>& keys, const int istate, const bool subdia) const {

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  auto out  = make_shared<RDM<2>>(nactA+nactB);

  // alpha-alpha
  auto gamma_AA_alpha = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
  auto gamma_BB_alpha = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});

  // beta-beta
  auto gamma_AA_beta = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
  auto gamma_BB_beta = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});

  auto rdmAA = make_shared<Matrix>(*gamma_AA_alpha % gamma_BB_alpha); //a'a|a'a
  auto rdmBB = make_shared<Matrix>(*gamma_AA_beta  % gamma_BB_beta);  //b'b|b'b

  auto rdmAB = make_shared<Matrix>(*gamma_AA_alpha % gamma_BB_beta);  //a'a|b'b
  auto rdmBA = make_shared<Matrix>(*gamma_AA_beta  % gamma_BB_alpha); //b'b|a'a

  {//d_aBBa
    auto rdmt = rdmAA->clone();
    sort_indices<0,3,2,1, 0,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //p'qR'S = -p'R'qS => -d_pSRq (0321)
    sort_indices<0,3,2,1, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    if (!subdia) {
      sort_indices<1,2,3,0, 1,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //q'pS'R = -q'S'pR => -d_qRSp (1230)
      sort_indices<1,2,3,0, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    }
    auto low = {    0, nactA, nactA,     0};
    auto up  = {nactA, nactT, nactT, nactA};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {//d_aaBB
    auto rdmt = rdmAA->clone();
    sort_indices<0,1,2,3, 0,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //p'qR'S = p'R'Sq => d_pqRS (0123)
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    if (!subdia) {
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //q'S'Rp => d_qpSR (1032)
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB);
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB);
    }
    auto low = {    0,     0, nactA, nactA};
    auto up  = {nactA, nactA, nactT, nactT};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_aET(const array<MonomerKey,4>& keys, const int istate) const {

  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out1 = make_shared<RDM<1>>(nactA+nactB);
  auto out2 = make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  {//d_aB
    auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha}, statetensor_, istate);
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    auto rdm = make_shared<Matrix>(*gamma_A % gamma_B); //a'|a
    auto rdmt = rdm->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {//d_aBBB
    auto gamma_A  = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha}, statetensor_, istate);
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A % gamma_B1); //a'|a'aa
    auto rdm2 = make_shared<Matrix>(*gamma_A % gamma_B2); //a'|b'ab
    auto rdmt = rdm1->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB); //p'Q'RS => d_pSQR (0312)
    sort_indices<0,2,1,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB); //p'Q'SR => -d_pRQS (0213)
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {//d_aBaa
    auto gamma_A1 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
    auto gamma_A2 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    auto rdm1 = make_shared<Matrix>(*gamma_A1 % gamma_B); //a'a'a|a
    auto rdm2 = make_shared<Matrix>(*gamma_A2 % gamma_B); //a'b'b|a
    auto rdmt = rdm1->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,3,1,2, 0,1, 1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB); //p'q'rS => d_pSqr (0312)
    sort_indices<0,3,1,2, 1,1, 1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA,     0,     0};
    auto up  = {nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(out1, out2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_bET(const array<MonomerKey,4>& keys, const int istate) const {

  auto& Ap = keys[2];

  auto& B  = keys[1];
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out1 = make_shared<RDM<1>>(nactA+nactB);
  auto out2 = make_shared<RDM<2>>(nactA+nactB);

  const int neleA = Ap.nelea() + Ap.neleb();

  {//d_aB
    auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta}, statetensor_, istate);
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    auto rdm = make_shared<Matrix>(*gamma_A % gamma_B); //b'|b
    auto rdmt = rdm->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {//d_aBBB
    auto gamma_A  = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta}, statetensor_, istate);
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A % gamma_B1); //b'|a'ba
    auto rdm2 = make_shared<Matrix>(*gamma_A % gamma_B2); //b'|b'bb
    auto rdmt = rdm1->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,2,1,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB); //see aET
    sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {//d_aBaa
    auto gamma_A1 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
    auto gamma_A2 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A1 % gamma_B); //b'a'a|b
    auto rdm2 = make_shared<Matrix>(*gamma_A2 % gamma_B); //b'b'b|b
    auto rdmt = rdm1->clone();

    const int fac = neleA%2 == 0 ? 1 : -1;
    sort_indices<0,3,1,2, 0,1, 1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB); //see aET
    sort_indices<0,3,1,2, 1,1, 1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA,     0,     0};
    auto up  = {nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(out1,out2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_abFlip(const array<MonomerKey,4>& keys, const int istate) const {

  auto& B = keys[1];
  auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}));

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  auto out = make_shared<RDM<2>>(nactA+nactB);

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm = make_shared<Matrix>(*gamma_A % gamma_B); //b'a|a'b
  auto rdmt = rdm->clone();

  //d_aBBa
  sort_indices<0,3,2,1, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //see inter_2e
  sort_indices<1,2,3,0, 1,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //lower-triangular

  auto low = {    0, nactA, nactA,     0};
  auto up  = {nactA, nactT, nactT, nactA};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_abET(const array<MonomerKey,4>& keys, const int istate) const {

  auto& B = keys[1];
  auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B); //a'b'|ab
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  sort_indices<0,2,1,3, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //p'q'RS = -p'q'SR => d_pRqS (0213)
  sort_indices<1,3,0,2, 1,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_aaET(const array<MonomerKey,4>& keys, const int istate) const {

  auto& B = keys[1];
  auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B); //a'a'|aa
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //p'q'RS => d_pSqR (0312)

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_bbET(const array<MonomerKey,4>& keys, const int istate) const {

  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::CreateBeta}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateBeta}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B); //b'b'|bb
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //see aaET

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


void ASD_base::print_rdm_info(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2, const int istate) const {
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  const int neleA = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int neleB = 2*(dimer_->isolated_refs().second->nclosed() - dimer_->active_refs().second->nclosed());
  const int nelec = neleA + neleB - charge_;

  cout << "=== RDM information: state(" << istate << ") ===" << endl;

  cout << "Nelectron A : " << neleA << endl;
  cout << "Nelectron B : " << neleB << endl;
  cout << "Charge      : " << charge_ << endl;
  cout << "Total elec  : " << nelec << endl;

  {//1RDM
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
      sum += rdm1->element(i,i);
    cout << "1RDM Trace = " << setw(20) << setprecision(8) <<sum << endl;
    assert(fabs((sum - static_cast<double>(nelec)) / static_cast<double>(nelec)) < 1.0e-8);
  }

  {//2RDM Trace: Gamma_ij,kl = <0|E_ij,kl|0> = <0|(k1)'(i2)'(j2)(l1)|0>  1,2 = spin
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
      for (int j = 0; j != nactT; ++j)
      sum += rdm2->element(i,i,j,j);
    cout << "2RDM Trace = " << setw(20) << setprecision(8) <<sum << endl;
    assert(fabs((sum - static_cast<double>(nelec*(nelec-1))) / static_cast<double>(nelec*(nelec-1))) < 1.0e-8);
  }

  {//Partial trace: Gamma_ij,kk
    auto debug = make_shared<RDM<1>>(*rdm1);
    for (int i = 0; i != nactT; ++i)
      for (int j = 0; j != nactT; ++j)
        for (int k = 0; k != nactT; ++k)
          debug->element(i,j) -= 1.0/(nelec-1) * rdm2->element(i,j,k,k);

    for (int i = 0; i !=nactT; ++i)
      for (int j = 0; j !=nactT; ++j)
        if(fabs(debug->element(j,i)) > 1.0e-8) {
          cout << j << " " << i << " : " << debug->element(j,i) << endl;
          throw runtime_error("Partial trace check failed");
        }
  }
}


void ASD_base::print_energy_info(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2, const int istate) const {
  const int nclosedA = dimer_->active_refs().first->nclosed();
  const int nclosedB = dimer_->active_refs().second->nclosed();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  cout << "=== Energy calculated using RDM: state(" << istate << ") ===" << endl;
  cout << "Number of closed orbitals: A(" << nclosedA << "), B(" << nclosedB << ")" << endl;
  cout << "Number of active orbitals: A(" << nactA << "), B(" << nactB << ")" << endl;

  shared_ptr<const Matrix> ha = jop_->monomer_jop<0>()->mo1e()->matrix(); //h_AA
  shared_ptr<const Matrix> hb = jop_->monomer_jop<1>()->mo1e()->matrix(); //h_BB
  //                                                    CSymMatrix -> Matrix conversion
  shared_ptr<const Matrix> hc = jop_->cross_mo1e(); //h_AB

  auto int1 = make_shared<Matrix>(nactT,nactT);
  int1->zero();
  int1->copy_block(0,0,ha->ndim(),ha->mdim(),ha);
  int1->copy_block(nactA,nactA,hb->ndim(),hb->mdim(),hb);
  int1->copy_block(0,nactA,hc->ndim(),hc->mdim(),hc);
  int1->copy_block(nactA,0,hc->mdim(),hc->ndim(),hc->transpose());

  auto rdm1_mat = rdm1->rdm1_mat(/*nclosed*/0);

  const double e1 = blas::dot_product(int1->element_ptr(0,0), nactT*nactT, rdm1_mat->element_ptr(0,0));

  double e1_aa;
  {//AA
    shared_ptr<const Matrix> i = int1->get_submatrix(0, 0, nactA, nactA);
    shared_ptr<const Matrix> r = rdm1_mat->get_submatrix(0, 0, nactA, nactA);
    e1_aa = blas::dot_product(i->element_ptr(0,0), nactA*nactA, r->element_ptr(0,0));
  }

  double e1_bb;
  {//BB
    shared_ptr<const Matrix> i = int1->get_submatrix(nactA, nactA, nactB, nactB);
    shared_ptr<const Matrix> r = rdm1_mat->get_submatrix(nactA, nactA, nactB, nactB);
    e1_bb = blas::dot_product(i->element_ptr(0,0), nactB*nactB, r->element_ptr(0,0));
  }

  double e1_ab;
  {//AB
    shared_ptr<const Matrix> i = int1->get_submatrix(0, nactA, nactA, nactB);
    shared_ptr<const Matrix> r = rdm1_mat->get_submatrix(0, nactA, nactA, nactB);
    e1_ab = blas::dot_product(i->element_ptr(0,0), nactA*nactB, r->element_ptr(0,0));
  }

  double e1_ba;
  {//BA
    shared_ptr<const Matrix> i = int1->get_submatrix(nactA, 0, nactB, nactA);
    shared_ptr<const Matrix> r = rdm1_mat->get_submatrix(nactA, 0, nactB, nactA);
    e1_ba = blas::dot_product(i->element_ptr(0,0), nactB*nactA, r->element_ptr(0,0));
  }

  double e2_aaaa;
  {//AAAA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,0>();
    auto int2 = make_shared<VectorB>(nactA*nactA*nactA*nactA);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactA); //conver to chemist not.

    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AAAA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactA*nactA); //empty d_AAAA (note: the dimension specification actually do not matter)
    copy(view.begin(), view.end(), rdm2->data()); //d_AAAA filled
    e2_aaaa = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactA*nactA, rdm2->data());
  }

  double e2_bbbb;
  //BBBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,1>();
    auto int2 = make_shared<VectorB>(nactB*nactB*nactB*nactB);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactB, nactB, nactB, nactB); //conver to chemist not.

    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBBB sector of d
    auto rdm2 = make_shared<VectorB>(nactB*nactB*nactB*nactB); //empty d_BBBB
    copy(view.begin(), view.end(), rdm2->data()); //d_BBBB filled
    e2_bbbb = 0.5 * blas::dot_product(int2->data(), nactB*nactB*nactB*nactB, rdm2->data());
  }

  double e2_aaab;
  {//AAAB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,1>(); // <pq|rs'> in (pqr,s') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactA*nactB);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.

    auto low = {    0,    0,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AAAB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactA*nactB); //empty d_AAAB
    copy(view.begin(), view.end(), rdm2->data()); //d_AAAB filled
    e2_aaab = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactA*nactB, rdm2->data());
  }

  double e2_aaba;
  {//AABA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,0>(); // <pq'|rs> in (prs,q') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactA);
    sort_indices<0,1,3,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [ij|k'l]

    auto low = {    0,    0,nactA,    0};
    auto up  = {nactA,nactA,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AABA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactA); //empty d_AABA
    copy(view.begin(), view.end(), rdm2->data()); //d_AABA filled
    e2_aaba = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactA*nactB, rdm2->data());
  }

  double e2_abaa;
  {//ABAA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,0>(); // <pq|r's> in (pqs,r') format
    auto int2 = make_shared<VectorB>(nactA*nactB*nactA*nactA);
    sort_indices<0,3,1,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [pr'|qs]=[ij'|kl]

    auto low = {    0,nactA,    0,    0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactB*nactA*nactA); //empty d_ABAA
    copy(view.begin(), view.end(), rdm2->data()); //d_ABAA filled
    e2_abaa = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactA*nactB, rdm2->data());
  }

  double e2_baaa;
  {//BAAA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,0>(); // <p'q|rs> in (qrs,p') format
    auto int2 = make_shared<VectorB>(nactB*nactA*nactA*nactA);
    sort_indices<3,1,0,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.

    auto low = {nactA,    0,    0,    0};
    auto up  = {nactT,nactA,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAA sector of d
    auto rdm2 = make_shared<VectorB>(nactB*nactA*nactA*nactA); //empty d_BAAA
    copy(view.begin(), view.end(), rdm2->data()); //d_BAAA filled
    e2_baaa = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactA*nactB, rdm2->data());
  }

  double e2_abbb;
  {//ABBB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,1>(); // <pq'|r's'> in (p,q'r's') format
    auto int2 = make_shared<VectorB>(nactA*nactB*nactB*nactB);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [pr'|q's']=[ij'|k'l']

    auto low = {    0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactB*nactB*nactB); //empty d_ABBB
    copy(view.begin(), view.end(), rdm2->data()); //d_ABBB filled
    e2_abbb = 0.5 * blas::dot_product(int2->data(), nactA*nactB*nactB*nactB, rdm2->data());
  }

  double e2_babb;
  {//BABB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,1>(); // <p'q'|rs'> in (r,p'q's') format
    auto int2 = make_shared<VectorB>(nactA*nactB*nactB*nactB);
    sort_indices<1,0,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r|q's']=[i'j|k'l']

    auto low = {nactA,    0,nactA,nactA};
    auto up  = {nactT,nactA,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BABB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactB*nactB*nactB); //empty d_BABB
    copy(view.begin(), view.end(), rdm2->data()); //d_ABBB filled
    e2_babb = 0.5 * blas::dot_product(int2->data(), nactA*nactB*nactB*nactB, rdm2->data());
  }

  double e2_bbab;
  {//BBAB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,1>(); // <p'q|r's'> in (q,p'r's') format
    auto int2 = make_shared<VectorB>(nactA*nactB*nactB*nactB);
    sort_indices<1,2,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|qs']=[i'j'|kl']

    auto low = {nactA,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactB*nactB*nactB); //empty d_BBAB
    copy(view.begin(), view.end(), rdm2->data()); //d_ABBB filled
    e2_bbab = 0.5 * blas::dot_product(int2->data(), nactA*nactB*nactB*nactB, rdm2->data());
  }

  double e2_bbba;
  {//BBBA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,0>(); // <p'q'|r's> in (s,p'q'r') format
    auto int2 = make_shared<VectorB>(nactA*nactB*nactB*nactB);
    sort_indices<1,3,2,0, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|q's]=[i'j'|k'l]

    auto low = {nactA,nactA,nactA,    0};
    auto up  = {nactT,nactT,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBBA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactB*nactB*nactB); //empty d_BBBA
    copy(view.begin(), view.end(), rdm2->data()); //d_BBBA filled
    e2_bbba = 0.5 * blas::dot_product(int2->data(), nactA*nactB*nactB*nactB, rdm2->data());
  }

  double e2_aabb;
  {//AABB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,1>(); // <pq'|rs'> in (pr,q's') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<0,1,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr|q's']=[ij|k'l']

    auto low = {    0,    0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AABB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_AABB
    copy(view.begin(), view.end(), rdm2->data()); //d_AABB filled
    e2_aabb = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  double e2_bbaa;
  {//BBAA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,0>(); // <p'q|r's> in (qs,p'r') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<2,3,0,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r'|qs]=[i'j'|kl]

    auto low = {nactA,nactA,    0,    0};
    auto up  = {nactT,nactT,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_BBAA
    copy(view.begin(), view.end(), rdm2->data()); //d_BBAA filled
    e2_bbaa = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  double e2_abab;
  {//ABAB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,1>(); // <pq|r's'> in (pq,r's') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|qs']=[ij'|kl']

    auto low = {    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_ABAB
    copy(view.begin(), view.end(), rdm2->data()); //d_ABAB filled
    e2_abab = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  double e2_baba;
  {//BABA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,0>(); // <p'q'|rs> in (rs,p'q') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<2,0,3,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|q's]=[i'j|k'l]

    auto low = {nactA,    0,nactA,    0};
    auto up  = {nactT,nactA,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BABA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_BABA
    copy(view.begin(), view.end(), rdm2->data()); //d_BABA filled
    e2_baba = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  double e2_abba;
  {//ABBA
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,0>(); // <pq'|r's> in (ps,q'r') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<0,3,2,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|q's]=[ij'|k'l]

    auto low = {    0,nactA,nactA,    0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBA sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_ABBA
    copy(view.begin(), view.end(), rdm2->data()); //d_ABBA filled
    e2_abba = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  double e2_baab;
  {//BAAB
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,1>(); // <p'q|rs'> in (qr,p's') format
    auto int2 = make_shared<VectorB>(nactA*nactA*nactB*nactB);
    sort_indices<2,1,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|qs']=[i'j|kl']

    auto low = {nactA,    0,    0,nactA};
    auto up  = {nactT,nactA,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAB sector of d
    auto rdm2 = make_shared<VectorB>(nactA*nactA*nactB*nactB); //empty d_BAAB
    copy(view.begin(), view.end(), rdm2->data()); //d_BAAB filled
    e2_baab = 0.5 * blas::dot_product(int2->data(), nactA*nactA*nactB*nactB, rdm2->data());
  }

  const double e2 = e2_aaaa + e2_bbbb + e2_aabb + e2_bbaa +
                    e2_aaab + e2_aaba + e2_abaa + e2_baaa +
                    e2_abbb + e2_babb + e2_bbab + e2_bbba +
                    e2_abab + e2_baba + e2_abba + e2_baab;
  const double etot = dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() + e1 + e2;

  cout << "Nuc + Core energy   = " << dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() << endl;
  cout << "  Nuclear repulsion = " << dimer_->sref()->geom()->nuclear_repulsion() << endl;
  cout << "  Core energy       = " << jop_->core_energy() << endl;
  cout << "One-electron energy = " << e1 << endl;
  cout << "  2A       AA       = " << e1_aa << endl;
  cout << "  2B       BB       = " << e1_bb << endl;
  cout << "  1A1B     AB       = " << e1_ab << endl;
  cout << "           BA       = " << e1_ba << endl;
  cout << "Two-electron energy = " << e2 << endl;
  cout << "  4A       AAAA     = " << e2_aaaa << endl;
  cout << "  4B       BBBB     = " << e2_bbbb << endl;
  cout << "  3A/1B    AAAB     = " << e2_aaab << endl;
  cout << "           AABA     = " << e2_aaba << endl;
  cout << "           ABAA     = " << e2_abaa << endl;
  cout << "           BAAA     = " << e2_baaa << endl;
  cout << "  1A/3B    ABBB     = " << e2_abbb << endl;
  cout << "           BABB     = " << e2_babb << endl;
  cout << "           BBAB     = " << e2_bbab << endl;
  cout << "           BBBA     = " << e2_bbba << endl;
  cout << "  2A/2B    AABB     = " << e2_aabb << endl;
  cout << "           BBAA     = " << e2_bbaa << endl;
  cout << "           ABAB     = " << e2_abab << endl;
  cout << "           BABA     = " << e2_baba << endl;
  cout << "           ABBA     = " << e2_abba << endl;
  cout << "           BAAB     = " << e2_baab << endl;
  cout << "Total energy = " << etot << endl;

  assert(fabs((etot - energy(istate)) / etot) < 1.0e-10);

}


void ASD_base::symmetrize_rdm12(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2) const {
  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA + nactB;

  {//D(AB)
    auto low = {0, nactA};
    auto up  = {nactA, nactT};
    auto view = btas::make_view(rdm1->range().slice(low,up), rdm1->storage()); //D_AB sector of D (read ptr)
    auto matAB = make_shared<Matrix>(nactA,nactB); //D_AB empty
    copy(view.begin(), view.end(), matAB->begin()); //D_AB filled
    {//D(BA)
      auto matBA = make_shared<Matrix>(nactB,nactA); //D_BA empty
      sort_indices<1,0, 0,1, 1,1>(matAB->data(), matBA->data(), nactA, nactB); // transpose and fill D_BA
      auto low = {nactA, 0};
      auto up  = {nactT, nactA};
      auto outv = btas::make_rwview(rdm1->range().slice(low,up), rdm1->storage()); //D_BA sector of D (read & write ptr)
      copy(matBA->begin(), matBA->end(), outv.begin()); //copy D_BA -> D_BA sector of D
    }
  }

  //1RDM symmetry
  assert(rdm1->rdm1_mat(/*nclosed*/0)->is_symmetric());

  {//d(ABAA)
    auto low = {0,nactA,0,0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactA); //empty d_ABAA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAA filled
    { //d(AAAB)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactA*nactB); //empty d_AAAB
      sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AAAB
      auto low = {0,0,0,nactA};
      auto up  = {nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_AAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AAAB into d_AAAB sector of d
    }
    { //d(BAAA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactA); //empty d_BAAA
      sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_BAAA
      auto low = {nactA,0,0,0};
      auto up  = {nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAA into d_BAAA sector of d
    }
    { //d(AABA)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactB*nactA); //empty d_AABA
      sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AABA
      auto low = {0,0,nactA,0};
      auto up  = {nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_AABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AABA into d_AABA sector of d
    }
  }

  {//d(ABBB)
    auto low = {0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactB); //empty d_ABBB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBB filled
    { //d(BBAB)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactA*nactB); //empty d_BBAB
      sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBAB
      auto low = {nactA,nactA,0,nactA};
      auto up  = {nactT,nactT,nactA,nactT};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAB into d_BBAB sector of d
    }
    { //d(BABB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactB); //empty d_BABB
      sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BABB
      auto low = {nactA,0,nactA,nactA};
      auto up  = {nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BABB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BABB into d_BABB sector of d
    }
    { //d(BBBA)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactB*nactA); //empty d_BBBA
      sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBBA
      auto low = {nactA,nactA,nactA,0};
      auto up  = {nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BBBA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BBBA sector of d
    }
  }

  {//d(ABBA)
    auto low = {    0,nactA,nactA,    0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactA); //empty d_ABBA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBA filled
    { //d(BAAB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactB); //empty d_BAAB
      sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_BAAB
      auto low = {nactA,    0,    0,nactA};
      auto up  = {nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAB into d_BAAB sector of d
    }
  }

  { //d(AABB)
    auto low = {0,0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AABB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), inmat->begin()); //d_AABB filled
    { //d(BBAA)
      auto outmat = make_shared<Matrix>(nactB*nactB*nactA*nactA,1); //empty d_BBAA
      sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactB, nactB); //reorder and fill d_BBAA
      auto low = {nactA,nactA,0,0};
      auto up  = {nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAA into d_BBAA sector of d
    }
  }

  {//d(ABAB)
    auto low = {    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactB); //empty d_ABAB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAB filled
    { //d(BABA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactA); //empty d_BABA
      sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactB); //reorder and fill d_BABA
      auto low = {nactA,    0,nactA,    0};
      auto up  = {nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(rdm2->range().slice(low,up), rdm2->storage()); //d_BABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BABA sector of d
    }
  }

  {//2RDM symmetry
    for (int l = 0; l != nactT; ++l)
      for (int k = 0; k != nactT; ++k)
        for (int j = 0; j != nactT; ++j)
          for (int i = 0; i != nactT; ++i) {
            assert(fabs(rdm2->element(i,j,k,l) - rdm2->element(k,l,i,j)) < 1.0e-10);
            assert(fabs(rdm2->element(i,j,k,l) - rdm2->element(j,i,l,k)) < 1.0e-10);
            assert(fabs(rdm2->element(i,j,k,l) - rdm2->element(l,k,j,i)) < 1.0e-10);
      }
  }

}
