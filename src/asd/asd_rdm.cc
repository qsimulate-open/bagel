//
// BAGEL - Parallel electron correlation program.
// Filename: asd/asd_rdm.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Inkoo Kim <inkoo.kim@northwestern.edu>
// Maintainer: Shiozaki Group
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

#include <src/asd/asd.h>
#include <src/asd/state_tensor.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;
using namespace btas;

void ASD_base::compute_rdm12_dimer() {

  //TODO: temporary
  const int norbA = dimer_->active_refs().first->nact();
  const int norbB = dimer_->active_refs().second->nact();
  trdm1_ = make_shared<RDM<1>>(norbA+norbB);
  trdm2_ = make_shared<RDM<2>>(norbA+norbB);

  statetensor_ = make_shared<StateTensor>(adiabats_, subspaces_base());
  statetensor_->print();

  for (int i = 0; i != nstates_; ++i) {
    shared_ptr<RDM<1>> rdm1;
    shared_ptr<RDM<2>> rdm2;
    tie(rdm1,rdm2) = compute_rdm12_dimer(i);
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
    tie(r1,r2) = compute_diagonal_block(subspace, istate);
    if (r1) assert(false);
    if (r2) *rdm2 += *r2;
  }

  // off diagonal subspaces
  for (auto iAB = subspaces.begin(); iAB != subspaces.end(); ++iAB) {
    for (auto jAB = subspaces.begin(); jAB != iAB; ++jAB) {
      shared_ptr<RDM<1>> r1;
      shared_ptr<RDM<2>> r2;
      tie(r1,r2) = couple_blocks(*jAB, *iAB, istate); //Lower-triangular (i<->j)
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

  bool flip = (static_cast<int>(term_type) < 0);
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

  auto& A  = keys[0]; auto& B  = keys[1];
  auto& Ap = keys[2]; auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  auto out = make_shared<RDM<2>>(nactA+nactB);

  // alpha-alpha
  auto gamma_AA_alpha = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
  auto gamma_BB_alpha = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});

  // beta-beta
  auto gamma_AA_beta = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
  auto gamma_BB_beta = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});

  auto rdmAA = make_shared<Matrix>(*gamma_AA_alpha % gamma_BB_alpha);
  auto rdmBB = make_shared<Matrix>(*gamma_AA_beta  % gamma_BB_beta);

  auto rdmAB = make_shared<Matrix>(*gamma_AA_alpha % gamma_BB_beta);
  auto rdmBA = make_shared<Matrix>(*gamma_AA_beta  % gamma_BB_alpha);

  {// P(p,q',r',s) : p15
    auto rdmt = rdmAA->clone();
    sort_indices<0,3,2,1, 0,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa
    sort_indices<0,3,2,1, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb
    if (!subdia) {
      sort_indices<3,0,1,2, 1,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa of (N,M)
      sort_indices<3,0,1,2, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
    }
    auto low = {    0, nactA, nactA,     0};
    auto up  = {nactA, nactT, nactT, nactA};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {// d_pqr's' : p19
    auto rdmt = rdmAA->clone();
    sort_indices<0,1,2,3, 0,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa bb
    sort_indices<0,1,2,3, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb aa
    if (!subdia) {
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa of (N,M)
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
      sort_indices<1,0,3,2, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
    }
    auto low = {    0,     0, nactA, nactA};
    auto up  = {nactA, nactA, nactT, nactT};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  if (B == Bp) { //Monomer A
    cout << "A" << endl;
    const int n = B.nstates();
    Matrix delta(n,n);
    delta.unit();
    btas::CRange<2> range(n*n, 1);
    const MatView gamma_B(btas::make_view(range, delta.storage()), /*localized*/true);
    //1RDM
    auto rdmA = make_shared<Matrix>(*gamma_AA_alpha % gamma_B);
    auto rdmB = make_shared<Matrix>(*gamma_AA_beta  % gamma_B);
    auto rdmt = rdmA->clone();
    sort_indices<0,1, 0,1, 1,1>(rdmA->data(), rdmt->data(), nactA, nactA); //a'a
    sort_indices<0,1, 1,1, 1,1>(rdmB->data(), rdmt->data(), nactA, nactA); //b'b
    if (!subdia) {
      sort_indices<1,0, 1,1, 1,1>(rdmA->data(), rdmt->data(), nactA, nactA); //a'a
      sort_indices<1,0, 1,1, 1,1>(rdmB->data(), rdmt->data(), nactA, nactA); //b'b
    }
    auto TEMP = make_shared<RDM<1>>(nactA+nactB);
    auto low = {    0,     0};
    auto up  = {nactA, nactA};
    auto outv = make_rwview(TEMP->range().slice(low, up), TEMP->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    *trdm1_ += *TEMP;
    cout << "A*" << endl;
  }

  if (A == Ap) { //Monomer B
    cout << "B" << endl;
    auto gamma_A = statetensor_->contract_statetensor(keys, istate);
    //1RDM
    auto rdmA = make_shared<Matrix>(*gamma_A % gamma_BB_alpha);
    auto rdmB = make_shared<Matrix>(*gamma_A % gamma_BB_beta);
    auto rdmt = rdmA->clone();
    sort_indices<0,1, 0,1, 1,1>(rdmA->data(), rdmt->data(), nactB, nactB); //a'a
    sort_indices<0,1, 1,1, 1,1>(rdmB->data(), rdmt->data(), nactB, nactB); //b'b
    if (!subdia) {
      sort_indices<1,0, 1,1, 1,1>(rdmA->data(), rdmt->data(), nactB, nactB); //a'a
      sort_indices<1,0, 1,1, 1,1>(rdmB->data(), rdmt->data(), nactB, nactB); //b'b
    }
    auto TEMP = make_shared<RDM<1>>(nactA+nactB);
    auto low = {nactA, nactA};
    auto up  = {nactT, nactT};
    auto outv = make_rwview(TEMP->range().slice(low, up), TEMP->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
    *trdm1_ += *TEMP;
    cout << "B" << endl;
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

  //1RDM
  {
    auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha}, statetensor_, istate);
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    auto rdm = make_shared<Matrix>(*gamma_A % gamma_B);
    auto rdmt = rdm->clone();

    // P(p,q') : p10
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  //2RDM
  {
    auto gamma_A  = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha}, statetensor_, istate);
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A % gamma_B1);
    auto rdm2 = make_shared<Matrix>(*gamma_A % gamma_B2);
    auto rdmt = rdm1->clone();

    // P(p,q',r',s') : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    sort_indices<0,2,1,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }
  {
    auto gamma_A1 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
    auto gamma_A2 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    auto rdm1 = make_shared<Matrix>(*gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(*gamma_A2 % gamma_B);
    auto rdmt = rdm1->clone();

    //P(p,q',r,s) : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
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
  //RDM1
  {
    auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta}, statetensor_, istate);
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    auto rdm = make_shared<Matrix>(*gamma_A % gamma_B);
    auto rdmt = rdm->clone();

    // P(p,q') : p10
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  //RDM2
  {
    auto gamma_A  = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta}, statetensor_, istate);
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A % gamma_B1);
    auto rdm2 = make_shared<Matrix>(*gamma_A % gamma_B2);
    auto rdmt = rdm1->clone();

    // P(p,q',r',s') : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,2,1,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }
  {
    auto gamma_A1 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
    auto gamma_A2 = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta}, statetensor_, istate);
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(*gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(*gamma_A2 % gamma_B);
    auto rdmt = rdm1->clone();

    // P(p,q',r,s) : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA,     0,     0};
    auto up  = {nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(out1,out2);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_abFlip(const array<MonomerKey,4>& keys, const int istate) const {

// if ab-flip, account ba-flip arising from (N,M)
// if(M,N) is ba-flip then (N,M) is ab-flip and this will include ba-flip of (M,N) too.
  auto& B = keys[1];
  auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}));

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  auto out = make_shared<RDM<2>>(nactA+nactB);

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm = make_shared<Matrix>(*gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  // P(p,q',r',s) : p4  ab-flip of (M,N)
  sort_indices<0,3,2,1, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //ab-flip
  sort_indices<1,2,3,0, 1,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //ba-flip of (N,M) p15B

  auto low = {    0, nactA, nactA,     0};
  auto up  = {nactA, nactT, nactT, nactA};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_abET(const array<MonomerKey,4>& keys, const int istate) const {

// for (M,N)
// if inverse ab-ET / compute (N,M)
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p14B
  sort_indices<0,2,1,3, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);
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

//off-diagonal subspaces only!
// if(M,N) is inverse-aa-ET, swap M,N as (N,M) will be aa-ET and contribute to 2RDM
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p1B
  sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


tuple<shared_ptr<RDM<1>>, shared_ptr<RDM<2>>> ASD_base::compute_bbET(const array<MonomerKey,4>& keys, const int istate) const {

// cf. aaET
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::CreateBeta}));

  auto gamma_A = gammatensor_[0]->contract_block_with_statetensor(keys, {GammaSQ::CreateBeta, GammaSQ::CreateBeta}, statetensor_, istate);
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(*gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p14B
  sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr, out);
}


void ASD_base::debug_RDM(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2) const {
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  const int neleA = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int neleB = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int nelec = neleA+neleB;

  cout << "#of active electrons is A : " << neleA << endl;
  cout << "#of active electrons is B : " << neleB << endl;
  cout << "#of total active electrons: " << nelec << endl;

  auto rdm1A = make_shared<RDM<1>>(nactA);
  {
    auto low = {0,0};
    auto up  = {nactA,nactA};
    auto view = btas::make_view(rdm1->range().slice(low,up), rdm1->storage());
    copy(view.begin(), view.end(), rdm1A->begin());
  }
  auto rdm1B = make_shared<RDM<1>>(nactB);
  {
    auto low = {nactA,nactA};
    auto up  = {nactT,nactT};
    auto view = btas::make_view(rdm1->range().slice(low,up), rdm1->storage());
    copy(view.begin(), view.end(), rdm1B->begin());
  }

  auto rdm2A = make_shared<RDM<2>>(nactA);
  {
    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage());
    copy(view.begin(), view.end(), rdm2A->begin());
  }
  auto rdm2B = make_shared<RDM<2>>(nactB);
  {
    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage());
    copy(view.begin(), view.end(), rdm2B->begin());
  }

/*
  auto rdm3A = make_shared<RDM<3>>(nactA);
  {
    auto low = {0,0,0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA,nactA,nactA};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    copy(view.begin(), view.end(), rdm3A->begin());
  }
  auto rdm3B = make_shared<RDM<3>>(nactB);
  {
    auto low = {nactA,nactA,nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT,nactT,nactT};
    auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
    copy(view.begin(), view.end(), rdm3B->begin());
  }
*/

  //1RDM (Trace)
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i) {
      sum += rdm1->element(i,i);
    }
    cout << "1RDM(A)  Trace = " << sum << endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i) {
      sum += rdm1->element(i,i);
    }
    cout << "1RDM(B)  Trace = " << sum << endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i) {
      sum += rdm1->element(i,i);
    }
    cout << "1RDM(AB) Trace = " << sum << endl;
  }

  shared_ptr<Matrix> rdm1_mat = rdm1->rdm1_mat(/*nclose*/0);
  cout << "is 1RDM symmetric? " << rdm1_mat->is_symmetric() << endl;

  //2RDM check (A) Gamma_ij,kl = <0|E_ij,kl|0> = <0|(k1)'(i2)'(j2)(l1)|0>  1,2 = spin
  //diagonal: i=j, k=l
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j) {
      sum += rdm2->element(i,i,j,j);
    }
    cout << "2RDM(A)  Trace = " << sum << endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j) {
      sum += rdm2->element(i,i,j,j);
    }
    cout << "2RDM(B)  Trace = " << sum << endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += rdm2->element(i,i,j,j);
    }
    cout << "2RDM(AB) Trace = " << sum << endl;
  }

  {//RDM2 symmetric check
    for (int l = 0; l != nactT; ++l)
    for (int k = 0; k != nactT; ++k)
    for (int j = 0; j != nactT; ++j)
    for (int i = 0; i != nactT; ++i) {
      double ijkl = rdm2->element(i,j,k,l);
      double klij = rdm2->element(k,l,i,j);
      double jilk = rdm2->element(j,i,l,k);
      double lkji = rdm2->element(l,k,j,i);
      if (abs(ijkl-klij) > 1.0e-10) { //assert(false); //cout << "ERROR1" << endl;
        cout << "ERROR1: " << i << j << k << l << ":" << ijkl << " /= "
                       << k << l << i << j << ":" << klij << endl;
      }
      if (abs(ijkl-jilk) > 1.0e-10) { //assert(false); //cout << "ERROR2" << endl;
        cout << "ERROR2: " << i << j << k << l << ":" << ijkl << " /= "
                       << j << i << l << k << ":" << jilk << endl;
      }
      if (abs(ijkl-lkji) > 1.0e-10) { //assert(false); //cout << "ERROR3" << endl;
        cout << "ERROR3: " << i << j << k << l << ":" << ijkl << " " << lkji << endl;
      }
    }
  }

  { //Gamma_ij,kk
    cout << "2RDM(A) Partial Trace Sum_k (i,j,k,k)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1A);
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      debug->element(i,j) -= 1.0/(neleA-1) * rdm2->element(i,j,k,k);
    }
    debug->print(1.0e-10);
  }
  { //Gamma_ij,kk
    cout << "2RDM(B) Partial Trace Sum_k (i,j,k,k)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1B);
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      debug->element(i-nactA,j-nactA) -= 1.0/(neleB-1) * rdm2->element(i,j,k,k);
    }
    debug->print(1.0e-10);
  }
  { //Gamma_ij,kk
    cout << "2RDM(AB) Partial Trace Sum_k (i,j,k,k)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      debug->element(i,j) -= 1.0/(nelec-1) * rdm2->element(i,j,k,k);
    }
    debug->print(1.0e-10);
  }

/*
  //construct approx twordm
  cout << "Build approx 2RDM:" << endl;
  for (int i = 0; i != nactA; ++i)
  for (int j = nactA; j != nactT; ++j) {
    approx2rdm_->element(i,i,j,j) = rdm1->element(i,i) * rdm1->element(j,j);
    approx2rdm_->element(j,j,i,i) = approx2rdm_->element(i,i,j,j);
    approx2rdm_->element(i,j,j,i) = -0.5*(rdm1->element(i,i) * rdm1->element(j,j));// - rdm1->element(i,i);
    approx2rdm_->element(j,i,i,j) = approx2rdm_->element(i,j,j,i); //-rdm1->element(i,i) * rdm1->element(j,j);// - rdm1->element(j,j);
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += approx2rdm_->element(i,i,j,j);
    }
    cout << "APPROX2RDM(AB) Trace = " << sum << endl;
  }
  { //difference
    auto debug = make_shared<RDM<2>>(*rdm2);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l) {
      debug->element(i,j,k,l) -= approx2rdm_->element(i,j,k,l);
    }
    debug->print(1.0e-8);
    auto low = {0,0,0,0};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(debug->range().slice(low,up), debug->storage());
    auto rdm2 = make_shared<Matrix>(nactT*nactT*nactT,nactT,1);
    copy(view.begin(), view.end(), rdm2->begin());
    cout << "Norm of approx. 2RDM =" << rdm2->norm() << endl;
  }
*/

#if 0
  //3RDM: Gamma_ij,kl,mn
  cout << "-------------- 3RDM --------------" << endl;
  { //Trace A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    cout << "3RDM Trace (A)  = " << sum << endl;
  }
  { //Trace B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    cout << "3RDM Trace (B)  = " << sum << endl;
  }
  { //Trace AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    cout << "3RDM Trace (AB) = " << sum << endl;
  }


  { //Gamma_ij,kl,mm : p21
    cout << "3RDM Partial Trace Sum_m (i,j,k,l,m,m)" << endl;
    auto debug = make_shared<RDM<2>>(*rdm2);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l)
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(i,j,k,l,m,m);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,mm,kl : p21
    cout << "3RDM Partial Trace Sum_m (i,j,m,m,k,l)" << endl;
    auto debug = make_shared<RDM<2>>(*rdm2);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l)
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(i,j,m,m,k,l);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,kl,mm : p21
    cout << "3RDM Partial Trace Sum_m (m,m,i,j,k,l)" << endl;
    auto debug = make_shared<RDM<2>>(*rdm2);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l)
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j,k,l) -= 1.0/(nelec-2) * threerdm_->element(m,m,i,j,k,l);
    }
    debug->print(1.0e-12);
  }
  { //Gamma_ij,kk,mm : p21
    cout << "3RDM Partial Trace Sum_m (i,j,k,k,m,m)" << endl;
    auto debug = make_shared<RDM<1>>(*rdm1);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int m = 0; m != nactT; ++m) {
      debug->element(i,j) -= 1.0/((nelec-2)*(nelec-1)) * threerdm_->element(i,j,k,k,m,m);
    }
    debug->print(1.0e-12);
  }

  cout << "-------------- 4RDM --------------" << endl;

  { //Trace A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k)
    for (int l = 0; l != nactA; ++l) {
      sum += fourrdm_->element(i,i,j,j,k,k,l,l);
    }
    cout << "4RDM Trace (A)  = " << sum << endl;
  }

  { //Trace B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k)
    for (int l = nactA; l != nactT; ++l) {
      sum += fourrdm_->element(i,i,j,j,k,k,l,l);
    }
    cout << "4RDM Trace (B)  = " << sum << endl;
  }

  { //Trace AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k)
    for (int l = 0; l != nactT; ++l) {
      sum += fourrdm_->element(i,i,j,j,k,k,l,l);
    }
    cout << "4RDM Trace (AB)  = " << sum << endl;
  }

  {
    auto debug = make_shared<RDM<3>>(*rdm3A);
    cout << "4RDM(A) debug test" << endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
        debug->element(b,j,c,k,d,l) -= 1.0/(neleA-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
      }
    debug->print(1.0e-8);
  }
  {
    auto debug = make_shared<RDM<3>>(*rdm3B);
    cout << "4RDM(B) debug test" << endl;
    for (int l = 0; l != nactB; ++l)
      for (int d = 0; d != nactB; ++d)
        for (int k = 0; k != nactB; ++k)
          for (int c = 0; c != nactB; ++c)
            for (int j = 0; j != nactB; ++j)
              for (int b = 0; b != nactB; ++b)
      for (int i = 0; i != nactB; ++i) {
        debug->element(b,j,c,k,d,l) -= 1.0/(neleB-3) * fourrdm_->element(i+nactA,i+nactA,b+nactA,j+nactA,c+nactA,k+nactA,d+nactA,l+nactA);
      }
    debug->print(1.0e-8);
  }

  {
    auto debug = make_shared<RDM<3>>(*threerdm_);
    cout << "4RDM debug test 1" << endl;
    for (int l = 0; l != nactT; ++l)
      for (int d = 0; d != nactT; ++d)
        for (int k = 0; k != nactT; ++k)
          for (int c = 0; c != nactT; ++c)
            for (int j = 0; j != nactT; ++j)
              for (int b = 0; b != nactT; ++b)
      for (int i = 0; i != nactT; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelec-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelec-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelec-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelec-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
#endif

}

#if 0


  // Checking 4RDM by comparing with 3RDM
  {
    auto debug = make_shared<RDM<3>>(*threerdm_);
    cout << "4RDM debug test 2" << endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  {
    auto debug = make_shared<RDM<3>>(*threerdm_);
    cout << "4RDM debug test 3" << endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  {
    auto debug = make_shared<RDM<3>>(*threerdm_);
    cout << "4RDM debug test 4" << endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }



  assert(false);

  cout << "Monomer A 4RDM print" << endl;
  for (int i = 0; i != nactA; ++ i)
  for (int j = 0; j != nactA; ++ j)
  for (int k = 0; k != nactA; ++ k)
  for (int l = 0; l != nactA; ++ l)
  for (int m = 0; m != nactA; ++ m)
  for (int n = 0; n != nactA; ++ n)
  for (int o = 0; o != nactA; ++ o)
  for (int p = 0; p != nactA; ++ p) {
    double elem = fourrdm_->element(i,j,k,l,m,n,o,p);
    elem = abs(elem);
    if(elem > 1.0e-8) cout << "RDM4(" << i << j << k << l << m << n << o << p << ") " << fourrdm_->element(i,j,k,l,m,n,o,p) << endl ;
  }
  assert(false);
#endif

void ASD_base::debug_energy(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2) const {
  //Energy calculation
  cout << "!@# Energy calculated from RDM:" << endl;
  const int nclosedA = dimer_->active_refs().first->nclosed();
  const int nclosedB = dimer_->active_refs().second->nclosed();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

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
  int1->print("1e integral",nactT);

  auto rdm1_mat = rdm1->rdm1_mat(0);
  rdm1_mat->print("1RDM",nactT);
  {
    auto trdm1_mat = trdm1_->rdm1_mat(0);
    trdm1_mat->print("1RDM_CHECK",nactT);
  }

  double  e1 = ddot_(nactT*nactT, int1->element_ptr(0,0), 1, rdm1_mat->element_ptr(0,0), 1);
  cout << "1E energy = " << e1 << endl;

  double e2 = 0.0;
  //AAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,0>();
    auto int2 = make_shared<Matrix>(nactA*nactA*nactA*nactA,1);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactA); //conver to chemist not.

    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA,nactA,1); //empty d_AAAA (note: the dimension specification actually do not matter)
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAA filled
    cout << "2E energy (AAAA) = " << 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BBBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,1>();
    auto int2 = make_shared<Matrix>(1,nactB*nactB*nactB*nactB);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactB, nactB, nactB, nactB); //conver to chemist not.

    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBBB sector of d
    auto rdm2 = make_shared<Matrix>(1,nactB*nactB*nactB*nactB); //empty d_BBBB
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBBB filled
    cout << "2E energy (BBBB) = " << 0.5 * ddot_(nactB*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactB*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  cout << "(3A,1B) part:" << endl;
  //AAAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,1>(); // <pq|rs'> in (pqr,s') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactA*nactB,1);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.
    auto low = {    0,    0,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AAAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA*nactB,1); //empty d_AAAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAB filled
    cout << "2E energy (AAAB) = " << 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //AABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,0>(); // <pq'|rs> in (prs,q') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1);
    sort_indices<0,1,3,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [ij|k'l]

    auto low = {    0,    0,nactA,    0};
    auto up  = {nactA,nactA,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1); //empty d_AABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABA filled
    cout << "2E energy (AABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,0>(); // <p'q|rs> in (qrs,p') format
    auto int2 = make_shared<Matrix>(nactB,nactA*nactA*nactA);
    sort_indices<3,1,0,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.

    auto low = {nactA,    0,    0,    0};
    auto up  = {nactT,nactA,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactB,nactA*nactA*nactA); //empty d_BAAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BAAA filled
    cout << "2E energy (BAAA) = " << 0.5 * ddot_(nactB*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,0>(); // <pq|r's> in (pqs,r') format
    auto int2 = make_shared<Matrix>(nactA*nactB*nactA*nactA,1);
    sort_indices<0,3,1,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [pr'|qs]=[ij'|kl]

    auto low = {    0,nactA,    0,    0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactA*nactA,1); //empty d_ABAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABAA filled
    cout << "2E energy (ABAA) = " << 0.5 * ddot_(nactA*nactB*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }



  cout << "(1A,3B) part:" << endl;
  //ABBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,1>(); // <pq'|r's'> in (p,q'r's') format
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [pr'|q's']=[ij'|k'l']

    auto low = {    0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_ABBB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (ABBB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,1>(); // <p'q'|rs'> in (r,p'q's') format
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    sort_indices<1,0,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r|q's']=[i'j|k'l']

    auto low = {nactA,    0,nactA,nactA};
    auto up  = {nactT,nactA,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BABB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,1>(); // <p'q|r's'> in (q,p'r's') format
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    sort_indices<1,2,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|qs']=[i'j'|kl']

    auto low = {nactA,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BBAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BBAB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,0>(); // <p'q'|r's> in (s,p'q'r') format
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    sort_indices<1,3,2,0, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|q's]=[i'j'|k'l]

    auto low = {nactA,nactA,nactA,    0};
    auto up  = {nactT,nactT,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBBA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BBBA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBBA filled
    cout << "2E energy (BBBA) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  cout << "(2A,2B) part:" << endl;
  //AABB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,1>(); // <pq'|rs'> in (pr,q's') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<0,1,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr|q's']=[ij|k'l']

    auto low = {    0,    0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_AABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABB filled
    cout << "2E energy (AABB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,0>(); // <p'q|r's> in (qs,p'r') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<2,3,0,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r'|qs]=[i'j'|kl]

    auto low = {nactA,nactA,    0,    0};
    auto up  = {nactT,nactT,nactA,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BBAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BBAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBAA filled
    cout << "2E energy (BBAA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,1>(); // <pq|r's'> in (pq,r's') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|qs']=[ij'|kl']

    auto low = {    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABAB filled
    cout << "2E energy (ABAB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,0>(); // <p'q'|rs> in (rs,p'q') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<2,0,3,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|q's]=[i'j|k'l]

    auto low = {nactA,    0,nactA,    0};
    auto up  = {nactT,nactA,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BABA filled
    cout << "2E energy (BABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,0>(); // <pq'|r's> in (ps,q'r') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<0,3,2,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|q's]=[ij'|k'l]

    auto low = {    0,nactA,nactA,    0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_ABBA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABBA
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBA filled
    cout << "2E energy (ABBA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BAAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,1>(); // <p'q|rs'> in (qr,p's') format
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    sort_indices<2,1,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|qs']=[i'j|kl']

    auto low = {nactA,    0,    0,nactA};
    auto up  = {nactT,nactA,nactA,nactT};
    auto view = btas::make_view(rdm2->range().slice(low,up), rdm2->storage()); //d_BAAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BAAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_BAAB filled
    cout << "2E energy (BAAB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //Energy print
  cout << "nuclear repulsion= " << dimer_->sref()->geom()->nuclear_repulsion() << endl;
  cout << "core energy      = " << jop_->core_energy() << endl;
  cout << "nuc + core energy= " << dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() << endl;
  cout << "1E energy = " <<  e1 << endl;
  cout << "2E energy = " <<  e2 << endl;
  cout << "Total energy = " << dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy() + e1 + e2 << endl;
}


void ASD_base::symmetrize_rdm12(shared_ptr<RDM<1>>& rdm1, shared_ptr<RDM<2>>& rdm2) const {
  cout << "!@# Unsymmetrized 1RDM" << endl;
  rdm1->print(1.0e-6);

  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA + nactB;

  //Symmetrize: D_AB (calculated) D_BA (uncalc.& symmetrized here)
  auto matBA = make_shared<Matrix>(nactB,nactA); //D_BA empty
  {
    auto low = {0, nactA};
    auto up  = {nactA, nactT};
    auto view = btas::make_view(rdm1->range().slice(low,up), rdm1->storage()); //D_AB sector of D (read ptr)
    auto matAB = make_shared<Matrix>(nactA,nactB); //D_AB empty
    copy(view.begin(), view.end(), matAB->begin()); //D_AB filled
    sort_indices<1,0, 0,1, 1,1>(matAB->data(), matBA->data(), nactA, nactB); // transpose and fill D_BA
  }
  {
    auto low = {nactA, 0};
    auto up  = {nactT, nactA};
    auto outv = btas::make_rwview(rdm1->range().slice(low,up), rdm1->storage()); //D_BA sector of D (read & write ptr)
    copy(matBA->begin(), matBA->end(), outv.begin()); //copy D_BA -> D_BA sector of D
  }

  cout << "!@# Symmetrized 1RDM" << endl;
  rdm1->print(1.0e-6);

  //Symmetrize: d(ABAA) note p18B, 19B
  {
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

  //Symmetrize: d(ABBB) note p18B, 19B
  {
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


  //Symmetrize: d(ABBA) note p19
  {
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

  //Symmetrize: d(AABB) note 19
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

  //Symmetrize: d(ABAB) note p19
  {
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

}
