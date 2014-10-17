//
// BAGEL - Parallel electron correlation program.
// Filename: asd_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_base.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

ASD_base::ASD_base(const shared_ptr<const PTree> input, shared_ptr<const Dimer> dimer) : dimer_(dimer)
{
  nstates_ = input->get<int>("nstates", 10);
  max_iter_ = input->get<int>("max_iter", 50);
  davidson_subspace_ = input->get<int>("davidson_subspace", 10);
  nguess_ = input->get<int>("nguess", 10*nstates_);
  dipoles_ = input->get<bool>("dipoles", false);
  thresh_ = input->get<double>("thresh", 1.0e-7);
  print_thresh_ = input->get<double>("print_thresh", 0.01);
  store_matrix_ = input->get<bool>("store_matrix", false);
  charge_ = input->get<int>("charge", 0);
  nspin_ = input->get<int>("spin", 0);

  shared_ptr<const PTree> model_input = input->get_child_optional("models");
  if (model_input) {
    shared_ptr<const PTree> pruning_input = model_input->get_child_optional("pruned");
    vector<int> pruned = (pruning_input ? model_input->get_vector<int>("pruned")
                                             : vector<int>());
    pruned.push_back(-1);

    for (auto& p : pruned) {
      vector<ModelBlock> this_model;
      shared_ptr<const PTree> sblock_input = model_input->get_child("subblocks");
      for (auto& b : *sblock_input) {
        array<int,2> charges = b->get_array<int, 2>("charges");
        array<int,2> spins = b->get_array<int, 2>("spins");
        const int nstates = b->get<int>("nstates");
        this_model.emplace_back( make_pair(spins[0],spins[1]), make_pair(charges[0],charges[1]), make_pair(p,p), nstates );
      }
      models_to_form_.emplace_back(move(this_model));
    }
  }

  Timer timer;

  shared_ptr<const Reference> dimerref = dimer_->sref();

  jop_ = make_shared<DimerJop>(dimerref, dimerref->nclosed(), dimerref->nclosed() + dimer_->active_refs().first->nact(), dimerref->nclosed() + dimerref->nact(), dimerref->coeff());
  cout << "  o computing integrals: " << timer.tick() << endl;

  energies_ = vector<double>(nstates_, 0.0);
}


shared_ptr<Matrix> ASD_base::compute_intra(const DimerSubspace_base& AB, shared_ptr<const DimerJop> jop, const double diag) const {
  auto out = make_shared<Matrix>(AB.dimerstates(), AB.dimerstates());

  const int nstatesA = AB.nstates<0>();
  const int nstatesB = AB.nstates<1>();

  // first H^{AA}_{AA}
  for(int stateA = 0; stateA < nstatesA; ++stateA) {
    for(int stateAp = 0; stateAp < stateA; ++stateAp) {
      const double value = AB.sigma<0>()->element(stateAp, stateA);
      for(int stateB = 0; stateB < nstatesB; ++stateB) {
        const int stateApB = AB.dimerindex(stateAp, stateB);
        const int stateAB = AB.dimerindex(stateA, stateB);
        (*out)(stateAB, stateApB) += value;
        (*out)(stateApB, stateAB) += value;
      }
    }
    const double value = AB.sigma<0>()->element(stateA, stateA);
    for(int stateB = 0; stateB < nstatesB; ++stateB) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      (*out)(stateAB,stateAB) += value;
    }
  }

  // H^{BB}_{BB}
  for(int stateB = 0; stateB < nstatesB; ++stateB) {
    for(int stateBp = 0; stateBp < stateB; ++stateBp) {
      const double value = AB.sigma<1>()->element(stateBp, stateB);
      for(int stateA = 0; stateA < nstatesA; ++stateA) {
        const int stateAB = AB.dimerindex(stateA, stateB);
        const int stateABp = AB.dimerindex(stateA, stateBp);
        (*out)(stateAB, stateABp) += value;
        (*out)(stateABp, stateAB) += value;
      }
    }
    const double value = AB.sigma<1>()->element(stateB, stateB);
    for(int stateA = 0; stateA < nstatesA; ++stateA) {
      const int stateAB = AB.dimerindex(stateA, stateB);
      (*out)(stateAB,stateAB) += value;
    }
  }

  out->add_diag(diag);
  return out;
}



shared_ptr<Matrix> ASD_base::compute_diagonal_block_H(const DimerSubspace_base& subspace) const {
  const double core = dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy();

  auto out = compute_intra(subspace, jop_, core);
  array<MonomerKey,4> keys {{ subspace.monomerkey<0>(), subspace.monomerkey<1>(), subspace.monomerkey<0>(), subspace.monomerkey<1>() }};
  *out += *compute_inter_2e_H(keys);

  return out;
}


//***************************************************************************************************************

tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_diagonal_block_RDM(const DimerSubspace_base& subspace) const {
// 1e is not considered (TODO maybe there is contribution?)
//***************************************************************************************************************

  array<MonomerKey,4> keys {{ subspace.monomerkey<0>(), subspace.monomerkey<1>(), subspace.monomerkey<0>(), subspace.monomerkey<1>() }};
  auto out = compute_inter_2e_RDM(keys, /*subspace diagonal*/true);

  return out;
}



shared_ptr<Matrix> ASD_base::compute_offdiagonal_1e_H(const array<MonomerKey,4>& keys, shared_ptr<const Matrix> hAB) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  Coupling term_type = coupling_type(keys);

  GammaSQ operatorA;
  GammaSQ operatorB;
  int neleA = Ap.nelea() + Ap.neleb();

  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());

  switch(term_type) {
    case Coupling::aET :
      operatorA = GammaSQ::CreateAlpha;
      operatorB = GammaSQ::AnnihilateAlpha;
      break;
    case Coupling::inv_aET :
      operatorA = GammaSQ::AnnihilateAlpha;
      operatorB = GammaSQ::CreateAlpha;
      --neleA;
      break;
    case Coupling::bET :
      operatorA = GammaSQ::CreateBeta;
      operatorB = GammaSQ::AnnihilateBeta;
      break;
    case Coupling::inv_bET :
      operatorA = GammaSQ::AnnihilateBeta;
      operatorB = GammaSQ::CreateBeta;
      --neleA;
      break;
    default :
      return out;
  }

  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {operatorA});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {operatorB});
  Matrix tmp = gamma_A * (*hAB) ^ gamma_B;

  if ((neleA % 2) == 1) {
    // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }

  return out;
}


//***************************************************************************************************************

//tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
//ASD_base::compute_offdiagonal_1e_RDM(const array<MonomerKey,4>& keys, shared_ptr<const Matrix> hAB) const {
// Unsed function!
//***************************************************************************************************************
//assert(false);
//return make_tuple(nullptr,nullptr);
//}


// This term will couple off-diagonal blocks since it has no delta functions involved

shared_ptr<Matrix> ASD_base::compute_inter_2e_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  // alpha-alpha
  auto gamma_AA_alpha = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
  auto gamma_BB_alpha = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});

  // beta-beta
  auto gamma_AA_beta = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
  auto gamma_BB_beta = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});

  // build J and K matrices
  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,1,0,1>();
  shared_ptr<const Matrix> Kmatrix = jop_->coulomb_matrix<0,1,1,0>();

  Matrix tmp((gamma_AA_alpha + gamma_AA_beta) * (*Jmatrix) ^ (gamma_BB_alpha + gamma_BB_beta));

  tmp -= gamma_AA_alpha * (*Kmatrix) ^ gamma_BB_alpha;
  tmp -= gamma_AA_beta * (*Kmatrix) ^ gamma_BB_beta;

  // sort: (A,A',B,B') --> (A,B,A',B') + block(A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  return out;
}


//***************************************************************************************************************

tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_inter_2e_RDM(const array<MonomerKey,4>& keys, const bool subdia) const {
//***************************************************************************************************************
  auto& B  = keys[1]; 
  auto& Bp = keys[3];

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;
  auto out = make_shared<RDM<2>>(nactA+nactB);

  // alpha-alpha
  auto gamma_AA_alpha = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
  auto gamma_BB_alpha = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});

  // beta-beta
  auto gamma_AA_beta = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
  auto gamma_BB_beta = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});

  auto rdmAA = make_shared<Matrix>(gamma_AA_alpha % gamma_BB_alpha);
  auto rdmBB = make_shared<Matrix>(gamma_AA_beta  % gamma_BB_beta);

  auto rdmAB = make_shared<Matrix>(gamma_AA_alpha % gamma_BB_beta);
  auto rdmBA = make_shared<Matrix>(gamma_AA_beta  % gamma_BB_alpha);

  {// P(p,q',r',s) : p15
    auto rdmt = rdmAA->clone();
    SMITH::sort_indices<0,3,2,1, 0,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa
    SMITH::sort_indices<0,3,2,1, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb
    if (!subdia) {
      SMITH::sort_indices<1,2,3,0, 1,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa of (N,M)
      SMITH::sort_indices<1,2,3,0, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
    }
    auto low = {    0, nactA, nactA,     0};
    auto up  = {nactA, nactT, nactT, nactA};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  {// d_pqr's' : p19
    auto rdmt = rdmAA->clone();
    SMITH::sort_indices<0,1,2,3, 0,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa
    SMITH::sort_indices<0,1,2,3, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb
    SMITH::sort_indices<0,1,2,3, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa bb
    SMITH::sort_indices<0,1,2,3, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb aa
    if (!subdia) {
      SMITH::sort_indices<1,0,3,2, 1,1, 1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa of (N,M)
      SMITH::sort_indices<1,0,3,2, 1,1, 1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
      SMITH::sort_indices<1,0,3,2, 1,1, 1,1>(rdmAB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
      SMITH::sort_indices<1,0,3,2, 1,1, 1,1>(rdmBA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)
    }
    auto low = {    0,     0, nactA, nactA};
    auto up  = {nactA, nactA, nactT, nactT};
    auto outv = make_rwview(out->range().slice(low, up), out->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(nullptr,out);
}



shared_ptr<Matrix> ASD_base::compute_aET_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  Matrix tmp(A.nstates()*Ap.nstates(), B.nstates()*Bp.nstates());

  // One-body aET
  {
    auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha});
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    shared_ptr<const Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }

  //Two-body aET, type 1
  {
    auto gamma_A  = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha});
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    auto gamma_A1 = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_A2 = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta,  GammaSQ::AnnihilateBeta});
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  const int neleA = Ap.nelea() + Ap.neleb();
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  if ((neleA % 2) == 1) {
    // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  } else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  return out;
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_aET_RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************

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
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha});
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});
    
    auto rdm = make_shared<Matrix>(gamma_A % gamma_B);
    auto rdmt = rdm->clone();

    // P(p,q') : p10
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);
 
    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  //2RDM
  {
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha});
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1);
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2);
    auto rdmt = rdm1->clone();

    // P(p,q',r',s') : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    SMITH::sort_indices<0,2,1,3, 1,1, -1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }
  {
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha});

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B);
    auto rdmt = rdm1->clone();

    //P(p,q',r,s) : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    SMITH::sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA,     0,     0};
    auto up  = {nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }
  
  return make_tuple(out1,out2);
}



shared_ptr<Matrix> ASD_base::compute_bET_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  Matrix tmp(A.nstates()*Ap.nstates(), B.nstates()*Bp.nstates());

  // One-body bET
  {
    auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta});
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += gamma_A * (*Fmatrix) ^ gamma_B;
  }


  //Two-body bET, type 1
  {
    auto gamma_A  = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta});
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,1,1,1>();

    tmp -= gamma_A * (*Jmatrix) ^ (gamma_B1 + gamma_B2);
  }

  //Two-body aET, type 2
  {
    auto gamma_A1 = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_A2 = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,0>();

    tmp += (gamma_A1 + gamma_A2) * (*Jmatrix) ^ gamma_B;
  }

  const int neleA = Ap.nelea() + Ap.neleb();
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  if ((neleA % 2) == 1) {
    // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    SMITH::sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }

  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_bET_RDM(const array<MonomerKey,4>& keys) const {
//***************************************************************************************************************
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
    auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta});
    auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});
    
    auto rdm = make_shared<Matrix>(gamma_A % gamma_B);
    auto rdmt = rdm->clone();
    
    // P(p,q') : p10
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,1, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactB);
    rdmt->scale(fac);
    
    auto low = {    0, nactA};
    auto up  = {nactA, nactT};
    auto outv = make_rwview(out1->range().slice(low, up), out1->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  //RDM2
  {
    auto gamma_A  = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta});
    auto gamma_B1 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(gamma_A % gamma_B1);
    auto rdm2 = make_shared<Matrix>(gamma_A % gamma_B2);
    auto rdmt = rdm1->clone();

    // P(p,q',r',s') : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,2,1,3, 0,1, -1,1>(rdm1->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    SMITH::sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactB, nactB, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA, nactA, nactA};
    auto up  = {nactA, nactT, nactT, nactT};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }
  {
    auto gamma_A1 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha, GammaSQ::AnnihilateAlpha});
    auto gamma_A2 = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta, GammaSQ::AnnihilateBeta});
    auto gamma_B  = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta});

    auto rdm1 = make_shared<Matrix>(gamma_A1 % gamma_B);
    auto rdm2 = make_shared<Matrix>(gamma_A2 % gamma_B);
    auto rdmt = rdm1->clone();

    // P(p,q',r,s) : p15
    int fac = {neleA%2 == 0 ? 1 : -1};
    SMITH::sort_indices<0,3,1,2, 0,1,  1,1>(rdm1->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    SMITH::sort_indices<0,3,1,2, 1,1,  1,1>(rdm2->data(), rdmt->data(), nactA, nactA, nactA, nactB);
    rdmt->scale(fac);

    auto low = {    0, nactA,     0,     0};
    auto up  = {nactA, nactT, nactA, nactA};
    auto outv = make_rwview(out2->range().slice(low, up), out2->storage());
    copy(rdmt->begin(), rdmt->end(), outv.begin());
  }

  return make_tuple(out1,out2);
}



shared_ptr<Matrix> ASD_base::compute_abFlip_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Kmatrix = jop_->coulomb_matrix<0,1,1,0>();

  Matrix tmp = gamma_A * (*Kmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_abFlip_RDM(const array<MonomerKey,4>& keys) const {
// if ab-flip, account ba-flip arising from (N,M)
// if(M,N) is ba-flip then (N,M) is ab-flip and this will include ba-flip of (M,N) too.
//***************************************************************************************************************
  auto& B = keys[1];
  auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha}));

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  auto out = make_shared<RDM<2>>(nactA+nactB);

  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm = make_shared<Matrix>(gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  // P(p,q',r',s) : p4  ab-flip of (M,N)
  SMITH::sort_indices<0,3,2,1, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //ab-flip
  SMITH::sort_indices<1,2,3,0, 1,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB); //ba-flip of (N,M) p15B

  auto low = {    0, nactA, nactA,     0};
  auto up  = {nactA, nactT, nactT, nactA};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr,out);
}



shared_ptr<Matrix> ASD_base::compute_abET_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_abET_RDM(const array<MonomerKey,4>& keys) const {
// for (M,N)
// if inverse ab-ET / compute (N,M)
//***************************************************************************************************************
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateBeta}));

  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p14B
  SMITH::sort_indices<0,2,1,3, 0,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);
  SMITH::sort_indices<1,3,0,2, 1,1, -1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr,out);
}



shared_ptr<Matrix> ASD_base::compute_aaET_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<0,2,1,3,0,1,-1,2>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_aaET_RDM(const array<MonomerKey,4>& keys) const {
//off-diagonal subspaces only!
// if(M,N) is inverse-aa-ET, swap M,N as (N,M) will be aa-ET and contribute to 2RDM
//***************************************************************************************************************
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha}));

  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  auto rdm  = make_shared<Matrix>(gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p1B
  SMITH::sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr,out);
}



shared_ptr<Matrix> ASD_base::compute_bbET_H(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<0,2,1,3,0,1,-1,2>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> 
ASD_base::compute_bbET_RDM(const array<MonomerKey,4>& keys) const {
// cf. aaET
//***************************************************************************************************************
  auto& B = keys[1]; auto& Bp = keys[3];

  assert(gammatensor_[0]->exist(keys[0], keys[2], {GammaSQ::CreateBeta, GammaSQ::CreateBeta}));

  auto gamma_A = worktensor_->get_block_as_matview(B, Bp, {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  auto rdm  = make_shared<Matrix>(gamma_A % gamma_B);
  auto rdmt = rdm->clone();

  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  // P(p,q',r,s') : p14B
  SMITH::sort_indices<0,3,1,2, 0,1, 1,1>(rdm->data(), rdmt->data(), nactA, nactA, nactB, nactB);

  auto out = make_shared<RDM<2>>(nactA+nactB);
  auto low = {    0, nactA,     0, nactA};
  auto up  = {nactA, nactT, nactA, nactT};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  assert(rdmt->size() == outv.size());
  copy(rdmt->begin(), rdmt->end(), outv.begin());

  return make_tuple(nullptr,out);
}


void ASD_base::print_hamiltonian(const string title, const int nstates) const {
  hamiltonian_->print(title, nstates);
}


void ASD_base::print_states(const Matrix& cc, const vector<double>& energies, const double thresh, const string title) const {
  const int nstates = cc.mdim();
  shared_ptr<Matrix> spn = spin_->apply(cc);
  cout << endl << " ===== " << title << " =====" << endl;
  for (int istate = 0; istate < nstates; ++istate) {
    cout << "   state  " << setw(3) << istate << ": "
         << setprecision(8) << setw(17) << fixed << energies.at(istate)
         << "   <S^2> = " << setw(4) << setprecision(4) << fixed << ddot_(dimerstates_, spn->element_ptr(0,istate), 1, cc.element_ptr(0,istate), 1) << endl;
    const double *eigendata = cc.element_ptr(0,istate);
    double printed = 0.0;
    for (auto& subspace : subspaces_base()) {
      const int nA = subspace.nstates<0>();
      const int nB = subspace.nstates<1>();
      for (int i = 0; i < nA; ++i) {
        for (int j = 0; j < nB; ++j, ++eigendata) {
          if ( (*eigendata)*(*eigendata) > thresh ) {
            cout << "      " << subspace.string(i,j) << setprecision(12) << setw(20) << *eigendata << endl;
            printed += (*eigendata)*(*eigendata);
          }
        }
      }
    }
    cout << "    total weight of printed elements: " << setprecision(12) << setw(20) << printed << endl << endl;
  }
}


void ASD_base::print_property(const string label, shared_ptr<const Matrix> property , const int nstates) const {
  const string indent("   ");
  const int nprint = min(nstates, property->ndim());

  cout << indent << " " << label << "    |0>";
  for (int istate = 1; istate < nprint; ++istate) cout << "         |" << istate << ">";
  cout << endl;
  for (int istate = 0; istate < nprint; ++istate) {
    cout << indent << "<" << istate << "|";
    for (int jstate = 0; jstate < nprint; ++jstate) {
      cout << setw(12) << setprecision(6) << property->element(jstate, istate);
    }
    cout << endl;
  }
  cout << endl;
}


void ASD_base::print(const double thresh) const {
  print_states(*adiabats_, energies_, thresh, "Adiabatic States");
  if (dipoles_) {for (auto& prop : properties_) print_property(prop.first, prop.second, nstates_); }
}

//***************************************************************************************************************
shared_ptr<Matrix>
ASD_base::couple_blocks_H(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {
//***************************************************************************************************************

  Coupling term_type = coupling_type(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  shared_ptr<Matrix> out;
  std::array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = nullptr; break;
    case Coupling::diagonal :
      out = compute_inter_2e_H(keys); break;
    case Coupling::aET :
      out = compute_aET_H(keys); break;
    case Coupling::bET :
      out = compute_bET_H(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip_H(keys); break;
    case Coupling::abET :
      out = compute_abET_H(keys); break;
    case Coupling::aaET :
      out = compute_aaET_H(keys); break;
    case Coupling::bbET :
      out = compute_bbET_H(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }

  /* if we are computing the Hamiltonian and flip = true, then we tranpose the output (see above) */
//if (flip) out = out->transpose(); // same as below??
  if (flip) transpose_call(out);

  return out;
}

//***************************************************************************************************************
tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>>
ASD_base::couple_blocks_RDM(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {
//***************************************************************************************************************

  Coupling term_type = coupling_type(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }
  
  tuple<shared_ptr<RDM<1>>,shared_ptr<RDM<2>>> out;
  std::array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = make_tuple(nullptr,nullptr); break;
    case Coupling::diagonal :
      out = compute_inter_2e_RDM(keys, /*subspace diagonal*/false); break;
    case Coupling::aET :
      out = compute_aET_RDM(keys); break;
    case Coupling::bET :
      out = compute_bET_RDM(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip_RDM(keys); break;
    case Coupling::abET :
      out = compute_abET_RDM(keys); break;
    case Coupling::aaET :
      out = compute_aaET_RDM(keys); break;
    case Coupling::bbET :
      out = compute_bbET_RDM(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }
  
  return out;
}

//***************************************************************************************************************
void
ASD_base::debug_RDM() const {
//***************************************************************************************************************
  const int nactA = dimer_->embedded_refs().first->nact();
  const int nactB = dimer_->embedded_refs().second->nact();
  const int nactT = nactA+nactB;

  const int neleA = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int neleB = 2*(dimer_->isolated_refs().first->nclosed() - dimer_->active_refs().first->nclosed());
  const int nelec = neleA+neleB;

  cout << "#of active electrons is A : " << neleA << endl;
  cout << "#of active electrons is B : " << neleB << endl;
  cout << "#of total active electrons: " << nelec << endl;

  auto rdm1A = std::make_shared<RDM<1>>(nactA);
  {
    auto low = {0,0};
    auto up  = {nactA,nactA};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage());
    copy(view.begin(), view.end(), rdm1A->begin());
  }
  auto rdm1B = std::make_shared<RDM<1>>(nactB);
  {
    auto low = {nactA,nactA};
    auto up  = {nactT,nactT};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage());
    copy(view.begin(), view.end(), rdm1B->begin());
  }

  auto rdm2A = std::make_shared<RDM<2>>(nactA);
  {
    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage());
    copy(view.begin(), view.end(), rdm2A->begin());
  }
  auto rdm2B = std::make_shared<RDM<2>>(nactB);
  {
    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage());
    copy(view.begin(), view.end(), rdm2B->begin());
  }

//auto rdm3A = std::make_shared<RDM<3>>(nactA);
//{
//  auto low = {0,0,0,0,0,0};
//  auto up  = {nactA,nactA,nactA,nactA,nactA,nactA};
//  auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
//  copy(view.begin(), view.end(), rdm3A->begin());
//}
//auto rdm3B = std::make_shared<RDM<3>>(nactB);
//{
//  auto low = {nactA,nactA,nactA,nactA,nactA,nactA};
//  auto up  = {nactT,nactT,nactT,nactT,nactT,nactT};
//  auto view = btas::make_view(threerdm_->range().slice(low,up), threerdm_->storage());
//  copy(view.begin(), view.end(), rdm3B->begin());
//}

  //1RDM
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(A)  Trace = " << sum << std::endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(B)  Trace = " << sum << std::endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i) {
      sum += onerdm_->element(i,i);
    }
    std::cout << "1RDM(AB) Trace = " << sum << std::endl;
  }

  //2RDM check (A) Gamma_ij,kl = <0|E_ij,kl|0> = <0|(k1)'(i2)'(j2)(l1)|0>  1,2 = spin 
  //diagonal: i=j, k=l
  { //Monomer A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(A)  Trace = " << sum << std::endl;
  }
  { //Monomer B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(B)  Trace = " << sum << std::endl;
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += twordm_->element(i,i,j,j);
    }
    std::cout << "2RDM(AB) Trace = " << sum << std::endl;
  }

  { //Gamma_ij,kk
    std::cout << "2RDM(A) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1A);
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      debug->element(i,j) -= 1.0/(neleA-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }
  { //Gamma_ij,kk
    std::cout << "2RDM(B) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1B);
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      debug->element(i-nactA,j-nactA) -= 1.0/(neleB-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }
  { //Gamma_ij,kk
    std::cout << "2RDM(AB) Partial Trace Sum_k (i,j,k,k)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*onerdm_);
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      debug->element(i,j) -= 1.0/(nelec-1) * twordm_->element(i,j,k,k);
    }
    debug->print(1.0e-3);
  }


  //construct approx twordm
  cout << "Build approx 2RDM:" << endl;
  for (int i = 0; i != nactA; ++i)
  for (int j = nactA; j != nactT; ++j) {
    approx2rdm_->element(i,i,j,j) = onerdm_->element(i,i) * onerdm_->element(j,j);
    approx2rdm_->element(j,j,i,i) = approx2rdm_->element(i,i,j,j);
    approx2rdm_->element(i,j,j,i) = -0.5*(onerdm_->element(i,i) * onerdm_->element(j,j));// - onerdm_->element(i,i);
    approx2rdm_->element(j,i,i,j) = approx2rdm_->element(i,j,j,i); //-onerdm_->element(i,i) * onerdm_->element(j,j);// - onerdm_->element(j,j);
  }
  { //Dimer AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j) {
      sum += approx2rdm_->element(i,i,j,j);
    }
    std::cout << "APPROX2RDM(AB) Trace = " << sum << std::endl;
  }
  { //difference
    auto debug = make_shared<RDM<2>>(*twordm_);
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

  assert(false);

  //3RDM: Gamma_ij,kl,mn
  { //Trace A
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (A)  = " << sum << std::endl;
  }
  { //Trace B
    double sum = 0.0;
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (B)  = " << sum << std::endl;
  }
  { //Trace AB
    double sum = 0.0;
    for (int i = 0; i != nactT; ++i)
    for (int j = 0; j != nactT; ++j)
    for (int k = 0; k != nactT; ++k) {
      sum += threerdm_->element(i,i,j,j,k,k);
    }
    std::cout << "3RDM Trace (AB) = " << sum << std::endl;
  }

  { //Gamma_ij,kl,mm : p21
    std::cout << "3RDM(A) Partial Trace Sum_m (i,j,k,l,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<2>>(*rdm2A);
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) 
    for (int l = 0; l != nactA; ++l) 
    for (int m = 0; m != nactA; ++m) {
      debug->element(i,j,k,l) -= 1.0/(neleA-2) * threerdm_->element(i,j,k,l,m,m);
    }
    debug->print(1.0e-8);
  }
  { //Gamma_ij,kk,mm : p21
    std::cout << "3RDM(A) Partial Trace Sum_m (i,j,k,k,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1A);
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k) 
    for (int m = 0; m != nactA; ++m) {
      debug->element(i,j) -= 1.0/((neleA-2)*(neleA-1)) * threerdm_->element(i,j,k,k,m,m);
    }
    debug->print(1.0e-8);
  }

  { //Gamma_ij,kl,mm : p21
    std::cout << "3RDM(B) Partial Trace Sum_m (i,j,k,l,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<2>>(*rdm2B);
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) 
    for (int l = nactA; l != nactT; ++l) 
    for (int m = nactA; m != nactT; ++m) {
      debug->element(i-nactA,j-nactA,k-nactA,l-nactA) -= 1.0/(neleA-2) * threerdm_->element(i,j,k,l,m,m);
    }
    debug->print(1.0e-8);
  }
  { //Gamma_ij,kk,mm : p21
    std::cout << "3RDM(B) Partial Trace Sum_m (i,j,k,k,m,m)" << std::endl;
    auto debug = std::make_shared<RDM<1>>(*rdm1B);
    for (int i = nactA; i != nactT; ++i)
    for (int j = nactA; j != nactT; ++j)
    for (int k = nactA; k != nactT; ++k) 
    for (int m = nactA; m != nactT; ++m) {
      debug->element(i-nactA,j-nactA) -= 1.0/((neleA-2)*(neleA-1)) * threerdm_->element(i,j,k,k,m,m);
    }
    debug->print(1.0e-8);
  }



  assert(false);
}

#if 0

  //4RDM check (A) Gamma_ij,kl,mn,op
  {
    double sum = 0.0;
    for (int i = 0; i != nactA; ++i)
    for (int j = 0; j != nactA; ++j)
    for (int k = 0; k != nactA; ++k)
    for (int l = 0; l != nactA; ++l) {
      sum += fourrdm_->element(i,i,j,j,k,k,l,l);
    }
    std::cout << "4RDM Trace = " << sum << std::endl;
  }



  // Checking 4RDM by comparing with 3RDM
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 1" << std::endl;
    for (int l = 0; l != nactA; ++l)
      for (int d = 0; d != nactA; ++d)
        for (int k = 0; k != nactA; ++k)
          for (int c = 0; c != nactA; ++c)
            for (int j = 0; j != nactA; ++j)
              for (int b = 0; b != nactA; ++b)
      for (int i = 0; i != nactA; ++i) {
        debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(i,i,b,j,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,i,i,c,k,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,i,i,d,l);
  //    debug->element(b,j,c,k,d,l) -= 1.0/(nelea+neleb-3) * fourrdm_->element(b,j,c,k,d,l,i,i);
      }
    debug->print(1.0e-8);
  }
  { 
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 2" << std::endl;
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
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 3" << std::endl;
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
    auto debug = std::make_shared<RDM<3>>(*threerdm_);
    std::cout << "4RDM debug test 4" << std::endl;
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

  std::cout << "Monomer A 4RDM print" << std::endl;
  for (int i = 0; i != nactA; ++ i)
  for (int j = 0; j != nactA; ++ j)
  for (int k = 0; k != nactA; ++ k)
  for (int l = 0; l != nactA; ++ l)
  for (int m = 0; m != nactA; ++ m)
  for (int n = 0; n != nactA; ++ n)
  for (int o = 0; o != nactA; ++ o)
  for (int p = 0; p != nactA; ++ p) {
    double elem = fourrdm_->element(i,j,k,l,m,n,o,p);
    elem = std::abs(elem);
    if(elem > 1.0e-8) std::cout << "RDM4(" << i << j << k << l << m << n << o << p << ") " << fourrdm_->element(i,j,k,l,m,n,o,p) << std::endl ;
  }
  assert(false);
#endif

//***************************************************************************************************************
void
ASD_base::debug_energy() const {
//***************************************************************************************************************

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

  auto rdm1 = onerdm_->rdm1_mat(0);
  rdm1->print("1RDM",nactT);

  double  e1 = ddot_(nactT*nactT, int1->element_ptr(0,0), 1, rdm1->element_ptr(0,0), 1);
  cout << "1E energy = " << e1 << endl;

  double e2 = 0.0;
  //AAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,0,0>();
    auto int2 = make_shared<Matrix>(nactA*nactA*nactA*nactA,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactA); //conver to chemist not.

    auto low = {0,0,0,0};
    auto up  = {nactA,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA,nactA,1); //empty d_AAAA (note: the dimension specification actually do not matter)
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAA filled
    cout << "2E energy (AAAA) = " << 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BBBB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,1>();
    auto int2 = make_shared<Matrix>(1,nactB*nactB*nactB*nactB);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactB, nactB, nactB, nactB); //conver to chemist not.

    auto low = {nactA,nactA,nactA,nactA};
    auto up  = {nactT,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBB sector of d
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
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.
    auto low = {    0,    0,    0,nactA};
    auto up  = {nactA,nactA,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactA*nactB,1); //empty d_AAAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AAAB filled
    cout << "2E energy (AAAB) = " << 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //AABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,0,0>(); // <pq'|rs> in (prs,q') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1);
    SMITH::sort_indices<0,1,3,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [ij|k'l]

    auto low = {    0,    0,nactA,    0};
    auto up  = {nactA,nactA,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactA,1); //empty d_AABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABA filled
    cout << "2E energy (AABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //BAAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,0>(); // <p'q|rs> in (qrs,p') format 
    auto int2 = make_shared<Matrix>(nactB,nactA*nactA*nactA);
    SMITH::sort_indices<3,1,0,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not.

    auto low = {nactA,    0,    0,    0};
    auto up  = {nactT,nactA,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAA sector of d
    auto rdm2 = make_shared<Matrix>(nactB,nactA*nactA*nactA); //empty d_BAAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BAAA filled
    cout << "2E energy (BAAA) = " << 0.5 * ddot_(nactB*nactA*nactA*nactA, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactA*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,0>(); // <pq|r's> in (pqs,r') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactA*nactA,1);
    SMITH::sort_indices<0,3,1,2, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactA, nactB); //conver to chemist not. [pr'|qs]=[ij'|kl]

    auto low = {    0,nactA,    0,    0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAA sector of d
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
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [pr'|q's']=[ij'|k'l']

    auto low = {    0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_ABBB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (ABBB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,1>(); // <p'q'|rs'> in (r,p'q's') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,0,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r|q's']=[i'j|k'l']

    auto low = {nactA,    0,nactA,nactA};
    auto up  = {nactT,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BABB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,1>(); // <p'q|r's'> in (q,p'r's') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,2,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|qs']=[i'j'|kl']

    auto low = {nactA,nactA,    0,nactA};
    auto up  = {nactT,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1); //empty d_BBAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBB filled
    cout << "2E energy (BBAB) = " << 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactB*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,1,0>(); // <p'q'|r's> in (s,p'q'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactB*nactB*nactB,1);
    SMITH::sort_indices<1,3,2,0, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactB, nactB, nactB); //conver to chemist not. [p'r'|q's]=[i'j'|k'l]

    auto low = {nactA,nactA,nactA,    0};
    auto up  = {nactT,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBA sector of d
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
    SMITH::sort_indices<0,1,2,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr|q's']=[ij|k'l']

    auto low = {    0,    0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), rdm2->begin()); //d_AABB filled
    cout << "2E energy (AABB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BBAA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,1,0>(); // <p'q|r's> in (qs,p'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r'|qs]=[i'j'|kl]

    auto low = {nactA,nactA,    0,    0};
    auto up  = {nactT,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BBAA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BBAA filled
    cout << "2E energy (BBAA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,0,1,1>(); // <pq|r's'> in (pq,r's') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<0,2,1,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|qs']=[ij'|kl']

    auto low = {    0,nactA,    0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAB sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABAB
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABAB filled
    cout << "2E energy (ABAB) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BABA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,1,0,0>(); // <p'q'|rs> in (rs,p'q') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,0,3,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|q's]=[i'j|k'l]

    auto low = {nactA,    0,nactA,    0};
    auto up  = {nactT,nactA,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BABA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_BABA
    copy(view.begin(), view.end(), rdm2->begin()); //d_BABA filled
    cout << "2E energy (BABA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }

  //ABBA
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<0,1,1,0>(); // <pq'|r's> in (ps,q'r') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<0,3,2,1, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [pr'|q's]=[ij'|k'l]

    auto low = {    0,nactA,nactA,    0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBA sector of d
    auto rdm2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_ABBA
    copy(view.begin(), view.end(), rdm2->begin()); //d_ABBA filled
    cout << "2E energy (ABBA) = " << 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1) << endl;
    e2 += 0.5 * ddot_(nactA*nactA*nactB*nactB, int2->element_ptr(0,0), 1, rdm2->element_ptr(0,0), 1);
  }
  //BAAB
  {
    shared_ptr<const Matrix> pint2 = jop_->coulomb_matrix<1,0,0,1>(); // <p'q|rs'> in (qr,p's') format 
    auto int2 = make_shared<Matrix>(nactA*nactA*nactB*nactB,1);
    SMITH::sort_indices<2,1,0,3, 0,1, 1,1>(pint2->data(), int2->data(), nactA, nactA, nactB, nactB); //conver to chemist not. [p'r|qs']=[i'j|kl']

    auto low = {nactA,    0,    0,nactA};
    auto up  = {nactT,nactA,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAB sector of d
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

//***************************************************************************************************************
void
ASD_base::symmetrize_RDM() const {
//***************************************************************************************************************

  cout << "!@# Unsymmetrized 1RDM" << endl;
  onerdm_->print(1.0e-6);

  const int nactA = dimer_->active_refs().first->nact();
  const int nactB = dimer_->active_refs().second->nact();
  const int nactT = nactA + nactB;  

  //Symmetrize: D_AB (calculated) D_BA (uncalc.& symmetrized here)
  auto matBA = std::make_shared<Matrix>(nactB,nactA); //D_BA empty
  {
    auto low = {0, nactA};
    auto up  = {nactA, nactT};
    auto view = btas::make_view(onerdm_->range().slice(low,up), onerdm_->storage()); //D_AB sector of D (read ptr)
    auto matAB = std::make_shared<Matrix>(nactA,nactB); //D_AB empty
    std::copy(view.begin(), view.end(), matAB->begin()); //D_AB filled
    SMITH::sort_indices<1,0, 0,1, 1,1>(matAB->data(), matBA->data(), nactA, nactB); // transpose and fill D_BA
  }
  {
    auto low = {nactA, 0};
    auto up  = {nactT, nactA};
    auto outv = btas::make_rwview(onerdm_->range().slice(low,up), onerdm_->storage()); //D_BA sector of D (read & write ptr)
    std::copy(matBA->begin(), matBA->end(), outv.begin()); //copy D_BA -> D_BA sector of D
  }

  cout << "!@# Symmetrized 1RDM" << endl;
  onerdm_->print(1.0e-6);

  //Symmetrize: d(ABAA) note p18B, 19B
  {
    auto low = {0,nactA,0,0};
    auto up  = {nactA,nactT,nactA,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactA); //empty d_ABAA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAA filled
    { //d(AAAB)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactA*nactB); //empty d_AAAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AAAB
      auto low = {0,0,0,nactA};
      auto up  = {nactA,nactA,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AAAB into d_AAAB sector of d
    } 
    { //d(BAAA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactA); //empty d_BAAA
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_BAAA
      auto low = {nactA,0,0,0};
      auto up  = {nactT,nactA,nactA,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAA into d_BAAA sector of d
    } 
    { //d(AABA)
      auto outmat = make_shared<Matrix>(nactA*nactA,nactB*nactA); //empty d_AABA
      SMITH::sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactA); //reorder and fill d_AABA
      auto low = {0,0,nactA,0};
      auto up  = {nactA,nactA,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_AABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_AABA into d_AABA sector of d
    } 
  }
 
  //Symmetrize: d(ABBB) note p18B, 19B
  {
    auto low = {0,nactA,nactA,nactA};
    auto up  = {nactA,nactT,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactB); //empty d_ABBB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBB filled
    { //d(BBAB)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactA*nactB); //empty d_BBAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBAB
      auto low = {nactA,nactA,0,nactA};
      auto up  = {nactT,nactT,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAB into d_BBAB sector of d
    } 
    { //d(BABB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactB); //empty d_BABB
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BABB
      auto low = {nactA,0,nactA,nactA};
      auto up  = {nactT,nactA,nactT,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BABB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BABB into d_BABB sector of d
    } 
    { //d(BBBA)
      auto outmat = make_shared<Matrix>(nactB*nactB,nactB*nactA); //empty d_BBBA
      SMITH::sort_indices<3,2,1,0, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactB); //reorder and fill d_BBBA
      auto low = {nactA,nactA,nactA,0};
      auto up  = {nactT,nactT,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBBA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BBBA sector of d
    } 
  }


  //Symmetrize: d(ABBA) note p19
  {
    auto low = {0,nactA,nactA,0};
    auto up  = {nactA,nactT,nactT,nactA};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABBA sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactB*nactA); //empty d_ABBA
    copy(view.begin(), view.end(), inmat->begin()); //d_ABBA filled
    { //d(BAAB)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactA*nactB); //empty d_BAAB
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactB, nactA); //reorder and fill d_BAAB
      auto low = {nactA,0,0,nactA};
      auto up  = {nactT,nactA,nactA,nactT};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BAAB sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BAAB into d_BAAB sector of d
    } 
  }

  //Symmetrize: d(AABB) note 19
  { //d(AABB)
    auto low = {0,0,nactA,nactA};
    auto up  = {nactA,nactA,nactT,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_AABB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactA*nactB*nactB,1); //empty d_AABB
    copy(view.begin(), view.end(), inmat->begin()); //d_AABB filled
    { //d(BBAA)
      auto outmat = make_shared<Matrix>(nactB*nactB*nactA*nactA,1); //empty d_BBAA
      SMITH::sort_indices<2,3,0,1, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactA, nactB, nactB); //reorder and fill d_BBAA
      auto low = {nactA,nactA,0,0};
      auto up  = {nactT,nactT,nactA,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BBAA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBAA into d_BBAA sector of d
    } 
  }

  //Symmetrize: d(ABAB) note p19
  {
    auto low = {0,nactA,0,nactA};
    auto up  = {nactA,nactT,nactA,nactT};
    auto view = btas::make_view(twordm_->range().slice(low,up), twordm_->storage()); //d_ABAB sector of d
    auto inmat = make_shared<Matrix>(nactA*nactB,nactA*nactB); //empty d_ABAB
    copy(view.begin(), view.end(), inmat->begin()); //d_ABAB filled
    { //d(BABA)
      auto outmat = make_shared<Matrix>(nactB*nactA,nactB*nactA); //empty d_BABA
      SMITH::sort_indices<1,0,3,2, 0,1, 1,1>(inmat->data(), outmat->data(), nactA, nactB, nactA, nactB); //reorder and fill d_BABA
      auto low = {nactA,0,nactA,0};
      auto up  = {nactT,nactA,nactT,nactA};
      auto outv = btas::make_rwview(twordm_->range().slice(low,up), twordm_->storage()); //d_BABA sector of d
      copy(outmat->begin(), outmat->end(), outv.begin()); //copy d_BBBA into d_BABA sector of d
    } 
  }

}

//***************************************************************************************************************
tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>>
ASD_base::couple_blocks_RDM34(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {
//***************************************************************************************************************

  Coupling term_type = coupling_type_RDM34(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }
  
  tuple<shared_ptr<RDM<3>>,shared_ptr<RDM<4>>> out;
//std::array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = make_tuple(nullptr,nullptr); break;
    case Coupling::diagonal :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_inter_2e_RDM(keys, /*subspace diagonal*/false); break;
    case Coupling::aET :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_aET_RDM(keys); break;
    case Coupling::bET :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_bET_RDM(keys); break;
    case Coupling::abFlip :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_abFlip_RDM(keys); break;
    case Coupling::abET :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_abET_RDM(keys); break;
    case Coupling::aaET :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_aaET_RDM(keys); break;
    case Coupling::bbET :
      out = make_tuple(nullptr,nullptr); break;
    //out = compute_bbET_RDM(keys); break;
    default :
      throw std::logic_error("Asking for a coupling type that has not been written.");
  }
  
  return out;
}

