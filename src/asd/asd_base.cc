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
  auto out = compute_inter_2e_RDM(keys);

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
ASD_base::compute_inter_2e_RDM(const array<MonomerKey,4>& keys) const {
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

  // P(p,q',r',s) : p15
  auto rdmt = rdmAA->clone();
  SMITH::sort_indices<0,3,2,1, 0,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa
  SMITH::sort_indices<1,2,3,0, 1,1, -1,1>(rdmAA->data(), rdmt->data(), nactA, nactA, nactB, nactB); //aa of (N,M)
  SMITH::sort_indices<0,3,2,1, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb
  SMITH::sort_indices<1,2,3,0, 1,1, -1,1>(rdmBB->data(), rdmt->data(), nactA, nactA, nactB, nactB); //bb of (N,M)

  auto low = {    0, nactA, nactA,     0};
  auto up  = {nactA, nactT, nactT, nactA};
  auto outv = make_rwview(out->range().slice(low, up), out->storage());
  copy(rdmt->begin(), rdmt->end(), outv.begin());


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
      out = compute_inter_2e_RDM(keys); break;
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
