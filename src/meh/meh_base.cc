//
// BAGEL - Parallel electron correlation program.
// Filename: meh_base.cc
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

#include <src/meh/meh_base.h>
#include <src/smith/prim_op.h>

using namespace std;
using namespace bagel;

MEH_base::MEH_base(const shared_ptr<const PTree> input, shared_ptr<const Dimer> dimer) : dimer_(dimer)
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


// This term will couple off-diagonal blocks since it has no delta functions involved
template <>
shared_ptr<Matrix> MEH_base::compute_inter_2e<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  // alpha-alpha
  auto gamma_AA_alpha = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
  auto gamma_BB_alpha = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});

  // beta-beta
  auto gamma_AA_beta = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});
  auto gamma_BB_beta = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});

  // build J and K matrices
  shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,1,0,1>();
  shared_ptr<const Matrix> Kmatrix = jop_->template coulomb_matrix<0,1,1,0>();

  Matrix tmp((*gamma_AA_alpha + *gamma_AA_beta) * (*Jmatrix) ^ (*gamma_BB_alpha + *gamma_BB_beta));

  tmp -= *gamma_AA_alpha * (*Kmatrix) ^ *gamma_BB_alpha;
  tmp -= *gamma_AA_beta * (*Kmatrix) ^ *gamma_BB_beta;

  // sort: (A',A,B',B) --> (A,B,A',B') + block(A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());
  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_aET<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  Matrix tmp(A.nstates()*Ap.nstates(), B.nstates()*Bp.nstates());

  // One-body aET
  {
    auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateAlpha});
    auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha});

    shared_ptr<const Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += *gamma_A * (*Fmatrix) ^ *gamma_B;
  }

  //Two-body aET, type 1
  {
    auto gamma_A  = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateAlpha});
    auto gamma_B1 = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,1,1,1>();

    tmp -= *gamma_A * (*Jmatrix) ^ (*gamma_B1 + *gamma_B2);
  }

  //Two-body aET, type 2
  {
    auto gamma_A1 = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
    auto gamma_A2 = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateAlpha});
    auto gamma_B  = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha});

    shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,0,1,0>();

    tmp += (*gamma_A1 + *gamma_A2) * (*Jmatrix) ^ *gamma_B;
  }

  const int neleA = A.nelea() + A.neleb();
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  if ((neleA % 2) == 1) {
    // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());
  }
  else {
    // sort: (A',A,B',B) --> (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());
  }
  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_bET<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  Matrix tmp(A.nstates()*Ap.nstates(), B.nstates()*Bp.nstates());

  // One-body bET
  {
    auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateBeta});
    auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Fmatrix = jop_->cross_mo1e();

    tmp += *gamma_A * (*Fmatrix) ^ *gamma_B;
  }


  //Two-body bET, type 1
  {
    auto gamma_A  = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateBeta});
    auto gamma_B1 = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha});
    auto gamma_B2 = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,1,1,1>();

    tmp -= *gamma_A * (*Jmatrix) ^ (*gamma_B1 + *gamma_B2);
  }

  //Two-body aET, type 2
  {
    auto gamma_A1 = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateAlpha, GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
    auto gamma_A2 = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateBeta, GammaSQ::CreateBeta, GammaSQ::CreateBeta});
    auto gamma_B  = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta});

    shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,0,1,0>();

    tmp += (*gamma_A1 + *gamma_A2) * (*Jmatrix) ^ *gamma_B;
  }

  const int neleA = A.nelea() + A.neleb();
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  if ((neleA % 2) == 1) {
    // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());
  }
  else {
    // sort: (A',A,B',B) --> (A,B,A',B')
    SMITH::sort_indices<1,3,0,2,0,1,1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());
  }

  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_abFlip<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::AnnihilateAlpha, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::CreateAlpha});

  shared_ptr<const Matrix> Kmatrix = jop_->template coulomb_matrix<0,1,1,0>();

  Matrix tmp = *gamma_A * (*Kmatrix) ^ *gamma_B;

  // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());

  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_abET<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateAlpha});
  auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateAlpha});

  shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = *gamma_A * (*Jmatrix) ^ *gamma_B;

  // sort: (A',A,B',B) --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<1,3,0,2,0,1,-1,1>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());

  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_aaET<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
  auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = *gamma_A * (*Jmatrix) ^ *gamma_B;

  // sort: (A',A,B',B) --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<1,3,0,2,0,1,-1,2>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());

  return out;
}


template <>
shared_ptr<Matrix> MEH_base::compute_bbET<true>(const array<MonomerKey,4>& keys) {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Jmatrix = jop_->template coulomb_matrix<0,0,1,1>();

  Matrix tmp = *gamma_A * (*Jmatrix) ^ *gamma_B;

  // sort: (A',A,B',B) --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  SMITH::sort_indices<1,3,0,2,0,1,-1,2>(tmp.data(), out->data(), Ap.nstates(), A.nstates(), Bp.nstates(), B.nstates());

  return out;
}


