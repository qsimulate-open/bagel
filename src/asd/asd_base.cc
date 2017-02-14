//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: asd_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Shane Parker <shane.parker@u.northwestern.edu>
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

#include <src/asd/asd_base.h>
#include <src/util/prim_op.h>

using namespace std;
using namespace bagel;

ASD_base::ASD_base(const shared_ptr<const PTree> input, shared_ptr<const Dimer> dimer, bool rdm) : dimer_(dimer), compute_rdm_(rdm) {
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
  print_info_ = input->get<bool>("print_info", false);

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
        this_model.emplace_back(make_pair(spins[0],spins[1]), make_pair(charges[0],charges[1]), make_pair(p,p), nstates);
      }
      models_to_form_.emplace_back(move(this_model));
    }
  }

  Timer timer;

  shared_ptr<const Reference> dimerref = dimer_->sref();

  jop_ = make_shared<DimerJop>(dimerref, dimerref->nclosed(), dimerref->nclosed() + dimer_->active_refs().first->nact(), dimerref->nclosed() + dimerref->nact(), dimerref->coeff());
  cout << "  o computing integrals: " << timer.tick() << endl;

  energies_ = vector<double>(nstates_, 0.0);

  for (int i = 0; i != nstates_; ++i) weight_.push_back(1.0/static_cast<double>(nstates_));

  // resizing rdm vectors (with null pointers)
  rdm1_.resize(nstates_);
  rdm2_.resize(nstates_);
}


void ASD_base::update_dimer_and_fix_ci(shared_ptr<const Dimer> dimer) {
  Timer timer;
  //update reference & integrals
  dimer_ = make_shared<const Dimer>(*dimer);
  shared_ptr<const Reference> dimerref = dimer_->sref();
  jop_ = make_shared<DimerJop>(dimerref, dimerref->nclosed(), dimerref->nclosed() + dimer_->active_refs().first->nact(), dimerref->nclosed() + dimerref->nact(), dimerref->coeff());
  cout << "  o computing integrals: " << timer.tick() << endl;
  // initialize
  energies_ = vector<double>(nstates_, 0.0);
  rdm1_.clear();
  rdm2_.clear();
  // resizing rdm vectors (with null pointers)
  rdm1_.resize(nstates_);
  rdm2_.resize(nstates_);
  //fix ci coefficients for subsequent asd calculations
  fix_ci_ = true;
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


shared_ptr<Matrix> ASD_base::compute_diagonal_block(const DimerSubspace_base& subspace) const {
  const double core = dimer_->sref()->geom()->nuclear_repulsion() + jop_->core_energy();

  auto out = compute_intra(subspace, jop_, core);
  array<MonomerKey,4> keys {{ subspace.monomerkey<0>(), subspace.monomerkey<1>(), subspace.monomerkey<0>(), subspace.monomerkey<1>() }};
  *out += *compute_inter_2e(keys);

  return out;
}


shared_ptr<Matrix> ASD_base::compute_offdiagonal_1e(const array<MonomerKey,4>& keys, shared_ptr<const Matrix> hAB) const {
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
    sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }

  return out;
}


// This term will couple off-diagonal blocks since it has no delta functions involved
shared_ptr<Matrix> ASD_base::compute_inter_2e(const array<MonomerKey,4>& keys) const {
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
  sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  return out;
}


shared_ptr<Matrix> ASD_base::compute_aET(const array<MonomerKey,4>& keys) const {
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
    sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  } else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  return out;
}


shared_ptr<Matrix> ASD_base::compute_bET(const array<MonomerKey,4>& keys) const {
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
    sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }
  else {
    // sort: (A,A',B,B') --> (A,B,A',B')
    sort_indices<0,2,1,3,0,1,1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());
  }

  return out;
}


shared_ptr<Matrix> ASD_base::compute_abFlip(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::AnnihilateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::CreateAlpha, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Kmatrix = jop_->coulomb_matrix<0,1,1,0>();

  Matrix tmp = gamma_A * (*Kmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


shared_ptr<Matrix> ASD_base::compute_abET(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];

  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -1.0 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  sort_indices<0,2,1,3,0,1,-1,1>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


shared_ptr<Matrix> ASD_base::compute_aaET(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateAlpha, GammaSQ::CreateAlpha});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateAlpha, GammaSQ::AnnihilateAlpha});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  sort_indices<0,2,1,3,0,1,-1,2>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
}


shared_ptr<Matrix> ASD_base::compute_bbET(const array<MonomerKey,4>& keys) const {
  auto& A = keys[0]; auto& B = keys[1]; auto& Ap = keys[2]; auto& Bp = keys[3];
  auto gamma_A = gammatensor_[0]->get_block_as_matview(A, Ap, {GammaSQ::CreateBeta, GammaSQ::CreateBeta});
  auto gamma_B = gammatensor_[1]->get_block_as_matview(B, Bp, {GammaSQ::AnnihilateBeta, GammaSQ::AnnihilateBeta});

  shared_ptr<const Matrix> Jmatrix = jop_->coulomb_matrix<0,0,1,1>();

  Matrix tmp = gamma_A * (*Jmatrix) ^ gamma_B;

  // sort: (A,A',B,B') --> -0.5 * (A,B,A',B')
  auto out = make_shared<Matrix>(A.nstates()*B.nstates(), Ap.nstates()*Bp.nstates());
  sort_indices<0,2,1,3,0,1,-1,2>(tmp.data(), out->data(), A.nstates(), Ap.nstates(), B.nstates(), Bp.nstates());

  return out;
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
         << "   <S^2> = " << setw(4) << setprecision(4) << fixed << blas::dot_product(spn->element_ptr(0,istate), dimerstates_, cc.element_ptr(0,istate)) << endl;
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
  const int nprint = min(nstates, static_cast<int>(property->ndim()));

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


shared_ptr<Matrix> ASD_base::couple_blocks(const DimerSubspace_base& AB, const DimerSubspace_base& ApBp) const {

  Coupling term_type = coupling_type(AB, ApBp);

  const DimerSubspace_base* space1 = &AB;
  const DimerSubspace_base* space2 = &ApBp;

  const bool flip = (static_cast<int>(term_type) < 0);
  if (flip) {
    term_type = Coupling(-1*static_cast<int>(term_type));
    std::swap(space1,space2);
  }

  shared_ptr<Matrix> out;
  array<MonomerKey,4> keys {{space1->template monomerkey<0>(), space1->template monomerkey<1>(), space2->template monomerkey<0>(), space2->template monomerkey<1>()}};

  switch(term_type) {
    case Coupling::none :
      out = nullptr; break;
    case Coupling::diagonal :
      out = compute_inter_2e(keys); break;
    case Coupling::aET :
      out = compute_aET(keys); break;
    case Coupling::bET :
      out = compute_bET(keys); break;
    case Coupling::abFlip :
      out = compute_abFlip(keys); break;
    case Coupling::abET :
      out = compute_abET(keys); break;
    case Coupling::aaET :
      out = compute_aaET(keys); break;
    case Coupling::bbET :
      out = compute_bbET(keys); break;
    default :
      throw logic_error("Asking for a coupling type that has not been written.");
  }

  /* For the Hamiltonian with flip = true, we tranpose the output */
  if (flip) out = out->transpose();

  return out;
}
