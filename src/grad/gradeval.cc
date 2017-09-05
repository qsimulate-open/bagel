//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: gradeval.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/grad/gradeval.h>
#include <src/util/timer.h>

using namespace std;
using namespace bagel;

template<>
shared_ptr<GradFile> GradEval<RHF>::compute(const std::string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  assert(task_->dodf());
  Timer timer;

  //- One ELECTRON PART -//
  const MatView coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(coeff_occ * *ref_->rdm1_mat() ^ coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  const bool external_half = static_cast<bool>(task_->half());
  shared_ptr<const DFHalfDist> half;
  if (external_half) {
    half = task_->half();
    task_->discard_half();
  } else {
    half = geom_->df()->compute_half_transform(coeff_occ);
  }

  shared_ptr<const DFFullDist> qij = external_half ? half->compute_second_transform(coeff_occ)->apply_J()
                                                   : half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM();
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFHalfDist> qrs_1 = qijd->back_transform(coeff_occ);
  shared_ptr<const DFDist> qrs = qrs_1->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  dipole_ = task_->scf_dipole();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  energy_ = ref_->energy(0);

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<UHF>::compute(const std::string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  Timer timer;

  //- One ELECTRON PART -//
  const MatView coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(coeff_occ * *ref_->rdm1_mat(0) ^ coeff_occ);
  shared_ptr<const Matrix> erdm1 = task_->compute_erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(*ref_->rdm1(1), *ref_->rdm1(2)); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  dipole_ = task_->scf_dipole();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  energy_ = ref_->energy(0);

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<ROHF>::compute(const std::string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  Timer timer;

  //- One ELECTRON PART -//
  const MatView coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(coeff_occ * *ref_->rdm1_mat(0) ^ coeff_occ);
  shared_ptr<const Matrix> erdm1 = task_->compute_erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(*ref_->rdm1(1), *ref_->rdm1(2)); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  dipole_ = task_->scf_dipole();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  energy_ = ref_->energy(0);

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<KS>::compute(const std::string jobtitle, shared_ptr<const GradInfo> gradinfo) {
  Timer timer;

  //- One ELECTRON PART -//
  const MatView coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(coeff_occ * *ref_->rdm1_mat() ^ coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  // ... exchange needs to be scaled.
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM(task_->func()->scale_ex());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  //- Exchange-correlation part -//
  shared_ptr<const GradFile> ggrad = task_->grid()->compute_xcgrad(task_->func(), make_shared<Matrix>(coeff_occ));
  *grad += *ggrad;

  grad->print();

  dipole_ = task_->scf_dipole();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  energy_ = ref_->energy(0);

  return grad;
}
