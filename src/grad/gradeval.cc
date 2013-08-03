//
// BAGEL - Parallel electron correlation program.
// Filename: gradeval.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/grad/gradeval.h>
#include <src/util/timer.h>

using namespace std;
using namespace bagel;

template<>
shared_ptr<GradFile> GradEval<SCF>::compute() {
  assert(task_->dodf());
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM();
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<UHF>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<ROHF>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<KS>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
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
  shared_ptr<const GradFile> ggrad = task_->grid()->compute_xcgrad(task_->func(), coeff_occ);
  *grad += *ggrad;

  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<WernerKnowles>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}

template<>
shared_ptr<GradFile> GradEval<SuperCI>::compute() {
  Timer timer;

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1 = make_shared<Matrix>(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ);
#if 0
  Dipole d(geom_, rdm1);
  d.compute();
#endif
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(coeff_occ);
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ)->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(coeff_occ)->back_transform(coeff_occ);

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right << setw(10) << timer.tick() << endl << endl;

  return grad;
}
