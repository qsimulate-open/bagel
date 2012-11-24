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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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
#include <chrono>

using namespace std;
using namespace std::chrono;
using namespace bagel;

template<>
shared_ptr<GradFile> GradEval<SCF<1> >::compute() {
  auto tp0 = high_resolution_clock::now();

  //- One ELECTRON PART -//
#if 1
  mpi__->broadcast(ref_->coeff()->data(), ref_->coeff()->size(), 0);
#endif

  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1(new Matrix(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
  shared_ptr<const Matrix> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_closed_2RDM();
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
qq->print();
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  auto tp1 = high_resolution_clock::now();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right <<
          setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<UHF>::compute() {
  auto tp0 = high_resolution_clock::now();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1(new Matrix(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ));
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  auto tp1 = high_resolution_clock::now();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right <<
          setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<ROHF>::compute() {
  auto tp0 = high_resolution_clock::now();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1(new Matrix(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ));
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();
  assert(erdm1 != nullptr);

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  auto tp1 = high_resolution_clock::now();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right <<
          setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<WernerKnowles>::compute() {
  auto tp0 = high_resolution_clock::now();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1(new Matrix(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  auto tp1 = high_resolution_clock::now();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right <<
          setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl << endl;

  return grad;
}

template<>
shared_ptr<GradFile> GradEval<SuperCI>::compute() {
  auto tp0 = high_resolution_clock::now();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix> rdm1(new Matrix(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
#if 0
  Dipole d(ref_->geom(), rdm1);
  d.compute();
#endif
  shared_ptr<const Matrix> erdm1 = ref_->erdm1();

  //- TWO ELECTRON PART -//
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<const DFFullDist> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_JJ();
  shared_ptr<const DFFullDist> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  shared_ptr<const Matrix> qq  = qij->form_aux_2index(qijd, 1.0);
  shared_ptr<const DFDist> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  auto tp1 = high_resolution_clock::now();
  cout << setw(50) << left << "  * Gradient computed with " << setprecision(2) << right <<
          setw(10) << duration_cast<milliseconds>(tp1-tp0).count()*0.001 << endl << endl;

  return grad;
}
