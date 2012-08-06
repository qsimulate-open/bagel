//
// Newint - Parallel electron correlation program.
// Filename: gradeval.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the Newint package (to be renamed).
//
// The Newint package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The Newint package is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the Newint package; see COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//


#include <src/grad/gradeval.h>

using namespace std;

template<>
shared_ptr<GradFile> GradEval<SCF<1> >::compute() {
  const size_t start = ::clock();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
  shared_ptr<const Matrix1e> erdm1 = ref_->coeff()->form_weighted_density_rhf(ref_->nocc(), ref_->eig());

  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<DF_Full> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_closed_2RDM();
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<UHF>::compute() {
  const size_t start = ::clock();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ));
  shared_ptr<const Matrix1e> erdm1 = ref_->erdm1();
  assert(static_cast<bool>(erdm1));

  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<DF_Full> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<ROHF>::compute() {
  const size_t start = ::clock();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref_->rdm1_mat(0) ^ *coeff_occ));
  shared_ptr<const Matrix1e> erdm1 = ref_->erdm1();
  assert(static_cast<bool>(erdm1));

  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<DF_Full> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_uhf_2RDM(ref_->rdm1(1)->data(), ref_->rdm1(2)->data()); // 1 and 2: alpha and beta
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  return grad;
}


template<>
shared_ptr<GradFile> GradEval<WernerKnowles>::compute() {
  const size_t start = ::clock();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
  shared_ptr<const Matrix1e> erdm1 = ref_->erdm1(); 

  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<DF_Full> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  return grad;
}

template<>
shared_ptr<GradFile> GradEval<SuperCI>::compute() {
  const size_t start = ::clock();

  //- One ELECTRON PART -//
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref_->rdm1_mat() ^ *coeff_occ));
#if 0
  Dipole d(ref_->geom(), rdm1);
  d.compute();
#endif
  shared_ptr<const Matrix1e> erdm1 = ref_->erdm1(); 

  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(coeff_occ->data(), ref_->nocc());
  shared_ptr<DF_Full> qij  = half->compute_second_transform(coeff_occ->data(), ref_->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_2rdm(ref_->rdm2(0)->data(), ref_->rdm1(0)->data(), ref_->nclosed(), ref_->nact());
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref_->coeff()->data())->back_transform(ref_->coeff()->data());

  shared_ptr<GradFile> grad = contract_gradient(rdm1, erdm1, qrs, qq);
  grad->print();

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

  return grad;
}
