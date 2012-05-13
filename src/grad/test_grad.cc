//
// Newint - Parallel electron correlation program.
// Filename: test_grad.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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

#include <src/df/fit.h>
#include <src/wfn/reference.h>
#include <iostream>
#include <src/osint/kineticbatch.h>

#include <src/grad/gradeval_base.h>

using namespace std;

void test_grad(shared_ptr<Reference> ref) {

  const size_t start = ::clock();
  cout << "  testing grad.." << endl;

  // target quantity here ... ==========

  GradEval_base g1(ref->geom());

  shared_ptr<const Matrix1e> coeff_occ = ref->coeff()->slice(0,ref->nocc());
  shared_ptr<const Matrix1e> rdm1(new Matrix1e(*coeff_occ * *ref->rdm1() ^ *coeff_occ));
  shared_ptr<const Matrix1e> erdm1 = ref->coeff()->form_weighted_density_rhf(ref->nocc(), ref->eig());

  vector<double> grad = g1.contract_grad1e(rdm1, erdm1);


  //- TWO ELECTRON PART -//
  shared_ptr<DF_Half> half = ref->geom()->df()->compute_half_transform(coeff_occ->data(), ref->nocc());
  shared_ptr<DF_Full> qij  =  half->compute_second_transform(coeff_occ->data(), ref->nocc())->apply_J()->apply_J();
  shared_ptr<DF_Full> qijd = qij->apply_closed_2RDM();
  unique_ptr<double[]> qq  = qij->form_aux_2index(qijd);
  shared_ptr<DF_AO> qrs = qijd->back_transform(ref->coeff()->data())->back_transform(ref->coeff()->data());

  vector<double> gradtmp = g1.contract_grad2e(qrs); 
  for (auto i = grad.begin(), j = gradtmp.begin(); i != grad.end(); ++i, ++j) *i += *j;

  vector<double> gradtmp2 = g1.contract_grad2e_2index(qq);
  for (auto i = grad.begin(), j = gradtmp2.begin(); i != grad.end(); ++i, ++j) *i += *j;


  cout << endl << "  * Nuclear energy gradient" << endl << endl;
  for (int i = 0; i != ref->geom()->natom(); ++i) {
    cout << "    o Atom " << setw(3) << i << endl;
    cout << "        x  " << setprecision(10) << setw(20) << fixed << grad[3*i+0] << endl;
    cout << "        y  " << setprecision(10) << setw(20) << fixed << grad[3*i+1] << endl;
    cout << "        z  " << setprecision(10) << setw(20) << fixed << grad[3*i+2] << endl;
  }

  cout << setw(50) << left << "  * Gradient computed with " << setprecision(3) << right <<
          setw(10) << (::clock() - start)/static_cast<double>(CLOCKS_PER_SEC) << endl << endl;

}


