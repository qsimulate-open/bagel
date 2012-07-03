//
// Newint - Parallel electron correlation program.
// Filename: supercigrad.cc
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
#include <src/grad/cpcasscf.h>
#include <src/casscf/supercigrad.h>
#include <src/util/pairfile.h>

using namespace std;

template<typename T>
static string tostring(const T i) {
  stringstream ss;
  ss << i; 
  return ss.str();
};


template<>
std::shared_ptr<GradFile> GradEval<SuperCIGrad>::compute() {

  // TODO
  const int target = 0;

  // related to denominators
  const int nbasis = ref_->geom()->nbasis();
  // TODO TODO TODO temp... to have it compiled
  vector<double> eig(nbasis*nbasis);

  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  // TODO they are redundant, though...
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(ref_->coeff()->data(), nocc)->apply_J();
  shared_ptr<DF_Half> halfjj = half->apply_J();

  shared_ptr<FCI> fci(new FCI(multimap<string, string>(), ref_, nclosed, nact));

  int la, lb; tie(la, lb) = fci->len_string();
  // orbital derivative is nonzero
  shared_ptr<Matrix1e> g0(new Matrix1e(ref_->geom())); 
  // 1/2 Y_ri = hd_ri + 2 K^{kl}_{rj} D^{lk}_{ji}
  //          = hd_ri + 2 (kr|G)(G|jl) D(lj, ki)
  // 1) one-electron contribution 
  shared_ptr<Matrix1e> hmo(new Matrix1e(*ref_->coeff() % *ref_->hcore() * *ref_->coeff()));
  shared_ptr<Matrix1e> rdm1 = ref_->rdm1_mat(target);
  dgemm_("N", "N", nbasis, nocc, nocc, 2.0, hmo->data(), nbasis, rdm1->data(), nbasis, 0.0, g0->data(), nbasis);
  // 2) two-electron contribution
  shared_ptr<DF_Full> full  = half->compute_second_transform(ref_->coeff()->data(), nocc);
  shared_ptr<DF_Full> fulld = full->apply_2rdm(ref_->rdm2(target)->data(), ref_->rdm1(target)->data(), nclosed, nact);
  unique_ptr<double[]> buf = half->form_2index(fulld);
  dgemm_("T", "N", nbasis, nocc, nbasis, 2.0, ref_->coeff()->data(), nbasis, buf.get(), nbasis, 1.0, g0->data(), nbasis); 

  // CI derivative is zero
  shared_ptr<Dvec> g1(new Dvec(lb, la, ref_->nstate()));
  // combine gradient file
  shared_ptr<PairFile<Matrix1e, Dvec> > grad(new PairFile<Matrix1e, Dvec>(g0, g1));

  // solve CP-CASSCF
  shared_ptr<CPCASSCF> cp(new CPCASSCF(grad, eig, half, halfjj, ref_, fci));
  shared_ptr<PairFile<Matrix1e, Dvec> > zvec = cp->solve();

  // compute dipole...
  shared_ptr<Matrix1e> dtot = ref_->rdm1_mat(target);

  // computes dipole mements
  shared_ptr<const Matrix1e> coeff_occ = ref_->coeff()->slice(0,ref_->nocc());
  shared_ptr<Matrix1e> dtotao(new Matrix1e(*coeff_occ * *dtot ^ *coeff_occ));
  Dipole dipole(geom_, dtotao);
  dipole.compute();

  std::shared_ptr<GradFile> out(new GradFile(3*geom_->natom()));
  return out;
}

