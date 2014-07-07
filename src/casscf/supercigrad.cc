//
// BAGEL - Parallel electron correlation program.
// Filename: supercigrad.cc
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
#include <src/grad/cpcasscf.h>
#include <src/casscf/qvec.h>
#include <src/prop/multipole.h>

using namespace std;
using namespace bagel;

template<typename T>
static string tostring(const T i) {
  stringstream ss;
  ss << i;
  return ss.str();
};


template<>
std::shared_ptr<GradFile> GradEval<SuperCIGrad>::compute() {

  shared_ptr<const Coeff> coeff = ref_->coeff();
  assert(task_->coeff() == coeff);

  const int target = task_->target_state();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();

  const MatView ocoeff = ref_->coeff()->slice(0,nocc);

  // state-averaged density matrices
  shared_ptr<const RDM<1>> rdm1_av = task_->fci()->rdm1_av();
  shared_ptr<const RDM<2>> rdm2_av = task_->fci()->rdm2_av();

  // related to denominators
  const int nmobasis = coeff->mdim();
  assert(nmobasis == nclosed+nact+ref_->nvirt());

  // TODO they are redundant, though...
  shared_ptr<DFHalfDist> half  = geom_->df()->compute_half_transform(ocoeff)->apply_J();
  shared_ptr<DFHalfDist> halfjj = half->apply_J();

  // orbital derivative is nonzero
  auto g0 = make_shared<Matrix>(nmobasis, nmobasis);
  // 1/2 Y_ri = hd_ri + K^{kl}_{rj} D^{lk}_{ji}
  //          = hd_ri + (kr|G)(G|jl) D(lj, ki)
  // 1) one-electron contribution
  auto hmo = make_shared<const Matrix>(*ref_->coeff() % *ref_->hcore() * ocoeff);
  shared_ptr<const Matrix> rdm1 = ref_->rdm1_mat(target);
  assert(rdm1->ndim() == nocc && rdm1->mdim() == nocc);
  g0->add_block(2.0, 0, 0, nmobasis, nocc, *hmo * *rdm1);
  // 2) two-electron contribution
  shared_ptr<const DFFullDist> full  = half->compute_second_transform(ocoeff);
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(*ref_->rdm2(target), *ref_->rdm1(target), nclosed, nact);
  shared_ptr<const Matrix> buf = half->form_2index(fulld, 1.0);
  g0->add_block(2.0, 0, 0, nmobasis, nocc, *ref_->coeff() % *buf);

  // Recalculate the CI vectors (which can be avoided... TODO)
  shared_ptr<const Dvec> civ = task_->fci()->civectors();

  // CI derivative is zero
  auto g1 = make_shared<Dvec>(task_->fci()->det(), ref_->nstate());
  // combine gradient file
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  // compute unrelaxed dipole...
  shared_ptr<Matrix> dtot = ref_->rdm1_mat(target)->resize(nmobasis, nmobasis);
  {
    Dipole dipole(geom_, make_shared<Matrix>(*ref_->coeff() * *dtot ^ *ref_->coeff()), "Unrelaxed");
    dipole.compute();
  }

  // solve CP-CASSCF
  auto cp = make_shared<CPCASSCF>(grad, civ, half, halfjj, ref_, task_->fci());
  shared_ptr<const Matrix> zmat, xmat, dummy;
  shared_ptr<const Dvec> zvec;
  tie(zmat, zvec, xmat, dummy) = cp->solve(task_->thresh());

  // form Zd + dZ^+
  shared_ptr<const Matrix> dsa = rdm1_av->rdm1_mat(nclosed)->resize(nmobasis, nmobasis);
  auto dm = make_shared<Matrix>(*zmat * *dsa + (*dsa ^ *zmat));

  dtot->ax_plus_y(1.0, dm);

  // form zdensity
  auto detex = make_shared<Determinants>(task_->fci()->norb(), task_->fci()->nelea(), task_->fci()->neleb(), false, /*mute=*/true);
  shared_ptr<const RDM<1>> zrdm1;
  shared_ptr<const RDM<2>> zrdm2;
  tie(zrdm1, zrdm2) = task_->fci()->compute_rdm12_av_from_dvec(zvec, civ, detex);

  shared_ptr<Matrix> zrdm1_mat = zrdm1->rdm1_mat(nclosed, false)->resize(nmobasis, nmobasis);
  zrdm1_mat->symmetrize();
  dtot->ax_plus_y(1.0, zrdm1_mat);

  // here dtot is the relaxed 1RDM in the MO basis
  auto dtotao = make_shared<Matrix>(*ref_->coeff() * *dtot ^ *ref_->coeff());

  // compute relaxed dipole moment
  {
    Dipole dipole(geom_, dtotao, "Relaxed");
    dipole.compute();
  }

  // xmat in the AO basis
  auto xmatao = make_shared<Matrix>(*ref_->coeff() * *xmat ^ *ref_->coeff());

  //- TWO ELECTRON PART -//
  // half is computed long before
  shared_ptr<const DFFullDist> qij  = halfjj->compute_second_transform(ocoeff);
  shared_ptr<DFHalfDist> qri;
  {
    shared_ptr<const Matrix> ztrans = make_shared<Matrix>(*ref_->coeff() * zmat->slice(0,nocc));
    {
      const RDM<2> D(*ref_->rdm2(target)+*zrdm2);
      const RDM<1> dd(*ref_->rdm1(target)+*zrdm1);

      shared_ptr<DFFullDist> qijd = qij->apply_2rdm(D, dd, nclosed, nact);
      qijd->ax_plus_y(2.0, halfjj->compute_second_transform(*ztrans)->apply_2rdm(*rdm2_av, *rdm1_av, nclosed, nact));
      qri = qijd->back_transform(ocoeff);
    }
    {
      shared_ptr<const DFFullDist> qijd2 = qij->apply_2rdm(*rdm2_av, *rdm1_av, nclosed, nact);
      qri->ax_plus_y(2.0, qijd2->back_transform(*ztrans));
    }
  }

  shared_ptr<const Matrix> qq  = qri->form_aux_2index(halfjj, 1.0);
  shared_ptr<const DFDist> qrs = qri->back_transform(ocoeff);

  shared_ptr<GradFile> gradient = contract_gradient(dtotao, xmatao, qrs, qq);
  gradient->print();

  return gradient;
}
