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

  const int target = 0;
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nstate = ref_->nstate();
  const int nvirt = ref_->nvirt();

  // related to denominators
  const int nbasis = geom_->nbasis();
  auto eig = make_shared<Matrix>(nbasis, nbasis);
  {
    // as in Theor Chem Acc (1997) 97:88-95
    vector<double> occup_ = task_->fci()->rdm1(target)->diag();

    shared_ptr<Matrix> deninact = task_->ao_rdm1(task_->fci()->rdm1(target), true); // true means inactive_only
    auto f_inactao = make_shared<Matrix>(nbasis, nbasis);
    copy_n(task_->fci()->jop()->core_fock_ptr(), nbasis*nbasis, f_inactao->data()); // TODO copy construct?
    auto finact = make_shared<Matrix>(*coeff % *f_inactao * *coeff);

    shared_ptr<Matrix> denall = task_->ao_rdm1(task_->fci()->rdm1(target));
    auto denact = make_shared<Matrix>(*denall-*deninact);
    auto fact_ao = make_shared<Fock<1>>(geom_, task_->hcore(), denact, ref_->schwarz());
    auto f = make_shared<Matrix>(*finact+ *coeff%(*fact_ao-*task_->hcore())**coeff);

    auto fact = make_shared<Qvec>(nbasis, nact, geom_->df(), ref_->coeff(), nclosed, task_->fci(), task_->fci()->rdm2(target));
    for (int i = 0; i != nact; ++i)
      daxpy_(nbasis, occup_[i], finact->element_ptr(0,nclosed+i), 1, fact->data()+i*nbasis, 1);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nvirt; ++j)
        eig->element(j+nocc,i+nclosed) = eig->element(i+nclosed,j+nocc) = -fact->element(i,i) + occup_[i]*f->element(j+nocc, j+nocc);

    for (int i = 0; i != nclosed; ++i)
      for (int j = 0; j != nvirt; ++j)
         eig->element(j+nocc,i) = eig->element(i,j+nocc) = 2.0*f->element(j+nocc, j+nocc) - 2.0*f->element(i, i);

    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nclosed; ++j)
         eig->element(j,i+nclosed) = eig->element(i+nclosed,j)
                                   = (f->element(nclosed+i,nclosed+i)*2.0-fact->element(i+nclosed,i)) - f->element(j, j)*(2.0 - occup_[i]);
#if 1
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        eig->element(j+nclosed,i+nclosed) = eig->element(i+nclosed,j+nclosed) = 1.0e0;
#endif

  }

  // TODO they are redundant, though...
  shared_ptr<DFHalfDist> half = geom_->df()->compute_half_transform(ref_->coeff()->slice(0,nocc))->apply_J();
  shared_ptr<DFHalfDist> halfjj = half->apply_J();

  // orbital derivative is nonzero
  auto g0 = make_shared<Matrix>(nbasis, nbasis);
  // 1/2 Y_ri = hd_ri + K^{kl}_{rj} D^{lk}_{ji}
  //          = hd_ri + (kr|G)(G|jl) D(lj, ki)
  // 1) one-electron contribution
  auto hmo = make_shared<const Matrix>(*ref_->coeff() % *ref_->hcore() * *ref_->coeff());
  shared_ptr<const Matrix> rdm1 = ref_->rdm1_mat(target);
  dgemm_("N", "N", nbasis, nocc, nocc, 2.0, hmo->data(), nbasis, rdm1->data(), nbasis, 0.0, g0->data(), nbasis);
  // 2) two-electron contribution
  shared_ptr<const DFFullDist> full  = half->compute_second_transform(ref_->coeff()->slice(0,nocc));
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(ref_->rdm2(target)->data(), ref_->rdm1(target)->data(), nclosed, nact);
  shared_ptr<const Matrix> buf = half->form_2index(fulld, 1.0);
  *g0 += *ref_->coeff() % *buf;

  // Recalculate the CI vectors (which can be avoided... TODO)
  shared_ptr<const Dvec> civ = task_->fci()->civectors();

  // CI derivative is zero
  auto g1 = make_shared<Dvec>(task_->fci()->det(), ref_->nstate());
  // combine gradient file
  auto grad = make_shared<PairFile<Matrix, Dvec>>(g0, g1);

  // solve CP-CASSCF
  auto cp = make_shared<CPCASSCF>(grad, civ, eig, half, halfjj, ref_, task_->fci());
  shared_ptr<PairFile<Matrix, Dvec>> zvec = cp->solve();

  // form Zd + dZ^+
  shared_ptr<Matrix> dsa = ref_->rdm1_mat()->resize(nbasis, nbasis);
  shared_ptr<Matrix> zslice = zvec->first();
  auto dm = make_shared<Matrix>(*zslice * *dsa + (*dsa ^ *zslice));

  // compute dipole...
  shared_ptr<Matrix> dtot = ref_->rdm1_mat(target)->resize(nbasis, nbasis);
  dtot->ax_plus_y(1.0, dm);

  // form zdensity
  auto detex = make_shared<Determinants>(task_->fci()->norb(), task_->fci()->nelea(), task_->fci()->neleb(), false, /*mute=*/true);
  shared_ptr<const RDM<1>> zrdm1;
  shared_ptr<const RDM<2>> zrdm2;
  tie(zrdm1, zrdm2) = task_->fci()->compute_rdm12_av_from_dvec(civ, zvec->second(), detex);

  shared_ptr<Matrix> zrdm1_mat = zrdm1->rdm1_mat(geom_, nclosed, false)->resize(nbasis, nbasis);
  zrdm1_mat->symmetrize();
  dtot->ax_plus_y(1.0, zrdm1_mat);

  // computes dipole mements
  auto dtotao = make_shared<Matrix>(*ref_->coeff() * *dtot ^ *ref_->coeff());
  Dipole dipole(geom_, dtotao);
  dipole.compute();

  return make_shared<GradFile>(geom_->natom());
}

