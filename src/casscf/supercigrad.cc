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
#include <src/casscf/qvec.h>
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

  shared_ptr<const Coeff> coeff = ref_->coeff();

  const int target = 0;
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nstate = ref_->nstate();
  const int nvirt = ref_->nvirt();

  // related to denominators
  const int nbasis = ref_->geom()->nbasis();
  shared_ptr<Matrix1e> eig(new Matrix1e(ref_->geom()));
  {
    // as in Theor Chem Acc (1997) 97:88-95
    vector<double> occup_ = task_->fci()->rdm1(target)->diag();

    shared_ptr<Matrix1e> deninact = task_->ao_rdm1(task_->fci()->rdm1(target), true); // true means inactive_only
    shared_ptr<Matrix1e> f_inactao(new Matrix1e(geom_));
    dcopy_(nbasis*nbasis, task_->fci()->jop()->core_fock_ptr(), 1, f_inactao->data(), 1);
    shared_ptr<Matrix1e> finact (new Matrix1e(*coeff % *f_inactao * *coeff));

    shared_ptr<Matrix1e> denall = task_->ao_rdm1(task_->fci()->rdm1(target));
    shared_ptr<Matrix1e> denact (new Matrix1e(*denall-*deninact));
    shared_ptr<Fock<1> > fact_ao(new Fock<1>(geom_, task_->hcore(), denact, ref_->schwarz()));
    shared_ptr<Matrix1e> f      (new Matrix1e(*finact+ *coeff%(*fact_ao-*task_->hcore())**coeff));

    shared_ptr<Qvec> fact(new Qvec(nbasis, nact, ref_->geom()->df(), ref_->coeff(), nclosed, task_->fci(), task_->fci()->rdm2(target)));
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
#if 0
    for (int i = 0; i != nact; ++i)
      for (int j = 0; j != nact; ++j)
        eig->element(j+nclosed,i+nclosed) = eig->element(i+nclosed,j+nclosed) = 1.0e2; 
#endif

  }

  // TODO they are redundant, though...
  shared_ptr<DF_Half> half = ref_->geom()->df()->compute_half_transform(ref_->coeff()->data(), nocc)->apply_J();
  shared_ptr<DF_Half> halfjj = half->apply_J();

  // orbital derivative is nonzero
  shared_ptr<Matrix1e> g0(new Matrix1e(ref_->geom())); 
  // 1/2 Y_ri = hd_ri + K^{kl}_{rj} D^{lk}_{ji}
  //          = hd_ri + (kr|G)(G|jl) D(lj, ki)
  // 1) one-electron contribution 
  shared_ptr<const Matrix1e> hmo(new Matrix1e(*ref_->coeff() % *ref_->hcore() * *ref_->coeff()));
  shared_ptr<const Matrix1e> rdm1 = ref_->rdm1_mat(target);
  dgemm_("N", "N", nbasis, nocc, nocc, 2.0, hmo->data(), nbasis, rdm1->data(), nbasis, 0.0, g0->data(), nbasis);
  // 2) two-electron contribution
  shared_ptr<const DF_Full> full  = half->compute_second_transform(ref_->coeff()->data(), nocc);
  shared_ptr<const DF_Full> fulld = full->apply_2rdm(ref_->rdm2(target)->data(), ref_->rdm1(target)->data(), nclosed, nact);
  unique_ptr<double[]> buf = half->form_2index(fulld);
  dgemm_("T", "N", nbasis, nocc, nbasis, 2.0, ref_->coeff()->data(), nbasis, buf.get(), nbasis, 1.0, g0->data(), nbasis); 

  // Recalculate the CI vectors (which can be avoided... TODO)
  shared_ptr<const Dvec> civ = task_->fci()->civectors();

  // CI derivative is zero
  shared_ptr<Dvec> g1(new Dvec(task_->fci()->det(), ref_->nstate()));
  // combine gradient file
  shared_ptr<PairFile<Matrix1e, Dvec> > grad(new PairFile<Matrix1e, Dvec>(g0, g1));

  // solve CP-CASSCF
  shared_ptr<CPCASSCF> cp(new CPCASSCF(grad, civ, eig, half, halfjj, ref_, task_->fci()));
  shared_ptr<PairFile<Matrix1e, Dvec> > zvec = cp->solve();

  // form Zd + dZ^+
  shared_ptr<Matrix1e> dsa = ref_->rdm1_mat()->expand();
  shared_ptr<Matrix1e> zslice = zvec->first();
  shared_ptr<Matrix1e> dm(new Matrix1e(*zslice * *dsa + (*dsa ^ *zslice)));

  // compute dipole...
  shared_ptr<Matrix1e> dtot = ref_->rdm1_mat(target)->expand();

  dtot->daxpy(1.0, dm);

  // form zdensity 
  shared_ptr<Determinants> detex(new Determinants(task_->fci()->norb(), task_->fci()->nelea(), task_->fci()->neleb(), false));
  shared_ptr<const RDM<1> > zrdm1;
  shared_ptr<const RDM<2> > zrdm2;
  tie(zrdm1, zrdm2) = task_->fci()->compute_rdm12_av_from_dvec(civ, zvec->second(), detex);

  shared_ptr<Matrix1e> zrdm1_mat = zrdm1->rdm1_mat(ref_->geom(), nclosed, false)->expand(); 
  zrdm1_mat->symmetrize();
  dtot->daxpy(1.0, zrdm1_mat);

  // computes dipole mements
  shared_ptr<Matrix1e> dtotao(new Matrix1e(*ref_->coeff() * *dtot ^ *ref_->coeff()));
  Dipole dipole(geom_, dtotao);
  dipole.compute();

  std::shared_ptr<GradFile> out(new GradFile(3*geom_->natom()));
  return out;
}

