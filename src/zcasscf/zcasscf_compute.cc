//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_compute.cc
// Copyright (C) 2013 Toru Shiozaki
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

#include <src/zcasscf/zqvec.h>
#include <src/math/bfgs.h>
#include <src/rel/dfock.h>
#include <src/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;


void ZCASSCF::compute() {
  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.
  shared_ptr<BFGS<ZMatrix>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  auto x = make_shared<ZMatrix>(nbasis_, nbasis_);
  x->unit();
  shared_ptr<const ZMatrix> xstart;

  shared_ptr<const ZMatrix> hcore = make_shared<RelHcore>(geom_);

  for (int iter = 0; iter != max_iter_; ++iter) {
    // first perform CASCI to obtain RDMs
    mute_stdcout();
    if (iter) fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
    // get energy
    energy_ = fci_->energy();
    resume_stdcout();

    // closed Fock operator
    shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore, coeff_->slice(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore;
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);

    // active Fock operator
    shared_ptr<const ZMatrix> afockao = active_fock();
    shared_ptr<const ZMatrix> afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    assert(coeff_->mdim()== nbasis_*2);

    // qvec
    shared_ptr<const ZMatrix> qvec = make_shared<ZQvec>(nbasis_, nact_, geom_, coeff_, nclosed_, fci_, gaunt_, breit_);

    // compute orbital gradients

    // print energy
    const double gradient = 0.0; // TODO
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());
  }
}


shared_ptr<const ZMatrix> ZCASSCF::active_fock() const {
  assert(fci_);

  array<shared_ptr<const ZMatrix>,2> kcoeff = fci_->kramers_coeff(); 

  auto rdm1_tot = make_shared<ZMatrix>(nact_*2, nact_*2);
  rdm1_tot->copy_block(    0,     0, nact_, nact_, fci_->rdm1_av("00")->data());
  rdm1_tot->copy_block(nact_, nact_, nact_, nact_, fci_->rdm1_av("11")->data());
  rdm1_tot->copy_block(nact_,     0, nact_, nact_, fci_->rdm1_av("10")->data());
  rdm1_tot->copy_block(    0, nact_, nact_, nact_, rdm1_tot->get_submatrix(nact_, 0, nact_, nact_)->transpose_conjg());

  auto coeff_tot = make_shared<ZMatrix>(kcoeff[0]->ndim(), nact_*2);
  assert(nact_ == kcoeff[0]->mdim() && nact_ == kcoeff[1]->mdim() && kcoeff[0]->ndim() % 4 == 0); 
  coeff_tot->copy_block(0,     0, kcoeff[0]->ndim(), nact_, kcoeff[0]);
  coeff_tot->copy_block(0, nact_, kcoeff[1]->ndim(), nact_, kcoeff[1]);

  // form natural orbitals
  {
    unique_ptr<double[]> eig(new double[nact_*2]);
    rdm1_tot->diagonalize(eig.get());
    *coeff_tot *= *rdm1_tot;

    // scale using eigen values
    for (int i = 0; i != nact_*2; ++i) {
      assert(eig[i] >= -1.0e-14);
      const double fac = eig[i] > 0 ? sqrt(eig[i]) : 0.0;
      transform(coeff_tot->element_ptr(0, i), coeff_tot->element_ptr(0, i+1), coeff_tot->element_ptr(0, i), [&fac](complex<double> a) { return fac*a; }); 
    }
  }

  auto zero = make_shared<ZMatrix>(geom_->nbasis()*4, geom_->nbasis()*4);
  return make_shared<const DFock>(geom_, zero, coeff_tot, gaunt_, breit_, /*store half*/false, /*robust*/breit_);
}



void ZCASSCF::grad_vc(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<ZRotFile> sigma) const {

}


void ZCASSCF::grad_va(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> qxr,   shared_ptr<ZRotFile> sigma) const {

}


void ZCASSCF::grad_ca(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<ZRotFile> sigma) const {

}
