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
#include <src/rel/reloverlap.h>

using namespace std;
using namespace bagel;


void ZCASSCF::compute() {
  // equation numbers refer to Chaban, Schmidt and Gordon 1997 TCA 97, 88.
  shared_ptr<BFGS<ZRotFile>> bfgs;

  // ============================
  // macro iteration from here
  // ============================
  Timer timer;

  auto x = make_shared<ZMatrix>(nbasis_*2, nbasis_*2);
  x->unit();

  shared_ptr<const ZMatrix> hcore = make_shared<RelHcore>(geom_);
  shared_ptr<const RelOverlap> overlap = make_shared<RelOverlap>(geom_);

  // intialize coefficients
  init_kramers_coeff(hcore, overlap);

  // TODO for debug, we may rotate coefficients. The magnitude can be specified in the input
  const bool ___debug___break_kramers = true;
  if (___debug___break_kramers)
    ___debug___orbital_rotation(/*kramers*/ false);

  if (nact_)
    fci_->update(coeff_);

  cout << " See casscf.log file for further information on FCI output " << endl;
  for (int iter = 0; iter != max_iter_; ++iter) {
    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout();
      if (iter) fci_->update(coeff_);
      cout << " Executing FCI calculation " << endl;
      fci_->compute();
      cout << " Computing RDMs from FCI calculation " << endl;
      fci_->compute_rdm12();
      resume_stdcout();
    }

    // calculate 1RDM in an original basis set
    shared_ptr<const ZMatrix> rdm1 = nact_ ? transform_rdm1() : shared_ptr<const ZMatrix>();

    // closed Fock operator
    shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore, coeff_->slice(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore;
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    // cfock->print("cfock", 7);

    // active Fock operator
    shared_ptr<const ZMatrix> afock;
    if (nact_) {
      shared_ptr<const ZMatrix> afockao = active_fock(rdm1);
      afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    } else {
      afock = make_shared<ZMatrix>(nbasis_*2, nbasis_*2);
    }
    assert(coeff_->mdim()== nbasis_*2);

    // qvec
    shared_ptr<const ZMatrix> qvec;
    if (nact_) {
      qvec = make_shared<ZQvec>(nbasis_, nact_, geom_, coeff_, nclosed_, fci_, gaunt_, breit_);
    }

    // get energy
    if (nact_) {
      energy_ = fci_->energy();
    } else {
      assert(nstate_ == 1);
      energy_.resize(1);
      energy_[0] = geom_->nuclear_repulsion();
      auto mo = make_shared<ZMatrix>(*coeff_ % (*cfockao+*hcore) * *coeff_);
      for (int i = 0; i != nclosed_*2; ++i)
        energy_[0] += 0.5*mo->element(i,i).real();
    }

    if (iter == 0) {
      shared_ptr<const ZRotFile> denom = compute_denom(cfock, afock, qvec, rdm1);
      bfgs = make_shared<BFGS<ZRotFile>>(denom);
    }

    // compute orbital gradients
    shared_ptr<ZRotFile> grad = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2, /*superci*/false);
    grad_vc(cfock, afock, grad);
    grad_va(cfock, qvec, rdm1, grad);
    grad_ca(cfock, afock, qvec, rdm1, grad);
    *grad *= 2.0;
    if (!___debug___break_kramers)
      kramers_adapt(grad);

    if (___debug___break_kramers)
      ___debug___print_gradient(grad);

    auto xlog = make_shared<ZRotFile>(x->log(4), nclosed_*2, nact_*2, nvirt_*2, /*superci*/ false);
    shared_ptr<ZRotFile> a = bfgs->extrapolate(grad, xlog);
    if (!___debug___break_kramers)
      kramers_adapt(a);
    shared_ptr<ZMatrix> amat = a->unpack<ZMatrix>();

    const double gradient = amat->rms();

    // multiply -1 from the formula. multiply -i to make amat hermite (will be compensated)
    *amat *= -1.0 * complex<double>(0.0, -1.0);

    // restore the matrix from RotFile
    unique_ptr<double[]> teig(new double[amat->ndim()]);
    amat->diagonalize(teig.get());
    auto amat_sav = amat->copy();
    for (int i = 0; i != amat->ndim(); ++i) {
      complex<double> ex = exp(complex<double>(0.0, teig[i]));
      for_each(amat->element_ptr(0,i), amat->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
    }
    auto expa = make_shared<ZMatrix>(*amat ^ *amat_sav);
    if (!___debug___break_kramers)
      kramers_adapt(expa);

    coeff_ = make_shared<const ZMatrix>(*coeff_ * *expa);
    // for next BFGS extrapolation
    *x *= *expa;

    // print energy
    print_iteration(iter, 0, 0, energy_, gradient, timer.tick());

    if (gradient < thresh_) break;
  }
}


shared_ptr<const ZMatrix> ZCASSCF::transform_rdm1() const {
  assert(fci_);

  // RDM transform as D_rs = C*_ri D_ij (C*_rj)^+
  auto rdm1_tot = make_shared<ZMatrix>(nact_*2, nact_*2);
  rdm1_tot->copy_block(    0,     0, nact_, nact_, fci_->rdm1_av("00")->data());
  rdm1_tot->copy_block(nact_, nact_, nact_, nact_, fci_->rdm1_av("11")->data());
  rdm1_tot->copy_block(nact_,     0, nact_, nact_, fci_->rdm1_av("10")->data());
  rdm1_tot->copy_block(    0, nact_, nact_, nact_, rdm1_tot->get_submatrix(nact_, 0, nact_, nact_)->transpose_conjg());

  auto coeff_tot = fci_->coeff()->get_conjg();

  // RDM transform as D_ij = (C*_ri)^+ S_rr' D_r's' S_s's C*_sj
  // TODO compute RelOverlap only once (this is comptued also in zqvec)
  auto overlap = make_shared<const RelOverlap>(geom_);
  shared_ptr<const ZMatrix> ocoeff = coeff_->slice(nclosed_*2, nclosed_*2+nact_*2)->get_conjg();
  const ZMatrix co = *ocoeff % *overlap * *coeff_tot;
  return make_shared<ZMatrix>(co * *rdm1_tot ^ co);
}


shared_ptr<const ZMatrix> ZCASSCF::active_fock(shared_ptr<const ZMatrix> rdm1) const {
  // form natural orbitals
  // JEB :  cout << "presently forming natural orbitals via diagonalization" << endl;
  unique_ptr<double[]> eig(new double[nact_*2]);
  auto tmp = make_shared<ZMatrix>(*rdm1);
  tmp->diagonalize(eig.get());
  auto ocoeff = coeff_->slice(nclosed_*2, nclosed_*2+nact_*2);
  // D_rs = C*_ri D_ij (C*_rj)^+. Dij = U_ik L_k (U_jk)^+. So, C'_ri = C_ri * U*_ik
  auto natorb = make_shared<ZMatrix>(*ocoeff * *tmp->get_conjg());

  // scale using eigen values
  for (int i = 0; i != nact_*2; ++i) {
    assert(eig[i] >= -1.0e-14);
    const double fac = eig[i] > 0 ? sqrt(eig[i]) : 0.0;
    for_each(natorb->element_ptr(0, i), natorb->element_ptr(0, i+1), [&fac](complex<double>& a) { a *= fac; });
  }

  auto zero = make_shared<ZMatrix>(geom_->nbasis()*4, geom_->nbasis()*4);
  return make_shared<const DFock>(geom_, zero, natorb, gaunt_, breit_, /*store half*/false, /*robust*/breit_);
}




