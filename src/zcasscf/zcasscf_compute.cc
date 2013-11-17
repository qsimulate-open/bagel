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
  shared_ptr<const ZMatrix> overlap = make_shared<RelOverlap>(geom_);

  // coeff_ is a kramers adapted coefficient..
  {
    array<shared_ptr<const ZMatrix>,2> tmp = RelMOFile::kramers(coeff_, overlap, hcore);
    auto ctmp = coeff_->clone();
    int i = 0;
    ctmp->copy_block(0, i, coeff_->ndim(), nclosed_, tmp[0]->slice(0,nclosed_)); i += nclosed_;
    ctmp->copy_block(0, i, coeff_->ndim(), nclosed_, tmp[1]->slice(0,nclosed_)); i += nclosed_;
    ctmp->copy_block(0, i, coeff_->ndim(), nact_, tmp[0]->slice(nclosed_, nocc_)); i += nact_;
    ctmp->copy_block(0, i, coeff_->ndim(), nact_, tmp[1]->slice(nclosed_, nocc_)); i += nact_;
    ctmp->copy_block(0, i, coeff_->ndim(), nvirt_, tmp[0]->slice(nocc_, nocc_+nvirt_)); i += nvirt_;
    ctmp->copy_block(0, i, coeff_->ndim(), nvirt_, tmp[1]->slice(nocc_, nocc_+nvirt_));
    coeff_ = ctmp;
  }

  for (int iter = 0; iter != max_iter_; ++iter) {
    // first perform CASCI to obtain RDMs
    if (nact_) {
      mute_stdcout();
      if (iter) fci_->update(coeff_);
      fci_->compute();
      fci_->compute_rdm12();
      resume_stdcout();
    }

    // calculate 1RDM in an original basis set
    shared_ptr<const ZMatrix> rdm1 = nact_ ? transform_rdm1() : shared_ptr<const ZMatrix>();

    // closed Fock operator
    shared_ptr<const ZMatrix> cfockao = nclosed_ ? make_shared<const DFock>(geom_, hcore, coeff_->slice(0,nclosed_*2), gaunt_, breit_, /*store half*/false, /*robust*/breit_) : hcore;
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);

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

    auto xlog = make_shared<ZRotFile>(x->log(4), nclosed_*2, nact_*2, nvirt_*2, /*superci*/ false);
    shared_ptr<ZMatrix> amat = bfgs->extrapolate(grad, xlog)->unpack<ZMatrix>();

    const double gradient = amat->rms();

    // multiply -1 from the formula. multiply -i to make amat hermite (will be compensated)
    *amat *= -1.0 * complex<double>(0.0, -1.0);

    // restore the matrix from RotFile
    unique_ptr<double[]> teig(new double[amat->ndim()]);
    amat->diagonalize(teig.get());
    auto amat_sav = amat->copy();
    for (int i = 0; i != amat->ndim(); ++i) {
      complex<double> ex = exp(complex<double>(0.0, teig[i]));
      transform(amat->element_ptr(0,i), amat->element_ptr(0,i+1), amat->element_ptr(0,i), [&ex](complex<double> a) { return a*ex; });
    }
    auto expa = make_shared<const ZMatrix>(*amat ^ *amat_sav);

    coeff_ = make_shared<const ZMatrix>(*coeff_**expa);
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
  // TODO compute RelOverlap only once (this is comptued also in qzvec)
  auto overlap = make_shared<const RelOverlap>(geom_);
  shared_ptr<const ZMatrix> ocoeff = coeff_->slice(nclosed_*2, nclosed_*2+nact_*2)->get_conjg();
  const ZMatrix co = *ocoeff % *overlap * *coeff_tot;
  return make_shared<ZMatrix>(co * *rdm1_tot ^ co);
}


shared_ptr<const ZMatrix> ZCASSCF::active_fock(shared_ptr<const ZMatrix> rdm1) const {
  // form natural orbitals
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
    transform(natorb->element_ptr(0, i), natorb->element_ptr(0, i+1), natorb->element_ptr(0, i), [&fac](complex<double> a) { return fac*a; });
  }

  auto zero = make_shared<ZMatrix>(geom_->nbasis()*4, geom_->nbasis()*4);
  return make_shared<const DFock>(geom_, zero, natorb, gaunt_, breit_, /*store half*/false, /*robust*/breit_);
}


shared_ptr<const ZRotFile> ZCASSCF::compute_denom(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1) const {
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2, /*superci*/false);

  shared_ptr<ZMatrix> cfockd;
  if (nact_) {
    cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1);
    // TODO check
    cfockd->hermite();
  }

  // ia part (4.7a)
  if (nvirt_ && nclosed_) {
    complex<double>* target = out->ptr_vc();
    for (int i = 0; i != nclosed_*2; ++i) {
      for (int j = 0; j != nvirt_*2; ++j) {
        *target++ = cfock->element(j+nocc_*2, j+nocc_*2) + afock->element(j+nocc_*2, j+nocc_*2) - cfock->element(i,i) - afock->element(i,i);
      }
    }
  }
  // ra part (4.7b)
  if (nvirt_ && nact_) {
    complex<double>* target = out->ptr_va();
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nvirt_*2; ++j) {
        *target++ = rdm1->element(i, i)*(cfock->element(j+nocc_*2, j+nocc_*2)+afock->element(j+nocc_*2, j+nocc_*2))
                  - cfockd->element(i, i) - qxr->element(i+nclosed_*2, i);
      }
    }
  }
  // it part (4.7c)
  if (nclosed_ && nact_) {
    complex<double>* target = out->ptr_ca();
    for (int i = 0; i != nact_*2; ++i) {
      for (int j = 0; j != nclosed_*2; ++j) {
        *target++ = (cfock->element(i+nclosed_*2, i+nclosed_*2)+afock->element(i+nclosed_*2, i+nclosed_*2) - cfock->element(j,j) - afock->element(j,j))
                  + rdm1->element(i,i)*(cfock->element(j,j)+afock->element(j,j)) - cfockd->element(i,i) - qxr->element(i+nclosed_*2, i);
      }
    }
  }

  const double thresh = 1.0e-8;
  for (int i = 0; i != out->size(); ++i)
    if (fabs(out->data(i)) < thresh) {
      out->data(i) = 1.0e10;
    }
  return out;
}



// grad(a/i) (eq.4.3a): (cfock_ai+afock_ai)
void ZCASSCF::grad_vc(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<ZRotFile> sigma) const {
  if (!nvirt_ || !nclosed_) return;
  complex<double>* target = sigma->ptr_vc();
  for (int i = 0; i != nclosed_*2; ++i, target += nvirt_*2) {
    zaxpy_(nvirt_*2, 1.0, cfock->element_ptr(nocc_*2, i), 1, target, 1);
    zaxpy_(nvirt_*2, 1.0, afock->element_ptr(nocc_*2, i), 1, target, 1);
  }

  // symmetry adaptation
  for (int i = 0; i != nclosed_; ++i) {
    for (int j = 0; j != nvirt_; ++j) {
      sigma->ele_vc(j, i) = (sigma->ele_vc(j, i) + conj(sigma->ele_vc(j+nvirt_, i+nclosed_))) * 0.5;
      sigma->ele_vc(j+nvirt_, i+nclosed_) = conj(sigma->ele_vc(j, i));

      sigma->ele_vc(j+nvirt_, i) = (sigma->ele_vc(j+nvirt_, i) - conj(sigma->ele_vc(j, i+nclosed_))) * 0.5;
      sigma->ele_vc(j, i+nclosed_) = - conj(sigma->ele_vc(j+nvirt_, i));
    }
  }
}


// grad(a/t) (eq.4.3b): cfock_au gamma_ut + q_at
void ZCASSCF::grad_va(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1, shared_ptr<ZRotFile> sigma) const {
  if (!nvirt_ || !nact_) return;
  // TODO not sure about complex conjugation of rdm1
  zgemm3m_("N", "T", nvirt_*2, nact_*2, nact_*2, 1.0, cfock->element_ptr(nocc_*2, nclosed_*2), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_*2);
  complex<double>* target = sigma->ptr_va();
  for (int i = 0; i != nact_*2; ++i, target += nvirt_*2) {
    zaxpy_(nvirt_*2, 1.0, qxr->element_ptr(nocc_*2, i), 1, target, 1);
  }
  // symmetry adaptation
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nvirt_; ++j) {
      sigma->ele_va(j, i) = (sigma->ele_va(j, i) + conj(sigma->ele_va(j+nvirt_, i+nact_))) * 0.5;
      sigma->ele_va(j+nvirt_, i+nact_) = conj(sigma->ele_va(j, i));

      sigma->ele_va(j+nvirt_, i) = (sigma->ele_va(j+nvirt_, i) - conj(sigma->ele_va(j, i+nact_))) * 0.5;
      sigma->ele_va(j, i+nact_) = - conj(sigma->ele_va(j+nvirt_, i));
    }
  }
}


// grad(r/i) (eq.4.3c): (cfock_ri+afock_ri) - cfock_iu gamma_ur - qxr_ir
void ZCASSCF::grad_ca(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1, shared_ptr<ZRotFile> sigma) const {
  if (!nclosed_ || !nact_) return;
  // TODO check
  auto qxrc = qxr->get_conjg();
  complex<double>* target = sigma->ptr_ca();
  for (int i = 0; i != nact_*2; ++i, target += nclosed_*2) {
    zaxpy_(nclosed_*2, 1.0, afock->element_ptr(0,nclosed_*2+i), 1, target, 1);
    zaxpy_(nclosed_*2, 1.0, cfock->element_ptr(0,nclosed_*2+i), 1, target, 1);
    zaxpy_(nclosed_*2, -1.0, qxrc->element_ptr(0, i), 1, target, 1);
  }
  // "T" effectively makes complex conjugate of cfock
  zgemm3m_("T", "N", nclosed_*2, nact_*2, nact_*2, -1.0, cfock->element_ptr(nclosed_*2, 0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, sigma->ptr_ca(), nclosed_*2);

  // symmetry adaptation
  for (int i = 0; i != nact_; ++i) {
    for (int j = 0; j != nclosed_; ++j) {
      sigma->ele_ca(j, i) = (sigma->ele_ca(j, i) + conj(sigma->ele_ca(j+nclosed_, i+nact_))) * 0.5;
      sigma->ele_ca(j+nclosed_, i+nact_) = conj(sigma->ele_ca(j, i));

      sigma->ele_ca(j+nclosed_, i) = (sigma->ele_ca(j+nclosed_, i) - conj(sigma->ele_ca(j, i+nact_))) * 0.5;
      sigma->ele_ca(j, i+nact_) = - conj(sigma->ele_ca(j+nclosed_, i));
    }
  }
}
