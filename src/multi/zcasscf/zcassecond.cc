//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond.cc
// Copyright (C) 2016 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
// Maintainer: Shiozaki group
//
// This file is part of the BAGEL package.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <src/util/math/aughess.h>
#include <src/scf/dhf/dfock.h>
#include <src/multi/zcasscf/zqvec.h>
#include <src/multi/zcasscf/zcassecond.h>

using namespace std;
using namespace bagel;

void ZCASSecond::compute() {
  assert(nvirt_ && nact_);
  Timer timer;

  mute_stdcout();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    {   
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      fci_->compute();
      fci_->compute_rdm12();
//    auto natorb = fci_->natorb_convert();
//    coeff_ = update_coeff(coeff_, natorb.first);
//    occup_ = natorb.second;
      fci_time.tick_print("FCI and RDMs");
      energy_ = fci_->energy();
    }

    shared_ptr<const ZMatrix> cfockao = fci_->jop()->core_fock();
    shared_ptr<const ZMatrix> afockao = compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av());
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const ZMatrix> afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const ZMatrix> qxr = make_shared<ZQvec>(nbasis_*2, nact_, geom_, coeff_, coeff_->slice_copy(nclosed_*2,nocc_*2), nclosed_, fci_, gaunt_, breit_)->get_conjg();

    shared_ptr<ZRotFile> grad = compute_gradient(cfock, afock, qxr);
    kramers_adapt(grad, nclosed_, nact_, nvirt_);

    shared_ptr<ZRotFile> gradele = grad->copy();
    zero_positronic_elements(gradele);

    // check gradient and break if converged
    const double gradient = grad->rms();
    resume_stdcout();
    print_iteration(iter, energy_, gradient, timer.tick());
    mute_stdcout();
    if (gradient < thresh_) {
      resume_stdcout();
      cout << endl << "    * Second-order optimization converged. *   " << endl << endl;
      mute_stdcout();
      break;
    }

    // half-transformed integrals (with JJ)
    list<shared_ptr<RelDFHalf>> halfc_1j, halfg_1j, halfb_1j;
    if (nclosed_) {
      halfc_1j = dynamic_pointer_cast<const DFock>(cfockao)->half_coulomb();
      halfg_1j = dynamic_pointer_cast<const DFock>(cfockao)->half_gaunt();
      halfb_1j = dynamic_pointer_cast<const DFock>(cfockao)->half_breit();
    }
    list<shared_ptr<const RelDFHalf>> halfc, halfg, halfb;
    for (auto& i : halfc_1j)
      halfc.push_back(i->apply_J());
    // caution... Breit cases are done with one J.
    for (auto& i : halfg_1j)
      halfg.push_back(breit_ ? i : i->apply_J());
    for (auto& i : halfb_1j)
      halfb.push_back(i);

    list<shared_ptr<RelDFHalf>> halfac_0j, halfag_0j, halfab_0j;
    shared_ptr<const ZMatrix> acoeff = coeff_->slice_copy(nclosed_*2, nocc_*2);
    tie(halfac_0j, ignore) = RelMOFile::compute_half(geom_, acoeff, false, false);
    if (gaunt_)
      tie(halfag_0j, halfab_0j) = RelMOFile::compute_half(geom_, acoeff, gaunt_, breit_);

    list<shared_ptr<const RelDFHalf>> halfac, halfag, halfab;
    for (auto& i : halfac_0j)
      halfac.push_back(i->apply_JJ());
    // caution... Breit cases are done with J. Breit has picked up J in compute_half
    for (auto& i : halfag_0j)
      halfag.push_back(breit_ ? i->apply_J() : i->apply_JJ());

    assert(gaunt_ || (halfg.empty() && halfag.empty()));
    assert(breit_ || (halfb.empty() && halfab.empty()));

    // compute denominator...
    shared_ptr<const ZRotFile> denom = compute_denom(cfock, afock, qxr, fci_->rdm1_av());

    // augmented Hessian solver
    AugHess<ZRotFile> solver(max_micro_iter_, grad);
    // dummy solver to computer lamda
    AugHess<ZRotFile> solverele(max_micro_iter_, grad);
    // initial trial vector
    shared_ptr<ZRotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    kramers_adapt(trot, nclosed_, nact_, nvirt_);
    trot->normalize();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      Timer mtimer;
      shared_ptr<ZRotFile> sigma = compute_hess_trial(trot, halfc, halfg, halfb, halfac, halfag, halfab, cfock, afock, qxr);
      kramers_adapt(sigma, nclosed_, nact_, nvirt_);

      shared_ptr<ZRotFile> trotele = trot->copy();
      shared_ptr<ZRotFile> sigmaele = sigma->copy();
      zero_positronic_elements(trotele);
      zero_positronic_elements(sigmaele);
      const double scal = trotele->normalize();
      sigmaele->scale(1.0/scal);

      solverele.update(trotele, sigmaele);
      tuple<double,double> lam = solverele.compute_lambda();
      shared_ptr<ZRotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma, lam);
      kramers_adapt(residual, nclosed_, nact_, nvirt_);
      const double err = residual->norm() / lambda;
      resume_stdcout();
      if (!miter) cout << endl;
      cout << "         res : " << setw(8) << setprecision(2) << scientific << err
           <<       "   lamb: " << setw(8) << setprecision(2) << scientific << lambda
           <<       "   eps : " << setw(8) << setprecision(2) << scientific << epsilon
           <<       "   step: " << setw(8) << setprecision(2) << scientific << stepsize
           << setw(8) << fixed << setprecision(2) << mtimer.tick() << endl;
      mute_stdcout();
      if (err < max(thresh_micro_, stepsize*thresh_microstep_))
        break;

      trot = apply_denom(residual, denom, -epsilon, 1.0/lambda);
      kramers_adapt(trot, nclosed_, nact_, nvirt_);
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    }

    shared_ptr<ZRotFile> sol = solver.civec();
    kramers_adapt(sol, nclosed_, nact_, nvirt_);

    shared_ptr<ZMatrix> a = sol->unpack();
    a->scale(complex<double>(0.0,-1.0));
    VectorB eig(a->ndim());
    a->diagonalize(eig);
    auto atmp = a->copy();
    for (int i = 0; i != a->mdim(); ++i) {
      const complex<double> fac = exp(complex<double>(0.0,eig(i)));
      blas::scale_n(fac, a->element_ptr(0,i), a->ndim());
    }
    const ZMatrix R = *a ^ *atmp;
    coeff_ = make_shared<RelCoeff_Block>(*coeff_ * R, coeff_->nclosed(), coeff_->nact(), coeff_->nvirt_nr(), coeff_->nneg());

    resume_stdcout();
    if (iter == max_iter_-1)
      cout << endl << "    * Max iteration reached during the second optimization. *     " << endl << endl;
    mute_stdcout();
  }
  resume_stdcout();
}


shared_ptr<ZRotFile> ZCASSecond::apply_denom(shared_ptr<const ZRotFile> grad, shared_ptr<const ZRotFile> denom, const double shift, const double scale) const {
  shared_ptr<ZRotFile> out = grad->copy();
  for (int i = 0; i != out->size(); ++i)
    if (abs(denom->data(i)*scale+shift) > 1.0e-12)
      out->data(i) /= denom->data(i)*scale+shift;
  return out;
}


shared_ptr<ZRotFile> ZCASSecond::compute_hess_trial(shared_ptr<const ZRotFile> trot,
                                                    list<shared_ptr<const RelDFHalf>> halfc, list<shared_ptr<const RelDFHalf>> halfg, list<shared_ptr<const RelDFHalf>> halfb,
                                                    list<shared_ptr<const RelDFHalf>> halfac, list<shared_ptr<const RelDFHalf>> halfag, list<shared_ptr<const RelDFHalf>> halfab,
                                                    shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
  shared_ptr<ZRotFile> sigma = trot->clone();

  shared_ptr<const ZMatrix> va = trot->va_mat();
  shared_ptr<      ZMatrix> ca = nclosed_ ? trot->ca_mat()->get_conjg() : nullptr; // CAUTION!
  shared_ptr<const ZMatrix> vc = nclosed_ ? trot->vc_mat() : nullptr;

  shared_ptr<const ZMatrix> fcaa = cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2);
  shared_ptr<const ZMatrix> faaa = afock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2);
  shared_ptr<const ZMatrix> fcva = cfock->get_submatrix(nocc_*2, nclosed_*2, nvirt_*2, nact_*2);
  shared_ptr<const ZMatrix> fava = afock->get_submatrix(nocc_*2, nclosed_*2, nvirt_*2, nact_*2);
  shared_ptr<const ZMatrix> fcvv = cfock->get_submatrix(nocc_*2, nocc_*2, nvirt_*2, nvirt_*2);
  shared_ptr<const ZMatrix> favv = afock->get_submatrix(nocc_*2, nocc_*2, nvirt_*2, nvirt_*2);
  shared_ptr<const ZMatrix> fccc = nclosed_ ? cfock->get_submatrix(0, 0, nclosed_*2, nclosed_*2) : nullptr;
  shared_ptr<const ZMatrix> facc = nclosed_ ? afock->get_submatrix(0, 0, nclosed_*2, nclosed_*2) : nullptr;
  shared_ptr<const ZMatrix> fcca = nclosed_ ? cfock->get_submatrix(0, nclosed_*2, nclosed_*2, nact_*2) : nullptr;
  shared_ptr<const ZMatrix> faca = nclosed_ ? afock->get_submatrix(0, nclosed_*2, nclosed_*2, nact_*2) : nullptr;
  shared_ptr<const ZMatrix> fcvc = nclosed_ ? cfock->get_submatrix(nocc_*2, 0, nvirt_*2, nclosed_*2) : nullptr;
  shared_ptr<const ZMatrix> favc = nclosed_ ? afock->get_submatrix(nocc_*2, 0, nvirt_*2, nclosed_*2) : nullptr;

  shared_ptr<const ZMatrix> ccoeff = coeff_->slice_copy(0, nclosed_*2);
  shared_ptr<const ZMatrix> acoeff = coeff_->slice_copy(nclosed_*2, nocc_*2);
  shared_ptr<const ZMatrix> vcoeff = coeff_->slice_copy(nocc_*2, (nocc_+nvirt_)*2);

  ZMatrix rdm1(nact_*2, nact_*2);
  copy_n(fci_->rdm1_av()->data(), rdm1.size(), rdm1.data());
  blas::conj_n(rdm1.data(), rdm1.size());

  // lambda for computing g(D)
  auto compute_gd = [&,this](list<shared_ptr<const RelDFHalf>> halftc, list<shared_ptr<const RelDFHalf>> halfjjc,
                             shared_ptr<const ZMatrix> pcoeff, shared_ptr<const ZMatrix> tpcoeff,
                             const bool gaunt, const bool breit) {
    // TODO to be updated
    list<shared_ptr<const RelDFHalf>> dum;
    auto dfock =  make_shared<DFock>(geom_, cfock->clone(), pcoeff, tpcoeff,
                                     gaunt, breit, halfjjc, halftc, dum, dum, dum, dum, /*robust*/breit);
    dfock->hermite();
    return dfock;
  };

  // g(t_vc) operator and g(t_ac) operator
  if (nclosed_) {
    auto tcoeff = make_shared<ZMatrix>(*vcoeff * *vc + *acoeff * *ca->transpose_conjg());
    list<shared_ptr<RelDFHalf>> tmp;
    tie(tmp, ignore) = RelMOFile::compute_half(geom_, tcoeff, false, false);
    list<shared_ptr<const RelDFHalf>> halftc;
    for (auto& i : tmp)
      halftc.push_back(i);
    const ZMatrix gt = *compute_gd(halftc, halfc, ccoeff, tcoeff, false, false);
    sigma->ax_plus_y_ca(4.0, *ccoeff % gt * *acoeff);
    sigma->ax_plus_y_vc(4.0, *vcoeff % gt * *ccoeff);
    sigma->ax_plus_y_va(4.0, *vcoeff % gt * *acoeff * rdm1);
    sigma->ax_plus_y_ca(-4.0, *ccoeff % gt * *acoeff * rdm1);
  }

  // g(t_va - t_ca)
  auto tcoeff = make_shared<ZMatrix>(nclosed_ ? (*vcoeff * *va - *ccoeff * *ca) : *vcoeff * *va);
  list<shared_ptr<RelDFHalf>> tmp;
  tie(tmp, ignore) = RelMOFile::compute_half(geom_, tcoeff, false, false);
  list<shared_ptr<const RelDFHalf>> halftac;
  for (auto& i : tmp)
    halftac.push_back(i);
  if (nclosed_) {
    list<shared_ptr<const RelDFHalf>> halftacd;
    shared_ptr<const ZMatrix> rdm1p = rdm1.get_conjg();
    for (auto& i : halftac) {
       shared_ptr<RelDFHalf> tmp = i->copy();
       tmp->rotate_occ(rdm1p);
       halftacd.push_back(tmp);
    }
    const ZMatrix gt = *compute_gd(halftacd, halfac, acoeff, make_shared<ZMatrix>(*tcoeff * rdm1), false, false);
    sigma->ax_plus_y_ca(4.0, *ccoeff % gt * *acoeff);
    sigma->ax_plus_y_vc(4.0, *vcoeff % gt * *ccoeff);
  }
  // terms with Qvec
  {
    shared_ptr<const ZMatrix> qaa = qxr->cut(nclosed_*2, nocc_*2);
    sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
    sigma->ax_plus_y_va(-2.0, *va * *qaa);
    if (nclosed_) {
      shared_ptr<const ZMatrix> qva = qxr->cut(nocc_*2, (nocc_+nvirt_)*2);
      shared_ptr<const ZMatrix> qca = qxr->cut(0, nclosed_*2);
      sigma->ax_plus_y_vc(-2.0, *va ^ *qca);
      sigma->ax_plus_y_va(-2.0, *vc * *qca);
      sigma->ax_plus_y_ca(-2.0, *vc % *qva);
      sigma->ax_plus_y_vc(-2.0, *qva ^ *ca);
      sigma->ax_plus_y_ca(-2.0, *ca ^ *qaa);
      sigma->ax_plus_y_ca(-2.0, *ca * *qaa);
    }
  }
  // compute Q' and Q''
  {
    shared_ptr<const ListRelDFFull> fullaa = RelMOFile::compute_full(acoeff, halfac,  /*apply_j*/false);
    shared_ptr<      ListRelDFFull> fullta = RelMOFile::compute_full(acoeff, halftac, /*apply_j*/false);
    shared_ptr<const ListRelDFFull> fulltas = fullta->swap();
    fullta->ax_plus_y(1.0, fulltas);

    shared_ptr<const ListRelDFFull> fullaaD = fullaa->apply_2rdm(fci_->rdm2_av());
    shared_ptr<const ListRelDFFull> fulltaD = fullta->apply_2rdm(fci_->rdm2_av());
    {
      shared_ptr<const ListRelDFFull> fullva  = RelMOFile::compute_full(vcoeff, halfac,  /*apply_j*/false);
      shared_ptr<const ListRelDFFull> fullvta = RelMOFile::compute_full(vcoeff, halftac, /*apply_j*/false);
      shared_ptr<const ZMatrix> qp  = fullva->form_2index(fulltaD, 1.0, false)->get_conjg();
      shared_ptr<const ZMatrix> qpp = fullvta->form_2index(fullaaD, 1.0, false)->get_conjg();
      sigma->ax_plus_y_va( 4.0, *qp + *qpp);
    }
    if (nclosed_) {
      shared_ptr<const ListRelDFFull> fullca  = RelMOFile::compute_full(ccoeff, halfac,  /*apply_j*/false);
      shared_ptr<const ListRelDFFull> fullcta = RelMOFile::compute_full(ccoeff, halftac, /*apply_j*/false);
      shared_ptr<const ZMatrix> qp  = fullca->form_2index(fulltaD, 1.0, false)->get_conjg();
      shared_ptr<const ZMatrix> qpp = fullcta->form_2index(fullaaD, 1.0, false)->get_conjg();
      sigma->ax_plus_y_ca(-4.0, *qp + *qpp);
    }
  }

  // next 1-electron contribution...
  {
    sigma->ax_plus_y_va( 4.0, *fcvv * *va * rdm1);
    sigma->ax_plus_y_va(-2.0, *va * (rdm1 * *fcaa + *fcaa * rdm1));
    if (nclosed_) {
      sigma->ax_plus_y_ca( 4.0, *ca * (*fcaa + *faaa));
      sigma->ax_plus_y_ca( 4.0, *vc % (*fcva + *fava));
      sigma->ax_plus_y_vc(-4.0, *vc * (*fccc + *facc));
      sigma->ax_plus_y_va(-2.0, *vc * (*fcca + *faca));
      sigma->ax_plus_y_vc(-2.0, *va ^ (*fcca + *faca));
      sigma->ax_plus_y_ca(-2.0, *ca * (rdm1 % *fcaa + *fcaa * rdm1));
      sigma->ax_plus_y_vc( 4.0, (*fcvv + *favv) * *vc);
      sigma->ax_plus_y_ca(-4.0, (*fccc + *facc) * *ca);
      sigma->ax_plus_y_va( 2.0, (*fcvc + *favc) * *ca);
      sigma->ax_plus_y_ca( 2.0, (*fcvc + *favc) % *va);
      sigma->ax_plus_y_vc( 2.0, (*fcva + *fava) ^ *ca);
      sigma->ax_plus_y_vc( 2.0, (*fcva + *fava) ^ *ca);
      sigma->ax_plus_y_ca( 4.0, *fccc * *ca * rdm1);
      sigma->ax_plus_y_ca(-4.0, *fcvc % *va * rdm1);
      sigma->ax_plus_y_va(-4.0, *fcvc * *ca * rdm1);
      sigma->ax_plus_y_vc(-2.0, *fcva * rdm1 ^ *ca);
      sigma->ax_plus_y_vc(-2.0, *va * rdm1 ^ *fcca);
      sigma->ax_plus_y_ca(-2.0, *vc % *fcva * rdm1);
      sigma->ax_plus_y_va(-2.0, *vc * *fcca * rdm1);
    }
  }

  blas::conj_n(sigma->ptr_ca(), nclosed_*nact_*4);
  sigma->scale(0.5);
  return sigma;
}


shared_ptr<ZRotFile> ZCASSecond::compute_gradient(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
  auto sigma = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
  shared_ptr<const ZMatrix> rdm1 = fci_->rdm1_av();
  if (nvirt_ && nclosed_) {
    complex<double>* target = sigma->ptr_vc();
    for (int i = 0; i != nclosed_*2; ++i, target += nvirt_*2) {
      zaxpy_(nvirt_*2, 1.0, cfock->element_ptr(nocc_*2, i), 1, target, 1);
      zaxpy_(nvirt_*2, 1.0, afock->element_ptr(nocc_*2, i), 1, target, 1);
    }
  }

  if (nvirt_ && nact_) {
    zgemm3m_("N", "T", nvirt_*2, nact_*2, nact_*2, 1.0, cfock->element_ptr(nocc_*2, nclosed_*2), cfock->ndim(), rdm1->data(), rdm1->ndim(), 0.0, sigma->ptr_va(), nvirt_*2);
    complex<double>* target = sigma->ptr_va();
    for (int i = 0; i != nact_*2; ++i, target += nvirt_*2) {
      zaxpy_(nvirt_*2, 1.0, qxr->element_ptr(nocc_*2, i), 1, target, 1);
    }
  }

  if (nclosed_ && nact_) {
    auto qxrc = qxr->get_conjg();
    auto afockc = afock->get_conjg();
    auto cfockc = cfock->get_conjg();
    complex<double>* target = sigma->ptr_ca();
    for (int i = 0; i != nact_*2; ++i, target += nclosed_*2) {
      zaxpy_(nclosed_*2, 1.0, afockc->element_ptr(0,nclosed_*2+i), 1, target, 1);
      zaxpy_(nclosed_*2, 1.0, cfockc->element_ptr(0,nclosed_*2+i), 1, target, 1);
      zaxpy_(nclosed_*2, -1.0, qxrc->element_ptr(0, i), 1, target, 1);
    }
    // "T" effectively makes complex conjugate of cfock
    zgemm3m_("T", "N", nclosed_*2, nact_*2, nact_*2, -1.0, cfock->element_ptr(nclosed_*2, 0), cfock->ndim(), rdm1->data(), rdm1->ndim(), 1.0, sigma->ptr_ca(), nclosed_*2);
  }
  *sigma *= 2.0;
  return sigma;
}


// TODO this is an approximate denominator. We will replace this with exact diagonal.
shared_ptr<ZRotFile> ZCASSecond::compute_denom(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1) const {
  auto out = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2);
  auto cfockd = make_shared<ZMatrix>(*cfock->get_submatrix(nclosed_*2, nclosed_*2, nact_*2, nact_*2) * *rdm1->get_conjg());

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

  const double thresh = 1.0e-15;
  for (int i = 0; i != out->size(); ++i)
    if (abs(out->data(i)) < thresh)
      out->data(i) = 1.0;
  return out;
}

