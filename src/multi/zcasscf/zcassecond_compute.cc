//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: zcassecond_compute.cc
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

#include <src/scf/dhf/dfock.h>
#include <src/scf/dhf/population_analysis.h>
#include <src/multi/zcasscf/zqvec.h>
#include <src/multi/zcasscf/zcassecond.h>
#include <src/prop/pseudospin/pseudospin.h>
#include <src/ci/zfci/reljop.h>
#include <src/util/math/quatmatrix.h>
#include <src/util/math/aughess.h>
#include <src/util/io/dfpcmo.h>

using namespace std;
using namespace bagel;

// TODO batch size should be automatically determined by the memory size etc. (like in DFock)
const static int batchsize = 250;

void ZCASSecond_base::compute() {
  assert(nvirt_ && nact_);
  Timer timer;

  const bool only_electrons = idata_->get<bool>("only_electrons", false);

  muffle_->mute();
  for (int iter = 0; iter != max_iter_; ++iter) {

    // first perform CASCI to obtain RDMs
    {
      if (iter) fci_->update(coeff_);
      Timer fci_time(0);
      if (external_rdm_.empty()) {
        fci_->compute();
        fci_->compute_rdm12();
      } else {
        if (iter != 0) throw runtime_error("\"external_rdm\" should be used with maxiter == 1");
        fci_->read_external_rdm12_av(external_rdm_);
      }
      // natural orbitals will be used in microiteration
      trans_natorb();
      fci_time.tick_print("FCI and RDMs");
      energy_ = fci_->energy();
    }

    shared_ptr<const ZMatrix> cfockao = fci_->jop()->core_fock();
    shared_ptr<const ZMatrix> afockao = compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av());
    shared_ptr<const ZMatrix> cfock = make_shared<ZMatrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const ZMatrix> afock = make_shared<ZMatrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const ZMatrix> qxr = make_shared<ZQvec>(nbasis_*2, nact_, geom_, coeff_, coeff_->slice_copy(nclosed_*2,nocc_*2), nclosed_, fci_, gaunt_, breit_);

    shared_ptr<ZRotFile> grad = compute_gradient(cfock, afock, qxr);
    if (only_electrons)
      zero_positronic_elements(grad);
    impose_symmetry(grad);

    // check gradient and break if converged
    const double gradient = grad->rms();
    print_iteration(iter, energy_, gradient, timer.tick());
    if (gradient < thresh_) {
      muffle_->unmute();
      cout << endl << "    * Second-order optimization converged. *   " << endl << endl;
      if (nclosed_)
        dynamic_pointer_cast<const DFock>(cfockao)->discard_half();
      break;
    }

    // half-transformed integrals (with JJ)
    list<shared_ptr<const RelDFHalf>> halfc;
    {
      list<shared_ptr<RelDFHalf>> halfc_1j;
      if (nclosed_)
        halfc_1j = dynamic_pointer_cast<const DFock>(cfockao)->half_coulomb();
      for (auto& i : halfc_1j)
        halfc.push_back(i->apply_J());
      if (nclosed_)
        dynamic_pointer_cast<const DFock>(cfockao)->discard_half();
    }

    list<shared_ptr<const RelDFHalf>> halfac;
    {
      list<shared_ptr<RelDFHalf>> halfac_0j;
      tie(halfac_0j, ignore) = RelJop::compute_half(geom_, coeff_->slice_copy(nclosed_*2, nocc_*2), false, false);
      for (auto& i : halfac_0j)
        halfac.push_back(i->apply_JJ());
    }

    // Fock and Q vector with Coulomb only
    shared_ptr<const ZMatrix> cfock_c = cfock;
    shared_ptr<const ZMatrix> afock_c = afock;
    shared_ptr<const ZMatrix> qxr_c = qxr;
    if (gaunt_) {
      cfock_c = nclosed_ ? make_shared<ZMatrix>(*coeff_ % DFock(geom_, hcore_, coeff_->slice(0, nclosed_*2), false, false, false, false) * *coeff_)
                         : make_shared<ZMatrix>(*coeff_ % *hcore_ * *coeff_);
      afock_c = make_shared<ZMatrix>(*coeff_ % *compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av(), /*coulomb only*/true) * *coeff_);
      qxr_c = make_shared<ZQvec>(nbasis_*2, nact_, geom_, coeff_, coeff_->slice_copy(nclosed_*2,nocc_*2), nclosed_, fci_, false, false);
    }

    // compute denominator...
    shared_ptr<const ZRotFile> denom = compute_denom(cfock_c, afock_c, qxr_c, fci_->rdm1_av());

    // augmented Hessian solver
    AugHess<ZRotFile> solver(max_micro_iter_, grad);
    // initial trial vector
    shared_ptr<ZRotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    impose_symmetry(trot);
    trot->normalize();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      Timer mtimer;
      shared_ptr<ZRotFile> sigma = compute_hess_trial(trot, halfc, halfac, cfock_c, afock_c, qxr_c);
      if (only_electrons)
        zero_positronic_elements(sigma);
      impose_symmetry(sigma);

      shared_ptr<ZRotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma);
      impose_symmetry(residual);
      const double err = residual->norm() / lambda;
      muffle_->unmute();
      if (!miter) cout << endl;
      cout << "         res : " << setw(8) << setprecision(2) << scientific << err
           <<       "   lamb: " << setw(8) << setprecision(2) << scientific << lambda
           <<       "   eps : " << setw(8) << setprecision(2) << scientific << epsilon
           <<       "   step: " << setw(8) << setprecision(2) << scientific << stepsize
           << setw(8) << fixed << setprecision(2) << mtimer.tick() << endl;
      muffle_->mute();
      if (err < max(thresh_micro_, stepsize*thresh_microstep_))
        break;

      trot = apply_denom(residual, denom, -epsilon, 1.0/lambda);
      impose_symmetry(trot);
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    }

    shared_ptr<ZRotFile> sol = solver.civec();
    impose_symmetry(sol);

    shared_ptr<ZMatrix> a = sol->unpack();
    a->scale(complex<double>(0.0,-1.0));
    VectorB eig(a->ndim());
    a->diagonalize(eig);
    auto atmp = a->copy();
    for (int i = 0; i != a->mdim(); ++i) {
      const complex<double> fac = exp(complex<double>(0.0,eig(i)));
      blas::scale_n(fac, a->element_ptr(0,i), a->ndim());
    }
    auto R = make_shared<ZMatrix>(*a ^ *atmp);
    impose_symmetry(R);
    coeff_ = make_shared<ZCoeff_Block>(*coeff_ * *R, nclosed_, nact_, nvirtnr_, nneg_);

    if (iter == max_iter_-1) {
      if (external_rdm_.empty() && !conv_ignore_) {
        throw runtime_error("Max iteration reached during the second-order optimization.");
      } else {
        muffle_->unmute();
        cout << endl << "    * Max iteration reached during the second-order optimization.  Convergence not reached! *   " << endl << endl;
      }
    }

#ifndef DISABLE_SERIALIZATION
    if (restart_cas_) {
      stringstream ss; ss << "zcasscf_" << iter;
      OArchive archive(ss.str());
      tie(coeff_,eig_,eigB_,occup_,occupB_) = semi_canonical_orb(kramers());
      shared_ptr<const Reference> ref = conv_to_ref();
      archive << ref;
    }
#endif
  }

  if (max_iter_ > 0) {
    tie(coeff_,eig_,eigB_,occup_,occupB_) = semi_canonical_orb(kramers());
    if (dfpcmo_) {
      auto scoeff = make_shared<ZCoeff_Striped>(*coeff_->striped_format(), nneg_/2, 0, 0, (nclosed_+nact_+nvirtnr_)*2, true);
      const ZMatrix fock = *fci_->jop()->core_fock() + *compute_active_fock(coeff_->slice(nclosed_*2, nocc_*2), fci_->rdm1_av());
      shared_ptr<const VectorB> eig = (*scoeff % fock * *scoeff).diag().get_real_part();
      DFPCMO dfpcmo(scoeff, eig, energy_[0], (nclosed_+nact_+nvirtnr_)*2, nneg_, scoeff->ndim());
      dfpcmo.print();
    }
  } else {
    eig_ = VectorB(coeff_->mdim());
    eigB_ = VectorB(coeff_->mdim());
    occup_ = VectorB(coeff_->mdim());
    occupB_ = VectorB(coeff_->mdim());
  }

  // update construct Jop from scratch
  muffle_->unmute();
  if (nact_ && external_rdm_.empty()) {
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
  }

  // print out orbital populations, if needed
  if (idata_->get<bool>("pop", false)) {
    Timer pop_timer;
    cout << " " << endl;
    cout << "    * Printing out population analysis to casscf.log" << endl;
    muffle_->mute();
    population_analysis(geom_, coeff_->striped_format()->slice(0, 2*(nclosed_+nact_+nvirtnr_)), overlap_);
    muffle_->unmute();
    pop_timer.tick_print("population analysis");
  }

  // TODO When the Property class is implemented, this should be one
  shared_ptr<const PTree> aniso_data = idata_->get_child_optional("aniso");
  if (aniso_data) {
    if (geom_->magnetism()) {
      cout << "  ** Magnetic anisotropy analysis is currently only available for zero-field calculations; sorry." << endl;
    } else {
      const int nspin = aniso_data->get<int>("nspin", (idata_->get_vector<int>("state", 0)).size()-1);
      Pseudospin ps(nspin, geom_, fci_->conv_to_ciwfn(), aniso_data);
      ps.compute(energy_, coeff_->active_part());
    }
  }
}


shared_ptr<ZRotFile> ZCASSecond_base::apply_denom(shared_ptr<const ZRotFile> grad, shared_ptr<const ZRotFile> denom, const double shift, const double scale) const {
  shared_ptr<ZRotFile> out = grad->copy();
  for (int i = 0; i != out->size(); ++i)
    if (abs(denom->data(i)*scale+shift) > 1.0e-12)
      out->data(i) /= denom->data(i)*scale+shift;
  return out;
}


shared_ptr<ZRotFile> ZCASSecond_base::compute_hess_trial(shared_ptr<const ZRotFile> trot, list<shared_ptr<const RelDFHalf>> halfc, list<shared_ptr<const RelDFHalf>> halfac,
                                                         shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
  shared_ptr<ZRotFile> sigma = trot->clone();

  shared_ptr<const ZMatrix> va = trot->va_mat();
  shared_ptr<const ZMatrix> ca = nclosed_ ? trot->ca_mat()->get_conjg() : nullptr; // CAUTION!
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

  // lambda for computing g(D) - only Coulomb
  auto compute_gd = [&,this](list<shared_ptr<const RelDFHalf>> halftc, list<shared_ptr<const RelDFHalf>> halfjjc,
                             shared_ptr<const ZMatrix> pcoeff, shared_ptr<const ZMatrix> tpcoeff) {
    auto aozero = make_shared<ZMatrix>(4*geom_->nbasis(), 4*geom_->nbasis());
    auto dfock =  make_shared<DFock>(geom_, aozero, pcoeff, tpcoeff, halfjjc, halftc);
    dfock->hermite();
    return dfock;
  };

  // g(t_vc) operator and g(t_ac) operator
  if (nclosed_) {
    auto tcoeff = make_shared<ZMatrix>(*vcoeff * *vc + *acoeff * *ca->transpose_conjg());

    const int noccupied = tcoeff->mdim();
    int nbatch = (noccupied-1) / batchsize+1;
    StaticDist dist(noccupied, nbatch);
    vector<pair<size_t, size_t>> table = dist.atable();

    for (auto& itable : table) {
      shared_ptr<const ZMatrix> tcoeff_slice = (nbatch == 1) ? tcoeff : tcoeff->slice_copy(itable.first, itable.first+itable.second);
      shared_ptr<const ZMatrix> ccoeff_slice = (nbatch == 1) ? ccoeff : ccoeff->slice_copy(itable.first, itable.first+itable.second);

      list<shared_ptr<RelDFHalf>> tmp;
      tie(tmp, ignore) = RelJop::compute_half(geom_, tcoeff_slice, false, false);
      list<shared_ptr<const RelDFHalf>> halftc_slice;
      for (auto& i : tmp)
        halftc_slice.push_back(i);

      list<shared_ptr<const RelDFHalf>> halfc_slice;
      if (nbatch == 1)
        halfc_slice = halfc;
      else
        for (auto& i : halfc)
          halfc_slice.push_back(i->slice_b1(itable.first, itable.second));

      const ZMatrix gt = *compute_gd(halftc_slice, halfc_slice, ccoeff_slice, tcoeff_slice);
      sigma->ax_plus_y_ca(4.0, *ccoeff % gt * *acoeff);
      sigma->ax_plus_y_vc(4.0, *vcoeff % gt * *ccoeff);
      sigma->ax_plus_y_va(4.0, *vcoeff % gt * *acoeff * rdm1);
      sigma->ax_plus_y_ca(-4.0, *ccoeff % gt * *acoeff * rdm1);
    }
  }

  // g(t_va - t_ca)
  auto tcoeff = make_shared<ZMatrix>(nclosed_ ? (*vcoeff * *va - *ccoeff * *ca) : *vcoeff * *va);
  list<shared_ptr<RelDFHalf>> tmp;
  tie(tmp, ignore) = RelJop::compute_half(geom_, tcoeff, false, false);
  list<shared_ptr<const RelDFHalf>> halftac;
  for (auto& i : tmp)
    halftac.push_back(i);
  if (nclosed_) {
    list<shared_ptr<const RelDFHalf>> halftacd;
    shared_ptr<const ZMatrix> rdm1p = rdm1.get_conjg();
    for (auto& i : halftac) {
       shared_ptr<RelDFHalf> tmp = i->copy();
       tmp = tmp->transform_occ(rdm1p);
       halftacd.push_back(tmp);
    }
    const ZMatrix gt = *compute_gd(halftacd, halfac, acoeff, make_shared<ZMatrix>(*tcoeff * rdm1));
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
    shared_ptr<const ListRelDFFull> fullaa = RelJop::compute_full(acoeff, halfac,  /*apply_j*/false);
    shared_ptr<      ListRelDFFull> fullta = RelJop::compute_full(acoeff, halftac, /*apply_j*/false);
    shared_ptr<const ListRelDFFull> fulltas = fullta->swap();
    fullta->ax_plus_y(1.0, fulltas);

    shared_ptr<const ListRelDFFull> fullaaD = fullaa->apply_2rdm(fci_->rdm2_av());
    shared_ptr<const ListRelDFFull> fulltaD = fullta->apply_2rdm(fci_->rdm2_av());
    {
      shared_ptr<const ListRelDFFull> fullva  = RelJop::compute_full(vcoeff, halfac,  /*apply_j*/false);
      shared_ptr<const ListRelDFFull> fullvta = RelJop::compute_full(vcoeff, halftac, /*apply_j*/false);
      shared_ptr<const ZMatrix> qp  = fullva->form_2index(fulltaD, 1.0, false)->get_conjg();
      shared_ptr<const ZMatrix> qpp = fullvta->form_2index(fullaaD, 1.0, false)->get_conjg();
      sigma->ax_plus_y_va( 4.0, *qp + *qpp);
    }
    if (nclosed_) {
      shared_ptr<const ListRelDFFull> fullca  = RelJop::compute_full(ccoeff, halfac,  /*apply_j*/false);
      shared_ptr<const ListRelDFFull> fullcta = RelJop::compute_full(ccoeff, halftac, /*apply_j*/false);
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
      sigma->ax_plus_y_vc( 4.0, (*fcva + *fava) ^ *ca);
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


shared_ptr<ZRotFile> ZCASSecond_base::compute_gradient(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr) const {
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
shared_ptr<ZRotFile> ZCASSecond_base::compute_denom(shared_ptr<const ZMatrix> cfock, shared_ptr<const ZMatrix> afock, shared_ptr<const ZMatrix> qxr, shared_ptr<const ZMatrix> rdm1) const {
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


void ZCASSecond::trans_natorb() {
  // first make natural orbitals
  shared_ptr<ZMatrix> rdm1 = fci_->rdm1_av();
  rdm1->scale(-1.0);
  rdm1->add_diag(1.0);

  VectorB occup(nact_*2);
  auto natorb = make_shared<QuatMatrix>(*rdm1);
  natorb->diagonalize(occup);
  for (int i = 0; i != nact_; ++i)
    occup[i] = 1.0-occup[i];

  // rotate RDMs stored in FCI object
  fci_->rotate_rdms(natorb);

  if (natocc_) {
    cout << "  ========       state-averaged       ======== " << endl;
    cout << "  ======== natural occupation numbers ======== " << endl;
    for (int i = 0; i != nact_; ++i)
      cout << setprecision(4) << "   Orbital " << i << " : " << (occup[i] < numerical_zero__ ? 0.0 : occup[i]) << endl;
    cout << "  ============================================ " << endl;
  }

  // D_rs = C*_ri D_ij (C*_rj)^+. Dij = U_ik L_k (U_jk)^+. So, C'_ri = C_ri * U*_ik ; hence conjugation needed
  auto cnew = make_shared<ZCoeff_Block>(*coeff_, nclosed_, nact_, nvirtnr_, nneg_);
  cnew->copy_block(0, nclosed_*2, cnew->ndim(), nact_*2, coeff_->slice(nclosed_*2, nocc_*2) * *natorb->get_conjg());
  coeff_ = cnew;
}


void ZCASSecond_London::trans_natorb() {
  // first make natural orbitals
  shared_ptr<ZMatrix> natorb = fci_->rdm1_av();
  natorb->scale(-1.0);
  natorb->add_diag(1.0);

  VectorB occup(nact_*2);
  natorb->diagonalize(occup);
  for (int i = 0; i != nact_*2; ++i)
    occup[i] = 1.0-occup[i];

  // reorder the orbitals into the block format
  shared_ptr<ZMatrix> trans = natorb->clone();
  VectorB roccup(nact_*2);
  for (int i = 0; i != nact_; ++i) {
    copy_n(natorb->element_ptr(0, 2*i), nact_*2, trans->element_ptr(0, i));
    copy_n(natorb->element_ptr(0, 2*i+1), nact_*2, trans->element_ptr(0, nact_+i));
    roccup[i] = occup[2*i];
    roccup[nact_+i] = occup[2*i+1];
  }

  // rotate RDMs stored in FCI object
  fci_->rotate_rdms(trans);

  if (natocc_) {
    cout << "  ========       state-averaged       ======== " << endl;
    cout << "  ======== natural occupation numbers ======== " << endl;
    for (int i = 0; i != roccup.size(); ++i)
      cout << setprecision(4) << "   Orbital " << i << " : " << (roccup[i] < numerical_zero__ ? 0.0 : roccup[i]) << endl;
    cout << "  ============================================ " << endl;
  }

  // D_rs = C*_ri D_ij (C*_rj)^+. Dij = U_ik L_k (U_jk)^+. So, C'_ri = C_ri * U*_ik ; hence conjugation needed
  auto cnew = make_shared<ZCoeff_Block>(*coeff_, nclosed_, nact_, nvirtnr_, nneg_);
  cnew->copy_block(0, nclosed_*2, cnew->ndim(), nact_*2, coeff_->slice(nclosed_*2, nocc_*2) * *trans->get_conjg());
  coeff_ = cnew;
}
