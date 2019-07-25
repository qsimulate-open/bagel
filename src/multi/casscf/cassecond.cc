//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cassecond.cc
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
#include <src/scf/hf/fock.h>
#include <src/multi/casscf/qvec.h>
#include <src/multi/casscf/cassecond.h>
#include <src/prop/hyperfine.h>
#include <src/prop/multipole.h>

using namespace std;
using namespace bagel;

void CASSecond::compute() {
  assert(nvirt_ && nact_);
  Timer timer;

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
        if (iter != 0)
          throw runtime_error("\"external_rdm\" should be used with maxiter == 1");
        fci_->read_external_rdm12_av(external_rdm_);
      }
      trans_natorb();
      fci_time.tick_print("FCI and RDMs");
      energy_ = fci_->energy();
    }

    shared_ptr<const Matrix> cfockao = fci_->jop()->core_fock();
    shared_ptr<const Matrix> afockao = compute_active_fock(coeff_->slice(nclosed_, nocc_), fci_->rdm1_av());
    shared_ptr<const Matrix> cfock = make_shared<Matrix>(*coeff_ % *cfockao * *coeff_);
    shared_ptr<const Matrix> afock = make_shared<Matrix>(*coeff_ % *afockao * *coeff_);
    shared_ptr<const Qvec> qxr = make_shared<Qvec>(coeff_->mdim(), nact_, coeff_, nclosed_, fci_, fci_->rdm2_av());

    shared_ptr<const RotFile> grad = compute_gradient(cfock, afock, qxr);

    // check gradient and break if converged
    const double gradient = grad->rms();
    print_iteration(iter, energy_, gradient, timer.tick());
    if (gradient < thresh_) {
      muffle_->unmute();
      cout << endl << "    * Second-order optimization converged. *   " << endl << endl;
      break;
    }

    // half-transformed integrals (with JJ)
    shared_ptr<const DFHalfDist> half_1j = nclosed_ ? dynamic_pointer_cast<const Fock<1>>(cfockao)->half() : nullptr;
    shared_ptr<const DFHalfDist> half = nclosed_ ? half_1j->apply_J() : nullptr;
    shared_ptr<const DFHalfDist> halfa = fci_->jop()->mo2e_1ext()->apply_JJ();

    // compute denominator...
    shared_ptr<const RotFile> denom = compute_denom(half, half_1j, halfa, cfock, afock);

    AugHess<RotFile> solver(max_micro_iter_, grad);
    // initial trial vector
    shared_ptr<RotFile> trot = apply_denom(grad, denom, 0.001, 1.0);
    trot->normalize();

    for (int miter = 0; miter != max_micro_iter_; ++miter) {
      Timer mtimer;
      shared_ptr<const RotFile> sigma = compute_hess_trial(trot, half, halfa, cfock, afock, qxr);
      shared_ptr<const RotFile> residual;
      double lambda, epsilon, stepsize;
      tie(residual, lambda, epsilon, stepsize) = solver.compute_residual(trot, sigma);
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
      for (int i = 0; i != 10; ++i) {
        const double norm = solver.orthog(trot);
        if (norm > 0.25) break;
      }
    }

    shared_ptr<const RotFile> sol = solver.civec();
    shared_ptr<const Matrix> a = sol->unpack();
    Matrix w(*a * *a);
    VectorB eig(a->ndim());
    w.diagonalize(eig);
    Matrix wc(w);
    Matrix ws(w);
    for (int i = 0; i != a->ndim(); ++i) {
      const double tau = sqrt(fabs(eig(i)));
      blas::scale_n(cos(tau), wc.element_ptr(0,i), wc.ndim());
      blas::scale_n(tau > 1.0e-15 ? sin(tau)/tau : 1.0, ws.element_ptr(0,i), ws.ndim());
    }
    const Matrix R = (wc ^ w) + (ws ^ w) * *a;

    coeff_ = make_shared<Coeff>(*coeff_ * R);

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
      stringstream ss; ss << "casscf_" << iter;
      OArchive archive(ss.str());
      auto ref = make_shared<const Reference>(geom_, coeff_, nclosed_, nact_, nvirt_, energy_);
      archive << ref;
    }
#endif
  }
  muffle_->unmute();

  // block diagonalize coeff_ in nclosed and nvirt
  if (max_iter_ > 0) {
    auto tmp = semi_canonical_orb();
    coeff_ = get<0>(tmp);
    eig_   = get<1>(tmp);
    occup_ = get<2>(tmp);
  }

  // this is not needed for energy, but for consistency we want to have this...
  // update construct Jop from scratch
  if (nact_ && external_rdm_.empty()) {
    fci_->update(coeff_);
    fci_->compute();
    fci_->compute_rdm12();
  }

  if (nstate_ == 1) {
    const MatView ocoeff = coeff_->slice(0, nocc_);
    Dipole dipole(geom_, make_shared<Matrix>(ocoeff * *fci_->rdm1(0)->rdm1_mat(nclosed_) ^ ocoeff));
    dipole.compute();
  }

  // calculate the HFCCs
  if (do_hyperfine_ && !geom_->external() && nstate_ == 1) {
    HyperFine hfcc(geom_, spin_density(), fci_->det()->nspin(), "CASSCF");
    hfcc.compute();
  }
}


shared_ptr<RotFile> CASSecond::apply_denom(shared_ptr<const RotFile> grad, shared_ptr<const RotFile> denom, const double shift, const double scale) const {
  shared_ptr<RotFile> out = grad->copy();
  for (int i = 0; i != out->size(); ++i)
    if (fabs(denom->data(i)*scale+shift) > 1.0e-12)
      out->data(i) /= denom->data(i)*scale+shift;
  return out;
}


shared_ptr<RotFile> CASSecond::compute_gradient(shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  auto sigma = make_shared<RotFile>(nclosed_, nact_, nvirt_);
  shared_ptr<const RDM<1>> rdm1 = fci_->rdm1_av();
  if (nclosed_) {
    double* target = sigma->ptr_vc();
    for (int i = 0; i != nclosed_; ++i, target += nvirt_) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(nocc_,i), nvirt_, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(nocc_,i), nvirt_, target);
    }
  }
  {
    double* target = sigma->ptr_va();
    for (int i = 0; i != nact_; ++i, target += nvirt_) {
      blas::ax_plus_y_n(2.0, qxr->element_ptr(nocc_, i), nvirt_, target);
      for (int j = 0; j != nact_; ++j)
        blas::ax_plus_y_n(2.0*rdm1->element(j,i), cfock->element_ptr(nocc_, nclosed_+j), nvirt_, target);
    }
  }
  if (nclosed_) {
    double* target = sigma->ptr_ca();
    for (int i = 0; i != nact_; ++i, target += nclosed_) {
      blas::ax_plus_y_n(4.0, cfock->element_ptr(0,nclosed_+i), nclosed_, target);
      blas::ax_plus_y_n(4.0, afock->element_ptr(0,nclosed_+i), nclosed_, target);
      blas::ax_plus_y_n(-2.0, qxr->element_ptr(0, i), nclosed_, target);
      for (int j = 0; j != nact_; ++j)
        blas::ax_plus_y_n(-2.0*rdm1->element(j,i), cfock->element_ptr(0,nclosed_+j), nclosed_, target);
    }
  }
  return sigma;
}


shared_ptr<RotFile> CASSecond::compute_denom(shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> half_1j, shared_ptr<const DFHalfDist> halfa,
                                             shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock) const {
  auto denom = make_shared<RotFile>(nclosed_, nact_, nvirt_);
  const MatView ccoeff = coeff_->slice(0, nclosed_);
  const MatView acoeff = coeff_->slice(nclosed_, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  Matrix rdm1(nact_, nact_);
  copy_n(fci_->rdm1_av()->data(), nact_*nact_, rdm1.data());
  {
    const Matrix fcd = *cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_) * rdm1;
    const Matrix fock = *cfock + *afock;
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nclosed_; ++j)
        denom->ele_ca(j, i) += 4.0 * fock(i+nclosed_, i+nclosed_) - 4.0 * fock(j, j) - 2.0 * fcd(i, i) + 2.0 * (*cfock)(j, j) * rdm1(i, i);
    for (int i = 0; i != nclosed_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_vc(j, i) += 4.0 * fock(j+nocc_, j+nocc_) - 4.0 * fock(i, i);
    for (int i = 0; i != nact_; ++i)
      for (int j = 0; j != nvirt_; ++j)
        denom->ele_va(j, i) += 2.0 * (*cfock)(j+nocc_, j+nocc_) * rdm1(i, i) - 2.0 * fcd(i, i);

    const int nao = coeff_->ndim();
    if (nclosed_) {
      auto vvc = half_1j->compute_second_transform(vcoeff)->form_4index_diagonal()->transpose();
      denom->ax_plus_y_vc(12.0, *vvc);

      shared_ptr<const DFFullDist> vgcc = half->compute_second_transform(ccoeff);
      const int nri = vgcc->block(0)->asize();
      Matrix tmp(nao, nao);
      for (int i = 0; i != nclosed_; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgcc->block(0)->data()+nri*(i+nclosed_*i), 1, 0.0, tmp.data(), 1);
        if (!vgcc->serial())
          tmp.allreduce();
        Matrix tmp0 = vcoeff % tmp * vcoeff;
        blas::ax_plus_y_n(-4.0, tmp0.diag().data(), nvirt_, denom->ptr_vc()+nvirt_*i);
      }
    }
    shared_ptr<const DFFullDist> vaa  = halfa->compute_second_transform(acoeff);
    const int nri = vaa->block(0)->asize();
    {
      shared_ptr<const DFFullDist> vgaa = vaa->apply_2rdm(*fci_->rdm2_av());
      Matrix tmp(nao, nao);
      for (int i = 0; i != nact_; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact_*i), 1, 0.0, tmp.data(), 1);
        if (!vgaa->serial())
          tmp.allreduce();
        Matrix tmp0 = vcoeff % tmp * vcoeff;
        blas::ax_plus_y_n(2.0, tmp0.diag().data(), nvirt_, denom->ptr_va()+nvirt_*i);
        if (nclosed_) {
          Matrix tmp1 = ccoeff % tmp * ccoeff;
          blas::ax_plus_y_n(2.0, tmp1.diag().data(), nclosed_, denom->ptr_ca()+nclosed_*i);
        }
      }
      shared_ptr<const DFFullDist> vaa_noj = fci_->jop()->mo2e_1ext()->compute_second_transform(acoeff);
      shared_ptr<const Matrix> mo2e = vaa->form_4index(vaa_noj, 1.0);
      for (int i = 0; i != nact_; ++i) {
        const double e2 = -2.0 * blas::dot_product(mo2e->element_ptr(0, i*nact_), nact_*nact_*nact_, fci_->rdm2_av()->element_ptr(0,0,0,i));
        for (int j = 0; j != nvirt_; ++j)
          denom->ele_va(j, i) += e2;
        for (int j = 0; j != nclosed_; ++j)
          denom->ele_ca(j, i) += e2;
      }

      Matrix rdmk(nact_*nact_, nact_);
      for (int i = 0; i != nact_; ++i)
        for (int j = 0; j != nact_; ++j)
          for (int k = 0; k != nact_; ++k)
            rdmk(k+nact_*j, i) = fci_->rdm2_av()->element(k, i, j, i) + fci_->rdm2_av()->element(k, i, i, j);
      shared_ptr<const DFFullDist> vav = fci_->jop()->mo2e_1ext()->compute_second_transform(vcoeff)->apply_J();
      denom->ax_plus_y_va(2.0, *(rdmk % *vav->form_4index_diagonal_part()).transpose());
      if (nclosed_) {
        shared_ptr<const DFFullDist> vac = fci_->jop()->mo2e_1ext()->compute_second_transform(ccoeff)->apply_J();
        shared_ptr<const Matrix> mcaa = vac->form_4index_diagonal_part()->transpose();
        denom->ax_plus_y_ca(2.0, *mcaa * rdmk);
        shared_ptr<Matrix> mcaad = mcaa->copy();
        dgemm_("N", "N", nclosed_*nact_, nact_, nact_, -1.0, mcaa->data(), nclosed_*nact_, rdm1.data(), nact_, 1.0, mcaad->data(), nclosed_*nact_);
        for (int i = 0; i != nact_; ++i)
          blas::ax_plus_y_n(12.0, mcaad->element_ptr(0, i+nact_*i), nclosed_, denom->ptr_ca()+i*nclosed_);
      }
    }
    if (nclosed_) {
      Matrix tmp(nao, nao);
      shared_ptr<DFFullDist> vgaa = vaa->copy();
      vgaa = vgaa->transform_occ1(make_shared<Matrix>(rdm1));
      vgaa->ax_plus_y(-1.0, vaa);
      for (int i = 0; i != nact_; ++i) {
        dgemv_("T", nri, nao*nao, 1.0, geom_->df()->block(0)->data(), nri, vgaa->block(0)->data()+nri*(i+nact_*i), 1, 0.0, tmp.data(), 1);
        if (!vgaa->serial())
          tmp.allreduce();
        Matrix tmp0 = ccoeff % tmp * ccoeff;
        blas::ax_plus_y_n(4.0, tmp0.diag().data(), nclosed_, denom->ptr_ca()+nclosed_*i);
      }
    }
  }
  return denom;
}


shared_ptr<RotFile> CASSecond::compute_hess_trial(shared_ptr<const RotFile> trot, shared_ptr<const DFHalfDist> half, shared_ptr<const DFHalfDist> halfa,
                                                  shared_ptr<const Matrix> cfock, shared_ptr<const Matrix> afock, shared_ptr<const Matrix> qxr) const {
  shared_ptr<RotFile> sigma = trot->clone();

  shared_ptr<const Matrix> va = trot->va_mat();
  shared_ptr<const Matrix> ca = nclosed_ ? trot->ca_mat() : nullptr;
  shared_ptr<const Matrix> vc = nclosed_ ? trot->vc_mat() : nullptr;

  shared_ptr<const Matrix> fcaa = cfock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
  shared_ptr<const Matrix> faaa = afock->get_submatrix(nclosed_, nclosed_, nact_, nact_);
  shared_ptr<const Matrix> fcva = cfock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
  shared_ptr<const Matrix> fava = afock->get_submatrix(nocc_, nclosed_, nvirt_, nact_);
  shared_ptr<const Matrix> fcvv = cfock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
  shared_ptr<const Matrix> favv = afock->get_submatrix(nocc_, nocc_, nvirt_, nvirt_);
  shared_ptr<const Matrix> fccc = nclosed_ ? cfock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
  shared_ptr<const Matrix> facc = nclosed_ ? afock->get_submatrix(0, 0, nclosed_, nclosed_) : nullptr;
  shared_ptr<const Matrix> fcca = nclosed_ ? cfock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
  shared_ptr<const Matrix> faca = nclosed_ ? afock->get_submatrix(0, nclosed_, nclosed_, nact_) : nullptr;
  shared_ptr<const Matrix> fcvc = nclosed_ ? cfock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;
  shared_ptr<const Matrix> favc = nclosed_ ? afock->get_submatrix(nocc_, 0, nvirt_, nclosed_) : nullptr;

  const MatView ccoeff = coeff_->slice(0, nclosed_);
  const MatView acoeff = coeff_->slice(nclosed_, nocc_);
  const MatView vcoeff = coeff_->slice(nocc_, nocc_+nvirt_);

  Matrix rdm1(nact_, nact_);
  copy_n(fci_->rdm1_av()->data(), nact_*nact_, rdm1.data());

  // lambda for computing g(D)
  auto compute_gd = [&,this](shared_ptr<const DFHalfDist> halft, shared_ptr<const DFHalfDist> halfjj, const MatView pcoeff) {
    shared_ptr<const Matrix> pcoefft = make_shared<Matrix>(pcoeff)->transpose();
    shared_ptr<Matrix> gd = geom_->df()->compute_Jop(halft, pcoefft);
    shared_ptr<Matrix> ex0 = halfjj->form_2index(halft, 1.0);
    ex0->symmetrize();
    gd->ax_plus_y(-0.5, ex0);
    return gd;
  };

  // g(t_vc) operator and g(t_ac) operator
  if (nclosed_) {
    const Matrix tcoeff = vcoeff * *vc + acoeff * *ca->transpose();
    auto halft = geom_->df()->compute_half_transform(tcoeff);
    const Matrix gt = *compute_gd(halft, half, ccoeff);
    sigma->ax_plus_y_ca(32.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_vc(32.0, vcoeff % gt * ccoeff);
    sigma->ax_plus_y_va(16.0, vcoeff % gt * acoeff * rdm1);
    sigma->ax_plus_y_ca(-16.0, ccoeff % gt * acoeff * rdm1);
  }
  // g(t_va - t_ca)
  const Matrix tcoeff = nclosed_ ? (vcoeff * *va - ccoeff * *ca) : vcoeff * *va;
  shared_ptr<const DFHalfDist> halfta = geom_->df()->compute_half_transform(tcoeff);
  if (nclosed_) {
    shared_ptr<DFHalfDist> halftad = halfta->copy();
    halftad = halftad->transform_occ(make_shared<Matrix>(rdm1));
    const Matrix gt = *compute_gd(halftad, halfa, acoeff);
    sigma->ax_plus_y_ca(16.0, ccoeff % gt * acoeff);
    sigma->ax_plus_y_vc(16.0, vcoeff % gt * ccoeff);
  }
  // terms with Qvec
  {
    shared_ptr<const Matrix> qaa = qxr->cut(nclosed_, nocc_);
    sigma->ax_plus_y_va(-2.0, *va ^ *qaa);
    sigma->ax_plus_y_va(-2.0, *va * *qaa);
    if (nclosed_) {
      shared_ptr<const Matrix> qva = qxr->cut(nocc_, nocc_+nvirt_);
      shared_ptr<const Matrix> qca = qxr->cut(0, nclosed_);
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
    shared_ptr<const DFFullDist> fullaa = halfa->compute_second_transform(acoeff);
    shared_ptr<DFFullDist> fullta = halfta->compute_second_transform(acoeff);
    shared_ptr<const DFFullDist> fulltas = fullta->swap();
    fullta->ax_plus_y(1.0, fulltas);
    shared_ptr<const DFFullDist> fullaaD = fullaa->apply_2rdm(*fci_->rdm2_av());
    shared_ptr<const DFFullDist> fulltaD = fullta->apply_2rdm(*fci_->rdm2_av());
    shared_ptr<const Matrix> qp  = halfa->form_2index(fulltaD, 1.0);
    shared_ptr<const Matrix> qpp = halfta->form_2index(fullaaD, 1.0);

    sigma->ax_plus_y_va( 4.0, vcoeff % (*qp + *qpp));
    if (nclosed_)
      sigma->ax_plus_y_ca(-4.0, ccoeff % (*qp + *qpp));
  }

  // next 1-electron contribution...
  {
    sigma->ax_plus_y_va( 4.0, *fcvv * *va * rdm1);
    sigma->ax_plus_y_va(-2.0, *va * (rdm1 * *fcaa + *fcaa * rdm1));
    if (nclosed_) {
      sigma->ax_plus_y_ca( 8.0, *ca * (*fcaa + *faaa));
      sigma->ax_plus_y_ca( 8.0, *vc % (*fcva + *fava));
      sigma->ax_plus_y_vc(-8.0, *vc * (*fccc + *facc));
      sigma->ax_plus_y_va(-4.0, *vc * (*fcca + *faca));
      sigma->ax_plus_y_vc(-4.0, *va ^ (*fcca + *faca));
      sigma->ax_plus_y_ca(-2.0, *ca * (rdm1 * *fcaa + *fcaa * rdm1));
      sigma->ax_plus_y_vc( 8.0, (*fcvv + *favv) * *vc);
      sigma->ax_plus_y_ca(-8.0, (*fccc + *facc) * *ca);
      sigma->ax_plus_y_va( 4.0, (*fcvc + *favc) * *ca);
      sigma->ax_plus_y_ca( 4.0, (*fcvc + *favc) % *va);
      sigma->ax_plus_y_vc( 8.0, (*fcva + *fava) ^ *ca);
      sigma->ax_plus_y_ca( 4.0, *fccc * *ca * rdm1);
      sigma->ax_plus_y_ca(-4.0, *fcvc % *va * rdm1);
      sigma->ax_plus_y_va(-4.0, *fcvc * *ca * rdm1);
      sigma->ax_plus_y_vc(-2.0, *fcva * rdm1 ^ *ca);
      sigma->ax_plus_y_vc(-2.0, *va * rdm1 ^ *fcca);
      sigma->ax_plus_y_ca(-2.0, *vc % *fcva * rdm1);
      sigma->ax_plus_y_va(-2.0, *vc * *fcca * rdm1);
    }
  }
  sigma->scale(0.5);
  return sigma;
}


void CASSecond::trans_natorb() {
  auto trans = make_shared<Matrix>(nact_, nact_);
  trans->add_diag(2.0);
  blas::ax_plus_y_n(-1.0, fci_->rdm1_av()->data(), nact_*nact_, trans->data());

  VectorB occup(nact_);
  trans->diagonalize(occup);

  if (natocc_) {
    cout << " " << endl;
    cout << "  ========       state-averaged       ======== " << endl;
    cout << "  ======== natural occupation numbers ======== " << endl;
    int cnt = 0;
    for (auto& i : occup)
      cout << setprecision(4) << "   Orbital " << cnt++ << " : " << (i < 2.0 ? 2.0 - i : 0.0) << endl;
    cout << "  ============================================ " << endl;
  }

  fci_->rotate_rdms(trans);

  auto cnew = make_shared<Coeff>(*coeff_);
  cnew->copy_block(0, nclosed_, cnew->ndim(), nact_, coeff_->slice(nclosed_, nocc_) * *trans);
  coeff_ = cnew;
}
