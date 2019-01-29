//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cpcasscf.cc
// Copyright (C) 2012 Toru Shiozaki
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

#include <src/grad/cpcasscf.h>
#include <src/util/math/linearRM.h>
#include <src/scf/hf/fock.h>

using namespace std;
using namespace bagel;
using namespace btas;

CPCASSCF::CPCASSCF(shared_ptr<const PairFile<Matrix, Dvec>> grad, shared_ptr<const Dvec> civ, shared_ptr<const DFHalfDist> h,
                   shared_ptr<const Reference> r, shared_ptr<FCI_base> f, const int ncore, const bool imag, shared_ptr<const Matrix> coeff)
: grad_(grad), civector_(civ), halfj_(h), ref_(r), geom_(r->geom()), fci_native_(f), ncore_(ncore), imag_(imag), coeff_(coeff ? coeff : ref_->coeff()) {

  // FCI object in CPCASSCF should be Knowles--Handy (due to form_sigma)
  if (ref_->nact())
    fci_ = make_shared<KnowlesHandy>(fci_native_->conv_to_ciwfn(), r);

  if (ref_->nact() && coeff_ != ref_->coeff())
    fci_->update(coeff_);

  qvec_ = ref_->nact() ? make_shared<Qvec>(coeff_->mdim(), ref_->nact(), coeff_, ref_->nclosed(), fci_, ref_->rdm2_av()) : nullptr;
}


tuple<shared_ptr<Matrix>,shared_ptr<Matrix>, shared_ptr<Matrix>> CPCASSCF::compute_orb_denom_and_fock() const {
  // denominator is calculated as in Theor Chem Acc (1997) 97:88-95
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nocc = ref_->nocc();
  const int nvirt = ref_->nvirt();
  const int nmobasis = coeff_->mdim();

  auto denom = make_shared<Matrix>(nmobasis, nmobasis);
  denom->fill(1.0e20);

  vector<double> occup;
  if (nact) occup = ref_->rdm1_av()->diag();

  // inactive fock in MO basis (when nact=0 this code is not efficient, but this case is for debugging)
  shared_ptr<const Matrix> core_fock = nact ? fci_->jop()->core_fock()
                                            : make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, coeff_->slice(0,nclosed), /*grad*/false, /*rhf*/true);
  auto finact = make_shared<Matrix>(*coeff_ % *core_fock * *coeff_);

  shared_ptr<Matrix> fock;
  if (nact) {
    auto rdm1av = make_shared<Matrix>(nact, nact);
    copy_n(ref_->rdm1_av()->data(), rdm1av->size(), rdm1av->data());
    rdm1av->sqrt();
    rdm1av->scale(1.0/sqrt(2.0));

    const MatView acoeff = coeff_->slice(nclosed, nclosed+nact);
    auto fact_ao = make_shared<Fock<1>>(geom_, fci_->jop()->core_fock()->clone(), nullptr, acoeff * *rdm1av, /*grad*/false, /*rhf*/true);
    fock = make_shared<Matrix>(*finact + *coeff_ % *fact_ao * *coeff_);
  } else {
    fock = finact;
  }

  shared_ptr<Matrix> fact = nact ? qvec_->copy() : nullptr;
  for (int i = 0; i != nact; ++i)
    blas::ax_plus_y_n(occup[i], finact->element_ptr(0,nclosed+i), nmobasis, fact->data()+i*nmobasis);

  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nvirt; ++j)
      denom->element(j+nocc,i+nclosed) = denom->element(i+nclosed,j+nocc) = -2.0*fact->element(i+nclosed,i) + 2.0*occup[i]*fock->element(j+nocc, j+nocc);

  for (int i = 0; i != nclosed; ++i)
    for (int j = 0; j != nvirt; ++j)
      denom->element(j+nocc,i) = denom->element(i,j+nocc) = 4.0*fock->element(j+nocc, j+nocc) - 4.0*fock->element(i, i);

  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nclosed; ++j)
      denom->element(j,i+nclosed) = denom->element(i+nclosed,j)
                                  = 4.0*fock->element(nclosed+i,nclosed+i) - 2.0*fact->element(i+nclosed,i) - fock->element(j, j)*(4.0 - 2.0*occup[i]);

  return make_tuple(denom, fock, finact);
}


tuple<shared_ptr<const Matrix>, shared_ptr<const Dvec>, shared_ptr<const Matrix>, shared_ptr<const Matrix>>
  CPCASSCF::solve(const double zthresh, const int zmaxiter, shared_ptr<const Matrix> additional_den, const bool project_all) {

  // helper function
  auto project_out = [&project_all](shared_ptr<Dvec> a, shared_ptr<const Dvec> b) {
    if (project_all)
      a->project_out_all(b);
    else
      a->project_out(b);
  };

  const size_t nocca = ref_->nocc();
  const int nmobasis = coeff_->mdim();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  const int nvirt = ref_->nvirt();
  assert(nact + nclosed == nocca);

  // RI determinant space
  auto detex = nact ? make_shared<Determinants>(fci_->norb(), fci_->nelea(), fci_->neleb(), false, /*mute=*/true)
                    : make_shared<Determinants>();
  assert(!nact || fci_->norb() == nact);

  const MatView ocoeff = coeff_->slice(0, nocca);

  // some DF vectors
  shared_ptr<const DFFullDist> fullb = halfj_->compute_second_transform(ocoeff);

  // making denominator...
  shared_ptr<PairFile<Matrix, Dvec>> denom;
  {
    shared_ptr<Matrix> d0;
    // setting fock_ and fockinact_ here
    tie(d0, fock_, fockinact_) = compute_orb_denom_and_fock();

    shared_ptr<Dvec> d1;
    if (nact) {
      d1 = make_shared<Dvec>(fci_->denom(), ref_->nstate());
      const double core_energy = geom_->nuclear_repulsion() + fci_->core_energy();
      for (int i = 0; i != ref_->nstate(); ++i)
        *d1->data(i) -= fci_->energy(i) - core_energy;
    } else {
      // creating Dvec with one element.
      d1 = make_shared<Dvec>(detex, ref_->nstate());
      d1->data(0)->element(0,0) = 1.0e10;
    }
    denom = make_shared<PairFile<Matrix, Dvec>>(d0, d1);
  }

  // frozen core contributions
  shared_ptr<Matrix> zcore, gzcore;
  if (ncore_ || imag_ || additional_den) {
    assert(ncore_ < nclosed);
    zcore = make_shared<Matrix>(nmobasis, nmobasis);

    if (imag_) {
      for (int i = 0; i != nclosed; ++i) {
        for (int j = 0; j != nclosed; ++j) {
          if (i == j) continue;
          const double fdiff = fock_->element(j,j) - fock_->element(i,i);
          zcore->element(j, i) = fabs(fdiff) < 1.0e-8 ? 0.0 : - .5 * (grad_->first()->element(j, i) - grad_->first()->element(i, j)) / fdiff;
        }
      }
      for (int i = nocca; i != nocca + nvirt; ++i) {
        for (int j = nocca; j != nocca + nvirt; ++j) {
          if (i == j) continue;
          const double fdiff = fock_->element(j,j) - fock_->element(i,i);
          zcore->element(j, i) = fabs(fdiff) < 1.0e-8 ? 0.0 : - .5 * (grad_->first()->element(j, i) - grad_->first()->element(i, j)) / fdiff;
        }
      }
    } else {
      for (int i = 0; i != ncore_; ++i)
        for (int j = ncore_; j != nclosed; ++j) {
          zcore->element(j, i) = - (grad_->first()->element(j, i) - grad_->first()->element(i, j)) / (fock_->element(j,j) - fock_->element(i,i));
          assert(abs(fock_->element(i, j)) < 1.0e-8);
        }
    }

    zcore->symmetrize();
    if (additional_den) {
      zcore->add_block(1.0, 0, 0, nocca, nocca, *additional_den);
    }

    shared_ptr<Matrix> rot;

    // find good lambda that makes zcore positive definite
    double lambda = 1.0;

    for (int i = 0; i != zmaxiter; ++i) {
      shared_ptr<Matrix> rotq;
      VectorB eig(nmobasis);
      if (nact) {
        rotq = make_shared<Matrix>(*zcore * lambda);
        for (int i = 0; i != nmobasis; ++i)
          rotq->element(i,i) += 2.0;
      } else {
        rotq = make_shared<Matrix>(*zcore * lambda);
        for (int i = 0; i != nmobasis; ++i)
          rotq->element(i,i) += 2.0;
      }
      rotq->diagonalize(eig);
      if (eig[0] < -numerical_zero__) lambda *= 0.5;
      else break;
    }

    if (nact) {
      rot = make_shared<Matrix>(*zcore * lambda);
      for (int i = 0; i != nmobasis; ++i) {
        rot->element(i,i) += 2.0;
      }
    } else {
      rot = make_shared<Matrix>(*zcore * lambda);
      for (int i = 0; i != nmobasis; ++i) {
        rot->element(i,i) += 2.0;
      }
    }
    rot->sqrt();
    rot->scale(1.0/sqrt(2.0));
    auto gzcoreao = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, *coeff_ * *rot, false, true);
    auto rottemp = make_shared<Matrix>(nmobasis, nmobasis);
    for (int i = 0; i != nmobasis; ++i) rottemp->element(i,i) = 1.0;
    auto fv = make_shared<Fock<1>>(geom_, ref_->hcore(), nullptr, *coeff_ * *rottemp, false, true);
    gzcore = make_shared<Matrix>(*coeff_ % (*gzcoreao - *fv) * *coeff_); // compensate for the "trick"
    gzcore->scale(1.0 / lambda);
  }

  // gradient Y and y
  auto source = make_shared<PairFile<Matrix, Dvec>>(*grad_);
  // divide by weight
  if (source->second()->ij() > 1) {
    for (int ij = 0; ij != source->second()->ij(); ++ij)
      source->second()->data(ij)->scale(1.0/fci_->weight(ij));
  }
  // patch frozen core contributions
  if (zcore) {
    // contributions to Y
    source->first()->ax_plus_y(2.0, *fock_ * *zcore + *gzcore * *ref_->rdm1_mat()->resize(nmobasis, nmobasis));
    // contributions to y
    if (nact) {
      for (int istate = 0; istate != ref_->nstate(); ++istate) {
        shared_ptr<const Dvec> rdm1deriv = fci_->rdm1deriv(istate);
        for (int i = 0; i != nact; ++i)
          for (int j = 0; j != nact; ++j)
            source->second()->data(istate)->ax_plus_y(gzcore->element(j+nclosed, i+nclosed), rdm1deriv->data(j+nact*i));
      }
    }
  }
  // antisymmetrize
  source->first()->antisymmetrize();

  // set zero to redundant rotation
  auto purify = [](Matrix& mat, const int nc, const int na, const int nv) {
    for (int g = 0; g != nc; ++g)
      for (int h = 0; h != nc; ++h)
        mat(h, g) = 0.0;
    for (int g = 0; g != na; ++g)
      for (int h = 0; h != na; ++h)
        mat(h+nc, g+nc) = 0.0;
    for (int g = 0; g != nv; ++g)
      for (int h = 0; h != nv; ++h)
        mat(h+nc+na, g+nc+na) = 0.0;
  };
  purify(*source->first(), nclosed, nact, nvirt);

  // project out Civector from the gradient
  project_out(source->second(), civector_);

  const double lambda = max(sqrt(static_cast<double>(source->second()->size()) / source->first()->size()), 1.0);
  source->second()->scale(1.0/lambda);
  denom->second()->scale(1.0/(lambda*lambda));

  auto solver = make_shared<LinearRM<PairFile<Matrix, Dvec>>>(zmaxiter, source);

  // initial guess
  shared_ptr<PairFile<Matrix, Dvec>> z = source->copy();
  *z /= *denom;
  project_out(z->second(), civector_);
  z->scale(1.0/z->norm());

  // inverse matrix of C
  auto ovl = make_shared<const Overlap>(geom_);
  auto cinv = make_shared<const Matrix>(*coeff_ % *ovl);

  cout << "  === CASSCF Z-vector iteration ===" << endl << endl;

  Timer timer;
  for (int iter = 0; iter != zmaxiter; ++iter) {
    // given z, computes sigma (before anti-symmetrization)
    shared_ptr<PairFile<Matrix, Dvec>> sigma = form_sigma(z, fullb, detex, cinv, /*antisym*/true, lambda);
    sigma->first()->antisymmetrize();
    purify(*sigma->first(), nclosed, nact, nvirt);
    project_out(sigma->second(), civector_);

    z = solver->compute_residual(z, sigma);

    cout << setw(10) <<  iter << " " << setw(17) << setprecision(10) << z->first()->rms()
                                     << setw(17) << z->second()->rms()*lambda  << setw(10) << setprecision(2) << timer.tick() << endl;
    if (z->first()->rms()+z->second()->rms()*lambda < zthresh) break;

    *z /= *denom;
    project_out(z->second(), civector_);
    z->scale(1.0/z->norm());
  }

  shared_ptr<PairFile<Matrix, Dvec>> result = solver->civec();
  result->second()->scale(1.0/lambda);
  shared_ptr<Matrix> xmat = form_sigma(result, fullb, detex, cinv, /*antisym*/false, 1.0)->first();

  *xmat += *grad_->first();
  if (zcore)
    xmat->ax_plus_y(2.0, *fock_ * *zcore + *gzcore * *ref_->rdm1_mat()->resize(nmobasis, nmobasis));
  xmat->symmetrize();
  xmat->scale(0.5); // due to convention
  return make_tuple(result->first(), result->second(), xmat, zcore);
}


shared_ptr<PairFile<Matrix,Dvec>>
  CPCASSCF::form_sigma(shared_ptr<const PairFile<Matrix,Dvec>> z,
                       shared_ptr<const DFFullDist> fullb, shared_ptr<const Determinants> detex, shared_ptr<const Matrix> cinv,
                       const bool antisym, const double lambda) const {

  const size_t nmobasis = coeff_->mdim();
  const size_t nocca = ref_->nocc();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  shared_ptr<const Matrix> z0 = z->first();
  shared_ptr<const Dvec>   z1 = z->second();
  shared_ptr<RDM<1>> rdm1_av = nact ? ref_->rdm1_av()->copy() : nullptr;
  shared_ptr<RDM<2>> rdm2_av = nact ? ref_->rdm2_av()->copy() : nullptr;

  // orbital part that depends on ci-z
  shared_ptr<Matrix> sigmaorb = nact ? compute_amat(z1, civector_, detex) : make_shared<Matrix>(nmobasis,nmobasis);
  sigmaorb->scale(1.0/lambda);

  // computation of Atilde.
  auto cz0 = make_shared<Matrix>(*coeff_ * *z0);
  auto cz0cinv = make_shared<Matrix>(*coeff_ * *z0 * *cinv);
  const MatView ocz0 = cz0->slice(0, nocca);

  // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
  // (D|jx) -> (D|jS)
  {
    shared_ptr<DFFullDist> tmp0 = halfj_->compute_second_transform(cz0cinv);
    shared_ptr<const DFHalfDist> tmp1 = geom_->df()->compute_half_transform(ocz0)->apply_J();
    tmp0->ax_plus_y(1.0, tmp1);
    shared_ptr<const DFFullDist> fulld = nact ? fullb->apply_2rdm(*rdm2_av, *rdm1_av, nclosed, nact)
                                              : fullb->apply_closed_2RDM();
    shared_ptr<const Matrix> buf = tmp0->form_2index(fulld, 2.0); // Factor of 2
    const Matrix cbuf(*coeff_ % *buf);
    sigmaorb->add_block(1.0, 0, 0, nmobasis, nocca, cbuf);
  }
  // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
  shared_ptr<DFFullDist> fullz = halfj_->compute_second_transform(ocz0);
  fullz->symmetrize();
  {
    shared_ptr<const DFFullDist> tmp = nact ? fullz->apply_2rdm(*rdm2_av, *rdm1_av, nclosed, nact)
                                            : fullz->apply_closed_2RDM();
    shared_ptr<const Matrix> buf = halfj_->form_2index(tmp, 2.0); // Factor of 2
    // mo transformation of s
    sigmaorb->add_block(1.0, 0, 0, nmobasis, nocca, *coeff_ % *buf);
  }

  // one electron part...
  auto htilde = make_shared<Matrix>(*coeff_ % *ref_->hcore() * *cz0);
  htilde->symmetrize();
  htilde->scale(2.0);
  const Matrix cbuf(htilde->slice(0, nocca) * *ref_->rdm1_mat());
  sigmaorb->add_block(2.0, 0, 0, nmobasis, nocca, cbuf);

  shared_ptr<Dvec> sigmaci;
  if (antisym) {
    // At this point
    // htilde = Z^daggerh + hZ
    // fullb  = (D|ij)
    // fullz  = (D|ir)Z_rj + (D|rj)Z_ri

    // internal core fock operator...
    // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

    // first form 4 index
    shared_ptr<Matrix> buf = fullz->form_4index(fullb, 1.0);
    buf->symmetrize();
    buf->scale(2.0);

    auto Htilde2 = make_shared<Matrix>(nact*nact, nact*nact);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii)
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj)
        for (int k = nclosed, kk = 0; k != nocca; ++k, ++kk)
          copy_n(buf->element_ptr(nclosed+nocca*k,j+nocca*i), nocca-nclosed, Htilde2->element_ptr(nact*kk, jj+nact*ii));

    auto Htilde1 = make_shared<Matrix>(nact,nact, true);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii) {
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj) {
        (*Htilde1)(jj, ii) = htilde->element(j,i);
        for (int k = 0; k != nclosed; ++k)
          (*Htilde1)(jj, ii) += 2.0*buf->element(k+nocca*k, j+nocca*i) - buf->element(k+nocca*i, j+nocca*k);
      }
    }
    // factor of 2 in the equation
    *Htilde1 *= 2.0;
    *Htilde2 *= 2.0;

    auto top = make_shared<Htilde>(ref_, 0, nact, Htilde1, Htilde2);
    vector<int> tmp(z1->ij(), 0);

    if (nact) {
      sigmaci = fci_->form_sigma(civector_, top, tmp);
      sigmaci->scale(lambda);
      *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);

      const double core_energy = geom_->nuclear_repulsion() + fci_->core_energy();
      for (int i = 0; i != z1->ij(); ++i)
        blas::ax_plus_y_n(-(fci_->energy(i)-core_energy), z1->data(i)->data(), z1->data(i)->size(), sigmaci->data(i)->data());
      sigmaci->scale(1.0/(lambda*lambda));
    } else {
      sigmaci = civector_->clone(); // always zero when nact=0
    }

  } else {
    // TODO I need to understand... (this is what Molpro does, and it works).
    Matrix qone(nocca, nocca);
    if (nclosed)
      qone.add_block(2.0, 0, 0, nocca, nclosed, fock_->get_submatrix(0, 0, nocca, nclosed));

    if (nact) {
      qone.add_block(1.0, 0, nclosed, nocca, nact, qvec_->get_submatrix(0, 0, nocca, nact));
      Matrix rdm1mat(nact, nact);
      copy_n(rdm1_av->data(), rdm1_av->size(), rdm1mat.data());
      qone.add_block(1.0, 0, nclosed, nocca, nact, *fockinact_->get_submatrix(0,nclosed,nocca,nact) * rdm1mat);
    }
    qone.symmetrize();
    const MatView oz0 = z0->slice(0, nocca);
    sigmaorb->add_block(4.0, 0, 0, nmobasis, nocca, oz0 * qone);
  }

  return make_shared<PairFile<Matrix, Dvec>>(sigmaorb, sigmaci);
}



// computes A matrix (scaled by 2 here)
shared_ptr<Matrix> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec, shared_ptr<const Dvec> dvec, shared_ptr<const Determinants> o) const {

  const size_t nmobasis = coeff_->mdim();

  auto amat = make_shared<Matrix>(nmobasis, nmobasis);
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  const MatView acoeff = coeff_->slice(nclosed, nclosed+nact);

  // compute RDMs
  shared_ptr<const RDM<1>> rdm1t;
  shared_ptr<const RDM<2>> rdm2t;
  tie(rdm1t, rdm2t) = fci_->compute_rdm12_av_from_dvec(zvec, dvec, o);

  // symmetrize
  auto rdm1mat = make_shared<Matrix>(nact, nact);
  copy_n(rdm1t->data(), rdm1t->size(), rdm1mat->data());
  rdm1mat->symmetrize();

  shared_ptr<RDM<2>> rdm2 = rdm2t->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j) {
      blas::ax_plus_y_n(0.5, rdm2t->element_ptr(0,0,j,i), nact*nact, rdm2->element_ptr(0,0,j,i));
      blas::ax_plus_y_n(0.5, rdm2t->element_ptr(0,0,i,j), nact*nact, rdm2->element_ptr(0,0,j,i));
    }

  // core Fock operator
  shared_ptr<const Matrix> core_fock = fci_->jop()->core_fock();
  amat->add_block(2.0, 0, nclosed, nmobasis, nact, *coeff_ % (*core_fock * acoeff * *rdm1mat));

  // Half transformed DF vector
  shared_ptr<const DFHalfDist> half = fci_->jop()->mo2e_1ext();
  shared_ptr<const DFFullDist> full = half->apply_JJ()->compute_second_transform(acoeff);
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(*rdm2);
  amat->add_block(2.0, 0, nclosed, nmobasis, nact, *coeff_ % *half->form_2index(fulld, 1.0));

  // additing f^z_ri contribution
  if (nclosed) {
    shared_ptr<const Matrix> aden = make_shared<Matrix>(acoeff * *rdm1mat ^ acoeff);
    shared_ptr<const Matrix> adenj = make_shared<Matrix>(*rdm1mat ^ acoeff);

    const MatView ocoeff = coeff_->slice(0, nclosed);
    // coulomb
    Matrix fockz(*geom_->df()->compute_Jop(half, adenj, /*only once*/false) * ocoeff);
    // exchange
    shared_ptr<DFFullDist> halfd = half->compute_second_transform(ocoeff);
    halfd = halfd->transform_occ1(rdm1mat);
    fockz += *half->form_2index(halfd->apply_JJ(), -0.5);
    // add to amat
    amat->add_block(4.0, 0, 0, nmobasis, nclosed, *coeff_ % fockz);
  }

  return amat;
}
