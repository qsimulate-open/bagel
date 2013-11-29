//
// BAGEL - Parallel electron correlation program.
// Filename: cpcasscf.cc
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


#include <src/grad/cpcasscf.h>
#include <src/math/linearRM.h>
#include <src/math/bfgs.h>

#define CPHF_MAX_ITER 100
#define CPHF_THRESH 1.0e-10

using namespace std;
using namespace bagel;

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix, Dvec>> grad, const shared_ptr<const Dvec> civ,
                   const shared_ptr<const Matrix> eig, const shared_ptr<const DFHalfDist> h,
                   const shared_ptr<const DFHalfDist> h2, const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: grad_(grad), civector_(civ), eig_(eig), half_(h), halfjj_(h2), ref_(r), geom_(r->geom()), fci_(f) {

  cout << "   CI vectors:" << endl;
  civector_->print(-1);

}


shared_ptr<PairFile<Matrix, Dvec>> CPCASSCF::solve() const {

  // RI determinant space
  auto detex = make_shared<Determinants>(fci_->norb(), fci_->nelea(), fci_->neleb(), false, /*mute=*/true);
  assert(fci_->norb() == ref_->nact());

  const size_t nmobasis = ref_->coeff()->mdim();
  const size_t naobasis = geom_->nbasis();
  const size_t naux = geom_->naux();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = nmobasis - nocca;

  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  assert(nact + nclosed == nocca);

  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(0, nocca);

  // some DF vectors
  shared_ptr<const DFDist> df = geom_->df();
  shared_ptr<const DFHalfDist> half = df->compute_half_transform(ocoeff)->apply_J();
  shared_ptr<const DFFullDist> fullb = half->compute_second_transform(ocoeff);

  // making denominator...
  shared_ptr<PairFile<Matrix, Dvec>> denom;
  const double core_energy = geom_->nuclear_repulsion() + fci_->core_energy();
  {
    shared_ptr<Matrix> d0 = eig_->copy();
    for (auto& i : *d0)
      if (fabs(i) < 1.0e-10) i = 1.0;
    auto d1_tmp = make_shared<const Civec>(*fci_->denom());
    auto d1 = make_shared<Dvec>(d1_tmp, ref_->nstate());
    for (int i = 0; i != ref_->nstate(); ++i)
      *d1->data(i) -= fci_->energy(i) - core_energy;
#if 1
    // TODO understand this factor of 2
    *d1 *= 2;
#endif
    denom = make_shared<PairFile<Matrix, Dvec>>(d0, d1);
  }

  // BFGS update of the denominator above
#if 0
  shared_ptr<BFGS<PairFile<Matrix, Dvec>>> bfgs(new BFGS<PairFile<Matrix, Dvec>>(denom, false));
#else
  auto bfgs = make_shared<BFGS<PairFile<Matrix, Dvec>>>(denom);
#endif


  // CI vector
  auto source = make_shared<PairFile<Matrix, Dvec>>(*grad_);
  // antisymmetrize
  source->first()->antisymmetrize();

  // divide by weight
  for (int ij = 0; ij != source->second()->ij(); ++ij) {
    source->second()->data(ij)->scale(1.0/fci_->weight(ij));
  }
  // project out Civector from the gradient
  source->second()->project_out(civector_);
  auto solver = make_shared<LinearRM<PairFile<Matrix, Dvec>>>(CPHF_MAX_ITER, source);

  // initial guess
  shared_ptr<PairFile<Matrix, Dvec>> z = source->clone();
  z->zero();

  z = bfgs->extrapolate(source, z);
// not needed as z->xecond() is zero
//z->second()->project_out(civector_);

  // inverse matrix of C
//shared_ptr<Matrix> cinv(new Matrix(*ref_->coeff())); cinv->inverse();
  auto ovl = make_shared<const Overlap>(geom_);
  auto cinv = make_shared<const Matrix>(*ref_->coeff() % *ovl);

  // State averaged density matrix
  shared_ptr<const Matrix> dsa = ref_->rdm1_mat()->resize(nmobasis, nmobasis); // TODO resize is just a waste of time..

  cout << "  === CPCASSCF iteration ===" << endl << endl;

  // during the iteration, we need to
  //  (i) antisymmetrize Z, as well as
  //  (ii) project out c from z

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
    const double norm = sqrt(z->dot_product(*z));
    cout << setw(4) <<  iter << " " << setprecision(14) << norm << endl;

    shared_ptr<const Matrix> z0 = z->first();
    shared_ptr<const Dvec>     z1 = z->second();

    // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
    // only here we need to have det_ instead of detex
    shared_ptr<Matrix> sigmaorb = compute_amat(z1, civector_, detex);

    // computation of Atilde. Will be separated.
    // TODO index transformation can be skipped by doing so at the very end...
    auto cz0 = make_shared<Matrix>(*ref_->coeff() * *z0);
    auto cz0cinv = make_shared<Matrix>(*ref_->coeff() * *z0 * *cinv);

    // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
    // (D|jx) -> (D|jS)
    {
      shared_ptr<DFFullDist> tmp0 = half->compute_second_transform(cz0cinv);
      shared_ptr<const DFHalfDist> tmp1 = df->compute_half_transform(cz0->slice(0,nocca))->apply_J();
      tmp0->ax_plus_y(1.0, tmp1);
      shared_ptr<const DFFullDist> fulld = fullb->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
      shared_ptr<const Matrix> buf = tmp0->form_2index(fulld, 2.0); // Factor of 2
      dgemm_("T", "N", nmobasis, nocca, naobasis, 1.0, ocoeff->data(), naobasis, buf->data(), naobasis, 1.0, sigmaorb->data(), nmobasis);
    }
    // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
    shared_ptr<DFFullDist> fullz = half->compute_second_transform(cz0->slice(0,nocca));
    fullz->symmetrize();
    {
      shared_ptr<const DFFullDist> tmp = fullz->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
      shared_ptr<const Matrix> buf = half->form_2index(tmp, 2.0); // Factor of 2
      // mo transformation of s
      dgemm_("T", "N", nmobasis, nocca, naobasis, 1.0, ocoeff->data(), naobasis, buf->data(), naobasis, 1.0, sigmaorb->data(), nmobasis);
    }

    // one electron part...
    auto htilde = make_shared<Matrix>(*cz0 % *ref_->hcore() * *ref_->coeff());
    htilde->symmetrize();
    *htilde *= 2.0;
    dgemm_("N", "N", nmobasis, nocca, nocca, 2.0, htilde->data(), nmobasis, dsa->data(), nmobasis, 1.0, sigmaorb->data(), nmobasis);

    sigmaorb->antisymmetrize();

    // At this point
    // htilde = Zh + hZ^dagger
    // fullb  = (D|ij)
    // fullz  = (D|ir)Z_rj + (D|rj)Z_ri

    // internal core fock operator...
    // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

    // TODO this is a reference implementation
    // first form 4 index
    shared_ptr<Matrix> buf = fullz->form_4index(fullb, 1.0);
    // TODO Awful code. To be updated. making the code that works in the quickest possible way
    // index swap
    unique_ptr<double[]> buf2(new double[nocca*nocca*nocca*nocca]);

    // bra ket symmetrization
    for (int i = 0; i != nocca*nocca; ++i)
      for (int j = 0; j != nocca*nocca; ++j)
        buf2[j+nocca*nocca*i] = buf->element(j, i) + buf->element(i, j);

    auto Htilde2 = make_shared<Matrix>(nact*nact, nact*nact);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii)
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj)
        for (int k = nclosed, kk = 0; k != nocca; ++k, ++kk)
          for (int l = nclosed, ll = 0; l != nocca; ++l, ++ll)
            Htilde2->element(ll+nact*kk, jj+nact*ii) = buf2[l+nocca*(k+nocca*(j+nocca*i))];

    auto Htilde1 = make_shared<Matrix>(nact,nact, true);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii) {
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj) {
        (*Htilde1)(jj, ii) = htilde->element(j,i);
        for (int k = 0; k != nclosed; ++k)
          (*Htilde1)(jj, ii) += 2.0*buf2[k+nocca*(k+nocca*(j+nocca*i))] - buf2[k+nocca*(i+nocca*(j+nocca*k))];
      }
    }
    // factor of 2 in the equation
    *Htilde1 *= 2.0;
    *Htilde2 *= 2.0;

    auto top = make_shared<Htilde>(ref_, nclosed, nocca, Htilde1, Htilde2);
    vector<int> tmp(z1->ij(), 0);
    shared_ptr<Dvec> sigmaci = fci_->form_sigma(civector_, top, tmp);

    *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);

    for (int i = 0; i != z1->ij(); ++i)
      for (int j = 0; j != z1->data(i)->size(); ++j)
        sigmaci->data(i)->data(j) -= (fci_->energy(i) - core_energy) * z1->data(i)->data(j);

    sigmaci->project_out(civector_);

    auto sigma = make_shared<PairFile<Matrix, Dvec>>(sigmaorb, sigmaci);

    z = solver->compute_residual(z, sigma);

    z = bfgs->extrapolate(z, solver->civec());
    z->second()->project_out(civector_);

    if (sqrt(z->dot_product(*z)) < CPHF_THRESH) break;

  }

solver->civec()->second()->print(-1);
  return solver->civec();

}


// computes A matrix (scaled by 2 here)
shared_ptr<Matrix> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec, shared_ptr<const Dvec> dvec, shared_ptr<const Determinants> o) const {

  const size_t naobasis = geom_->nbasis();
  const size_t nmobasis = ref_->coeff()->mdim();

  auto amat = make_shared<Matrix>(nmobasis, nmobasis);
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  const double* const coeff = ref_->coeff()->data();
  shared_ptr<const Matrix> acoeff = ref_->coeff()->slice(nclosed, nclosed+nact);

  // compute RDMs
  shared_ptr<const RDM<1>> rdm1t;
  shared_ptr<const RDM<2>> rdm2t;
  tie(rdm1t, rdm2t) = fci_->compute_rdm12_av_from_dvec(zvec, dvec, o);

  // symmetrize
  shared_ptr<RDM<1>> rdm1 = rdm1t->clone();
  shared_ptr<RDM<2>> rdm2 = rdm2t->clone();
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      rdm1->element(j,i) = 0.5*(rdm1t->element(j,i)+rdm1t->element(i,j));
  for (int i = 0; i != nact; ++i)
    for (int j = 0; j != nact; ++j)
      for (int k = 0; k != nact; ++k)
        for (int l = 0; l != nact; ++l)
//        rdm2->element(l,k,j,i) = 0.125*(rdm2t->element(l,k,j,i)+rdm2t->element(k,l,j,i)+rdm2t->element(l,k,i,j)+rdm2t->element(k,l,i,j)
//                                      + rdm2t->element(j,i,l,k)+rdm2t->element(j,i,k,l)+rdm2t->element(i,j,l,k)+rdm2t->element(i,j,k,l));
//        rdm2->element(l,k,j,i) = 0.25*(rdm2t->element(l,k,j,i)+rdm2t->element(k,l,j,i)+rdm2t->element(l,k,i,j)+rdm2t->element(k,l,i,j));
          rdm2->element(l,k,j,i) = 0.5*(rdm2t->element(l,k,j,i)+rdm2t->element(k,l,j,i));

  // prefactor
  const double prefactor = 2.0;

  // core Fock operator
  shared_ptr<const Matrix> core_fock = fci_->jop()->core_fock();
  unique_ptr<double[]> buf(new double[naobasis*nact]);
  unique_ptr<double[]> buf2(new double[naobasis*nact]);
  dgemm_("N", "N", naobasis, nact, naobasis, 1.0, core_fock->data(), naobasis, acoeff->data(), naobasis, 0.0, buf.get(), naobasis);
  dgemm_("N", "N", naobasis, nact, nact, 1.0, buf.get(), naobasis, rdm1->data(), nact, 0.0, buf2.get(), naobasis);
  dgemm_("T", "N", nmobasis, nact, naobasis, prefactor, coeff, naobasis, buf2.get(), naobasis, 0.0, amat->element_ptr(0,nclosed), nmobasis);

  // Half transformed DF vector
#if 0
  shared_ptr<const DFHalfDist> half = fci_->jop()->mo2e_1ext();
#else
  shared_ptr<const DFHalfDist> half = geom_->df()->compute_half_transform(acoeff);
#endif
  shared_ptr<const DFFullDist> full = half->compute_second_transform(acoeff)->apply_JJ();
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(rdm2->data());
  shared_ptr<const Matrix> jd = half->form_2index(fulld, 1.0);
  dgemm_("T", "N", nmobasis, nact, naobasis, prefactor, coeff, naobasis, jd->data(), naobasis, 1.0, amat->element_ptr(0,nclosed), nmobasis);

  return amat;
}
