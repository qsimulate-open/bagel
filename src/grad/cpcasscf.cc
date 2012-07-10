//
// Newint - Parallel electron correlation program.
// Filename: cpcasscf.cc
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


#include <src/grad/cpcasscf.h>
#include <cassert>
#include <src/util/f77.h>

#define CPHF_MAX_ITER 10
#define CPHF_THRESH 1.0e-8

using namespace std;

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix1e, Dvec> > grad, const shared_ptr<const Dvec> civ, 
                   const shared_ptr<const Matrix1e> eig, const shared_ptr<const DF_Half> h,
                   const shared_ptr<const DF_Half> h2, const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: grad_(grad), civector_(civ), eig_(eig), half_(h), halfjj_(h2), ref_(r), geom_(r->geom()), fci_(f) {

}


shared_ptr<PairFile<Matrix1e, Dvec> > CPCASSCF::solve() const {

  // RI determinant space
  shared_ptr<Determinants> detex(new Determinants(fci_->norb(), fci_->nelea(), fci_->neleb(), false));
  assert(fci_->norb() == ref_->nact());

  const size_t naux = geom_->naux();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = geom_->nbasis() - nocca;

  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();
  assert(nact + nclosed == nocca);

  const size_t nbasis = geom_->nbasis();

  const double* const ocoeff = ref_->coeff()->data();
  const double* const vcoeff = ocoeff + nocca*nbasis;

  // some DF vectors
  shared_ptr<const DensityFit> df = ref_->geom()->df();
  shared_ptr<const DF_Half> half = df->compute_half_transform(ocoeff, nocca)->apply_J();
  shared_ptr<const DF_Full> fullb = half->compute_second_transform(ocoeff, nocca);
  shared_ptr<const DF_Full> fulld = fullb->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);

  // making denominator...
  shared_ptr<PairFile<Matrix1e, Dvec> > denom;
  const double core_energy = ref_->geom()->nuclear_repulsion() + fci_->core_energy();
  {
    shared_ptr<Matrix1e> d0(new Matrix1e(*eig_));
    shared_ptr<Civec> d1_tmp(new Civec(*fci_->denom()));
    shared_ptr<Civec> d1_each(new Civec(d1_tmp, detex));
    shared_ptr<Dvec>  d1(new Dvec(d1_each, ref_->nstate()));
    for (int i = 0; i != ref_->nstate(); ++i)
      *d1->data(i) -= fci_->energy(i) - core_energy;
    denom = shared_ptr<PairFile<Matrix1e, Dvec> >(new PairFile<Matrix1e, Dvec>(d0, d1)); 
  }


  // CI vector
  civector_->set_det(detex);

  shared_ptr<PairFile<Matrix1e, Dvec> > source(new PairFile<Matrix1e, Dvec>(*grad_));
  // antisymmetrize
  source->first()->antisymmetrize();
  source->first()->print();
  shared_ptr<Linear<PairFile<Matrix1e, Dvec> > > solver(new Linear<PairFile<Matrix1e, Dvec> >(CPHF_MAX_ITER, source));

  // initial guess
  shared_ptr<PairFile<Matrix1e, Dvec> > z(new PairFile<Matrix1e, Dvec>(*source));

  cout << "  === CPCASSCF iteration ===" << endl << endl;

  // during the iteration, we need to
  //  (i) antisymmetrize Z, as well as
  //  (ii) project out c from z

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {

    apply_denom(z, denom); 
    z->second()->set_det(detex);
    shared_ptr<const Matrix1e> z0 = z->first();
    shared_ptr<const Dvec>     z1 = z->second();

    // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
    shared_ptr<Matrix1e> sigmaorb = compute_amat(z1, civector_);

    // computation of Atilde. Will be separated.
    // TODO index transformation can be skipped by doing so at the very end...
    shared_ptr<Matrix1e> cz0(new Matrix1e(*ref_->coeff() * *z0));

    // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
    shared_ptr<DF_Full> tmp0 = half->compute_second_transform(cz0->data(), nbasis); 
    shared_ptr<DF_Half> tmp1_1 = df->compute_half_transform(cz0->data(), nocca);
    shared_ptr<const DF_Full> tmp1 = tmp1_1->compute_second_transform(ocoeff, nbasis)->apply_J();
    tmp0->daxpy(1.0, tmp1);
    unique_ptr<double[]> term0 = tmp0->form_2index(fulld, 2.0);

    // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
    shared_ptr<DF_Full> tmp2_1 = half->compute_second_transform(cz0->data(), nocca);
    tmp2_1->symmetrize();
    shared_ptr<DF_Full> tmp2 = tmp2_1->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
    unique_ptr<double[]> term1 = half->form_2index(tmp2, 2.0);
    // mo transformation of s
    dgemm_("T", "N", nbasis, nocca, nbasis, 1.0, ocoeff, nbasis, term1.get(), nbasis, 1.0, term0.get(), nbasis); 

    // one electron part...
    shared_ptr<const Matrix1e> h(new Matrix1e(*ref_->coeff() % *ref_->hcore() * *ref_->coeff()));
    shared_ptr<const Matrix1e> htilde(new Matrix1e(*z0 % *h + *h * *z0)); 
    shared_ptr<const Matrix1e> dsa = ref_->rdm1_mat();
    dgemm_("N", "N", nbasis, nocca, nocca, 2.0, htilde->data(), nbasis, dsa->data(), nbasis, 1.0, term0.get(), nbasis); 

    daxpy_(sigmaorb->size(), 1.0, term0.get(), 1, sigmaorb->data(), 1);
    sigmaorb->antisymmetrize();

    // internal core fock operator...
    // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

    // TODO this is a reference implementation
    // first form 4 index
    unique_ptr<double[]> buf = tmp2_1->form_4index(fullb);
    // TODO Awful code. To be updated. making the code that works in the quickest possible way
    // index swap
    unique_ptr<double[]> buf2(new double[nocca*nocca*nocca*nocca]);
    for (int i = 0; i != nocca; ++i)
      for (int j = 0; j != nocca; ++j)
        for (int k = 0; k != nocca; ++k)
          for (int l = 0; l != nocca; ++l)
            buf2[l+nocca*(k+nocca*(j+nocca*i))] = buf[l+nocca*(k+nocca*(i+nocca*j))];

    // bra ket symmetrization
    for (int i = 0; i != nocca; ++i)
      for (int j = 0; j != nocca; ++j)
        for (int k = 0; k != nocca; ++k)
          for (int l = 0; l != nocca; ++l)
            buf[l+nocca*(k+nocca*(j+nocca*i))] = buf2[l+nocca*(k+nocca*(j+nocca*i))] + buf2[j+nocca*(i+nocca*(l+nocca*k))];

    unique_ptr<double[]> Htilde2(new double[nact*nact*nact*nact]);
    for (int i = nclosed; i != nocca; ++i)
      for (int j = nclosed; j != nocca; ++j)
        for (int k = nclosed; k != nocca; ++k)
          for (int l = nclosed; l != nocca; ++l)
            Htilde2[l-nclosed+nact*(k-nclosed+nact*(j-nclosed+nact*(i-nclosed)))] = buf[l+nocca*(k+nocca*(j+nocca*i))]; 

    unique_ptr<double[]> Htilde1(new double[nact*nact]);
    for (int i = nclosed; i != nocca; ++i) {
      for (int j = nclosed; j != nocca; ++j) {
        Htilde1[j-nclosed+nact*(i-nclosed)] = htilde->element(j,i);
        for (int k = 0; k != nclosed; ++k)
          Htilde1[j-nclosed+nact*(i-nclosed)] += 2.0*buf[k+nocca*(k+nocca*(j+nocca*i))] - buf[k+nocca*(i+nocca*(j+nocca*k))];
      }
    }
    // factor of 2 in the equation
    dscal_(nact*nact, 2.0, Htilde1.get(), 1);
    dscal_(nact*nact*nact*nact, 2.0, Htilde2.get(), 1);

    shared_ptr<MOFile> top(new MOFile(Htilde1, Htilde2));
    vector<int> tmp(z1->ij(), 0);
    shared_ptr<Dvec> sigmaci = fci_->form_sigma(civector_, top, tmp);

    // only here we need to have det_ instead of detex
    z1->set_det(fci_->det());
    *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);
    z1->set_det(detex);

    for (int i = 0; i != z1->ij(); ++i)
      for (int j = 0; j != z1->data(i)->size(); ++j)
        sigmaci->data(i)->data(j) -= (fci_->energy(i) - core_energy) * z1->data(i)->data(j);

sigmaci->zero();
    shared_ptr<PairFile<Matrix1e, Dvec> > sigma(new PairFile<Matrix1e, Dvec>(sigmaorb, sigmaci));

    z = solver->compute_residual(z, sigma);
z->first()->print();
    const double norm = sqrt(z->ddot(*z));
    z->scale(1.0/norm);
cout << "residual :   " << norm << endl;

  }

  return shared_ptr<PairFile<Matrix1e, Dvec> >();

}


void CPCASSCF::apply_denom(shared_ptr<PairFile<Matrix1e, Dvec> > o, shared_ptr<const PairFile<Matrix1e, Dvec> > d) const {
  shared_ptr<Matrix1e> o0 = o->first();
  shared_ptr<Dvec> o1 = o->second();
  shared_ptr<const Matrix1e> d0 = d->first();
  shared_ptr<const Dvec> d1 = d->second();

  for (size_t i = 0lu; i != o0->size(); ++i) {
    o0->data(i) /= max(0.01, d0->data(i));
  }
  for (size_t ij = 0lu; ij != o1->ij(); ++ij) {
    for (size_t k = 0lu; k != o1->lena()*o1->lenb(); ++k) {
      o1->data(ij)->data(k) /= max(0.01, d1->data(ij)->data(k));
    }
  }
}


// computes A matrix (scaled by 2 here)
shared_ptr<Matrix1e> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec, shared_ptr<const Dvec> dvec) const {
  shared_ptr<Matrix1e> amat(new Matrix1e(ref_->geom())); 

  const size_t nbasis = geom_->nbasis();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  const double* const coeff = ref_->coeff()->data(); 
  const double* const acoeff = coeff + nclosed*nbasis;

  // compute RDMs 
  shared_ptr<const RDM<1> > rdm1;
  shared_ptr<const RDM<2> > rdm2;
  tie(rdm1, rdm2) = fci_->compute_rdm12_av_from_dvec(zvec, dvec); 

  // core Fock operator
  shared_ptr<const Matrix1e> core_fock = fci_->jop()->core_fock();
  unique_ptr<double[]> buf(new double[nbasis*nact]);
  unique_ptr<double[]> buf2(new double[nbasis*nact]);
  dgemm_("N", "N", nbasis, nact, nbasis, 1.0, core_fock->data(), nbasis, acoeff, nbasis, 0.0, buf.get(), nbasis); 
  dgemm_("N", "N", nbasis, nact, nact, 1.0, buf.get(), nbasis, rdm1->data(), nact, 0.0, buf2.get(), nbasis); 
  dgemm_("T", "N", nbasis, nact, nbasis, 2.0, coeff, nbasis, buf2.get(), nbasis, 0.0, amat->element_ptr(0,nclosed), nbasis); 

  // Half transformed DF vector
  shared_ptr<const DF_Half> half = fci_->jop()->mo2e_1ext();
  shared_ptr<const DF_Full> full = half->compute_second_transform(acoeff, nact)->apply_JJ();
  shared_ptr<const DF_Full> fulld = full->apply_2rdm(rdm2->data());
  unique_ptr<double[]> jd = half->form_2index(fulld);
  dgemm_("T", "N", nbasis, nact, nbasis, 2.0, coeff, nbasis, jd.get(), nbasis, 1.0, amat->element_ptr(0,nclosed), nbasis); 

  return amat;
}
