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
#include <src/util/bfgs.h>

#define CPHF_MAX_ITER 40
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
    for (int i = 0; i != d0->size(); ++i)
      if (fabs(d0->data(i)) < 1.0e-15) d0->data(i) = 1.0e10;
    shared_ptr<Civec> d1_tmp(new Civec(*fci_->denom()));
    shared_ptr<Dvec>  d1(new Dvec(d1_tmp, ref_->nstate()));
    for (int i = 0; i != ref_->nstate(); ++i)
      *d1->data(i) -= fci_->energy(i) - core_energy;
    denom = shared_ptr<PairFile<Matrix1e, Dvec> >(new PairFile<Matrix1e, Dvec>(d0, d1)); 
  }


  // BFGS update of the denominator above
  shared_ptr<BFGS<PairFile<Matrix1e, Dvec> > > bfgs(new BFGS<PairFile<Matrix1e, Dvec> >(denom));


  // CI vector
  shared_ptr<PairFile<Matrix1e, Dvec> > source(new PairFile<Matrix1e, Dvec>(*grad_));
  // antisymmetrize
  source->first()->antisymmetrize();

  // divide by weight
  for (int ij = 0; ij != source->second()->ij(); ++ij) {
    source->second()->data(ij)->scale(1.0/fci_->weight(ij));
  }
  // project out Civector from the gradient
  source->second()->project_out(civector_);
  shared_ptr<Linear<PairFile<Matrix1e, Dvec> > > solver(new Linear<PairFile<Matrix1e, Dvec> >(CPHF_MAX_ITER, source));

  // initial guess
  shared_ptr<PairFile<Matrix1e, Dvec> > z = source->clone();
  z->zero();

  z = bfgs->extrapolate(source, z);

  // inverse matrix of C
  shared_ptr<Matrix1e> cinv(new Matrix1e(ref_->geom())); cinv->unit(); 
  {
    shared_ptr<Matrix1e> ctmp(new Matrix1e(*ref_->coeff()));
    unique_ptr<int[]> ipiv(new int[nbasis]);
    int info;
    dgesv_(nbasis, nbasis, ctmp->data(), nbasis, ipiv.get(), cinv->data(), nbasis, info);
    if (info) throw runtime_error("DGESV failed in CPCASSCF");
  }

  // State averaged density matrix
  shared_ptr<const Matrix1e> dsa = ref_->rdm1_mat();

  cout << "  === CPCASSCF iteration ===" << endl << endl;

  // during the iteration, we need to
  //  (i) antisymmetrize Z, as well as
  //  (ii) project out c from z

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
    const double norm = sqrt(z->ddot(*z));

cout << setw(4) <<  iter << " " << norm << endl;

    shared_ptr<const Matrix1e> z0 = z->first();
    shared_ptr<const Dvec>     z1 = z->second();

    // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
    // only here we need to have det_ instead of detex
#if 0
    z1->set_det(detex);
    civector_->set_det(detex);
    shared_ptr<Matrix1e> sigmaorb = compute_amat(z1, civector_);
    z1->set_det(fci_->det());
    civector_->set_det(fci_->det());
#else
shared_ptr<Matrix1e> sigmaorb = source->first()->clone(); sigmaorb->zero();
#endif

    // computation of Atilde. Will be separated.
    // TODO index transformation can be skipped by doing so at the very end...
    shared_ptr<Matrix1e> cz0(new Matrix1e(*ref_->coeff() * *z0));
    shared_ptr<Matrix1e> cz0cinv(new Matrix1e(*ref_->coeff() * *z0 * *cinv));

    // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
    // (D|jx) -> (D|jS)
    {
      shared_ptr<DF_Full> tmp0 = half->compute_second_transform(cz0cinv->data(), nbasis); 
      shared_ptr<const DF_Half> tmp1 = df->compute_half_transform(cz0->data(), nocca)->apply_J();
      tmp0->daxpy(1.0, tmp1);
      unique_ptr<double[]> buf = tmp0->form_2index(fulld, 2.0); // Factor of 2
      dgemm_("T", "N", nbasis, nocca, nbasis, 1.0, ocoeff, nbasis, buf.get(), nbasis, 1.0, sigmaorb->data(), nbasis); 
    }
    // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
    shared_ptr<DF_Full> fullz = half->compute_second_transform(cz0->data(), nocca);
    fullz->symmetrize();
    {
      shared_ptr<DF_Full> tmp = fullz->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
      unique_ptr<double[]> buf = half->form_2index(tmp, 2.0); // Factor of 2
      // mo transformation of s
      dgemm_("T", "N", nbasis, nocca, nbasis, 1.0, ocoeff, nbasis, buf.get(), nbasis, 1.0, sigmaorb->data(), nbasis); 
    }

    // one electron part...
    shared_ptr<Matrix1e> htilde(new Matrix1e(*cz0 % *ref_->hcore() * *ref_->coeff())); 
    htilde->symmetrize();
    dgemm_("N", "N", nbasis, nocca, nocca, 4.0, htilde->data(), nbasis, dsa->data(), nbasis, 1.0, sigmaorb->data(), nbasis); 

    sigmaorb->antisymmetrize();
    sigmaorb->purify_redrotation(nclosed, nact, nvirt);

    // internal core fock operator...
    // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

#if 0
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

    shared_ptr<MOFile> top(new Htilde(ref_, nclosed, nocca, move(Htilde1), move(Htilde2)));
    vector<int> tmp(z1->ij(), 0);
    shared_ptr<Dvec> sigmaci = fci_->form_sigma(civector_, top, tmp);

    *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);

    for (int i = 0; i != z1->ij(); ++i)
      for (int j = 0; j != z1->data(i)->size(); ++j)
        sigmaci->data(i)->data(j) -= (fci_->energy(i) - core_energy) * z1->data(i)->data(j);

    sigmaci->project_out(civector_);

sigmaci->zero();
#else
shared_ptr<Dvec> sigmaci = source->second()->clone(); sigmaci->zero();
#endif

    shared_ptr<PairFile<Matrix1e, Dvec> > sigma(new PairFile<Matrix1e, Dvec>(sigmaorb, sigmaci));

    z = solver->compute_residual(z, sigma);
    z = bfgs->extrapolate(z, solver->civec());

    z->first()->purify_redrotation(nclosed, nact, nvirt);
    z->second()->project_out(civector_);

z->second()->zero();
if (sqrt(z->ddot(*z)) < 1.0e-8) break;

  }

solver->civec()->first()->print("result", 12);
solver->civec()->second()->print(0);

  return solver->civec();

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
