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
// The BAGEL package is free software; you can redistribute it and\/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
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


#include <stddef.h>
#include <src/grad/cpcasscf.h>
#include <src/util/linearRM.h>
#include <cassert>
#include <src/util/f77.h>
#include <src/util/bfgs.h>

#define CPHF_MAX_ITER 100
#define CPHF_THRESH 1.0e-10

using namespace std;
using namespace bagel;

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix, Dvec> > grad, const shared_ptr<const Dvec> civ,
                   const shared_ptr<const Matrix> eig, const shared_ptr<const DFHalfDist> h,
                   const shared_ptr<const DFHalfDist> h2, const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: grad_(grad), civector_(civ), eig_(eig), half_(h), halfjj_(h2), ref_(r), geom_(r->geom()), fci_(f) {

  cout << "   CI vectors:" << endl;
  civector_->print(-1);

}


shared_ptr<PairFile<Matrix, Dvec> > CPCASSCF::solve() const {

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
  shared_ptr<const DFDist> df = ref_->geom()->df();
  shared_ptr<const DFHalfDist> half = df->compute_half_transform(ocoeff, nocca)->apply_J();
  shared_ptr<const DFFullDist> fullb = half->compute_second_transform(ocoeff, nocca);

  // making denominator...
  shared_ptr<PairFile<Matrix, Dvec> > denom;
  const double core_energy = ref_->geom()->nuclear_repulsion() + fci_->core_energy();
  {
    shared_ptr<Matrix> d0(new Matrix(*eig_));
    for (int i = 0; i != d0->size(); ++i)
      if (::fabs(d0->data(i)) < 1.0e-10) d0->data(i) = 1.0;
    shared_ptr<const Civec> d1_tmp(new Civec(*fci_->denom()));
    shared_ptr<Dvec>  d1(new Dvec(d1_tmp, ref_->nstate()));
    for (int i = 0; i != ref_->nstate(); ++i)
      *d1->data(i) -= fci_->energy(i) - core_energy;
#if 1
    // TODO understand this factor of 2
    *d1 *= 2;
#endif
    denom = shared_ptr<PairFile<Matrix, Dvec> >(new PairFile<Matrix, Dvec>(d0, d1));
  }

  // BFGS update of the denominator above
#if 0
  shared_ptr<BFGS<PairFile<Matrix, Dvec> > > bfgs(new BFGS<PairFile<Matrix, Dvec> >(denom, false));
#else
  shared_ptr<BFGS<PairFile<Matrix, Dvec> > > bfgs(new BFGS<PairFile<Matrix, Dvec> >(denom));
#endif


  // CI vector
  shared_ptr<PairFile<Matrix, Dvec> > source(new PairFile<Matrix, Dvec>(*grad_));
  // antisymmetrize
  source->first()->antisymmetrize();

  // divide by weight
  for (int ij = 0; ij != source->second()->ij(); ++ij) {
    source->second()->data(ij)->scale(1.0/fci_->weight(ij));
  }
  // project out Civector from the gradient
  source->second()->project_out(civector_);
  shared_ptr<LinearRM<PairFile<Matrix, Dvec> > > solver(new LinearRM<PairFile<Matrix, Dvec> >(CPHF_MAX_ITER, source));

  // initial guess
  shared_ptr<PairFile<Matrix, Dvec> > z = source->clone();
  z->zero();

  z = bfgs->extrapolate(source, z);
// not needed as z->xecond() is zero
//z->second()->project_out(civector_);

  // inverse matrix of C
//shared_ptr<Matrix> cinv(new Matrix(*ref_->coeff())); cinv->inverse();
  shared_ptr<const Matrix> ovl(new Overlap(geom_));
  shared_ptr<const Matrix> cinv(new Matrix(*ref_->coeff() % *ovl));

  // State averaged density matrix
  shared_ptr<const Matrix> dsa = ref_->rdm1_mat();

  cout << "  === CPCASSCF iteration ===" << endl << endl;

  // during the iteration, we need to
  //  (i) antisymmetrize Z, as well as
  //  (ii) project out c from z

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
    const double norm = sqrt(z->ddot(*z));
    cout << setw(4) <<  iter << " " << setprecision(14) << norm << endl;

    shared_ptr<const Matrix> z0 = z->first();
    shared_ptr<const Dvec>     z1 = z->second();

    // TODO duplicated operation of <I|H|z>. Should be resolved at the end.
    // only here we need to have det_ instead of detex
    shared_ptr<Matrix> sigmaorb = compute_amat(z1, civector_, detex);

    // computation of Atilde. Will be separated.
    // TODO index transformation can be skipped by doing so at the very end...
    shared_ptr<Matrix> cz0(new Matrix(*ref_->coeff() * *z0));
    shared_ptr<Matrix> cz0cinv(new Matrix(*ref_->coeff() * *z0 * *cinv));

    // [G_ij,kl (kl|D)] [(D|jS)+(D|Js)]   (capital denotes a Z transformed index)
    // (D|jx) -> (D|jS)
    {
      shared_ptr<DFFullDist> tmp0 = half->compute_second_transform(cz0cinv->data(), nbasis);
      shared_ptr<const DFHalfDist> tmp1 = df->compute_half_transform(cz0->data(), nocca)->apply_J();
      tmp0->daxpy(1.0, tmp1);
      shared_ptr<const DFFullDist> fulld = fullb->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
      shared_ptr<const Matrix> buf = tmp0->form_2index(fulld, 2.0); // Factor of 2
      dgemm_("T", "N", nbasis, nocca, nbasis, 1.0, ocoeff, nbasis, buf->data(), nbasis, 1.0, sigmaorb->data(), nbasis);
    }
    // [G_ij,kl (Kl|D)+(kL|D)] (D|sj)
    shared_ptr<DFFullDist> fullz = half->compute_second_transform(cz0->data(), nocca);
    fullz->symmetrize();
    {
      shared_ptr<const DFFullDist> tmp = fullz->apply_2rdm(ref_->rdm2_av()->data(), ref_->rdm1_av()->data(), nclosed, nact);
      shared_ptr<const Matrix> buf = half->form_2index(tmp, 2.0); // Factor of 2
      // mo transformation of s
      dgemm_("T", "N", nbasis, nocca, nbasis, 1.0, ocoeff, nbasis, buf->data(), nbasis, 1.0, sigmaorb->data(), nbasis);
    }

    // one electron part...
    shared_ptr<Matrix> htilde(new Matrix(*cz0 % *ref_->hcore() * *ref_->coeff()));
    htilde->symmetrize();
    *htilde *= 2.0;
    dgemm_("N", "N", nbasis, nocca, nocca, 2.0, htilde->data(), nbasis, dsa->data(), nbasis, 1.0, sigmaorb->data(), nbasis);

    sigmaorb->antisymmetrize();

    // At this point
    // htilde = Zh + hZ^dagger
    // fullb  = (D|ij)
    // fullz  = (D|ir)Z_rj + (D|rj)Z_ri

    // internal core fock operator...
    // [htilde + (kl|D)(D|ij) (2delta_ij - delta_ik)]_active

    // TODO this is a reference implementation
    // first form 4 index
    unique_ptr<double[]> buf = fullz->form_4index(fullb, 1.0);
    // TODO Awful code. To be updated. making the code that works in the quickest possible way
    // index swap
    unique_ptr<double[]> buf2(new double[nocca*nocca*nocca*nocca]);

    // bra ket symmetrization
    for (int i = 0; i != nocca*nocca; ++i)
      for (int j = 0; j != nocca*nocca; ++j)
        buf2[j+nocca*nocca*i] = buf[j+nocca*nocca*i] + buf[i+nocca*nocca*j];

    unique_ptr<double[]> Htilde2(new double[nact*nact*nact*nact]);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii)
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj)
        for (int k = nclosed, kk = 0; k != nocca; ++k, ++kk)
          for (int l = nclosed, ll = 0; l != nocca; ++l, ++ll)
            Htilde2[ll+nact*(kk+nact*(jj+nact*ii))] = buf2[l+nocca*(k+nocca*(j+nocca*i))];

    unique_ptr<double[]> Htilde1(new double[nact*nact]);
    for (int i = nclosed, ii = 0; i != nocca; ++i, ++ii) {
      for (int j = nclosed, jj = 0; j != nocca; ++j, ++jj) {
        Htilde1[jj+nact*ii] = htilde->element(j,i);
        for (int k = 0; k != nclosed; ++k)
          Htilde1[jj+nact*ii] += 2.0*buf2[k+nocca*(k+nocca*(j+nocca*i))] - buf2[k+nocca*(i+nocca*(j+nocca*k))];
      }
    }
    // factor of 2 in the equation
    dscal_(nact*nact,           2.0, Htilde1.get(), 1);
    dscal_(nact*nact*nact*nact, 2.0, Htilde2.get(), 1);

    shared_ptr<MOFile> top(new Htilde(ref_, nclosed, nocca, move(Htilde1), move(Htilde2)));
    vector<int> tmp(z1->ij(), 0);
    shared_ptr<Dvec> sigmaci = fci_->form_sigma(civector_, top, tmp);

    *sigmaci += *fci_->form_sigma(z1, fci_->jop(), tmp);

    for (int i = 0; i != z1->ij(); ++i)
      for (int j = 0; j != z1->data(i)->size(); ++j)
        sigmaci->data(i)->data(j) -= (fci_->energy(i) - core_energy) * z1->data(i)->data(j);

    sigmaci->project_out(civector_);

    shared_ptr<PairFile<Matrix, Dvec> > sigma(new PairFile<Matrix, Dvec>(sigmaorb, sigmaci));

    z = solver->compute_residual(z, sigma);

    z = bfgs->extrapolate(z, solver->civec());
    z->second()->project_out(civector_);

    if (sqrt(z->ddot(*z)) < CPHF_THRESH) break;

  }

solver->civec()->second()->print(-1);
  return solver->civec();

}


// computes A matrix (scaled by 2 here)
shared_ptr<Matrix> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec, shared_ptr<const Dvec> dvec, shared_ptr<const Determinants> o) const {

  const size_t nbasis = geom_->nbasis();
  shared_ptr<Matrix> amat(new Matrix(nbasis, nbasis));
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  const double* const coeff = ref_->coeff()->data();
  const double* const acoeff = coeff + nclosed*nbasis;

  // compute RDMs
  shared_ptr<const RDM<1> > rdm1t;
  shared_ptr<const RDM<2> > rdm2t;
  tie(rdm1t, rdm2t) = fci_->compute_rdm12_av_from_dvec(zvec, dvec, o);

  // symmetrize
  shared_ptr<RDM<1> > rdm1 = rdm1t->clone();
  shared_ptr<RDM<2> > rdm2 = rdm2t->clone();
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
  unique_ptr<double[]> buf(new double[nbasis*nact]);
  unique_ptr<double[]> buf2(new double[nbasis*nact]);
  dgemm_("N", "N", nbasis, nact, nbasis, 1.0, core_fock->data(), nbasis, acoeff, nbasis, 0.0, buf.get(), nbasis);
  dgemm_("N", "N", nbasis, nact, nact, 1.0, buf.get(), nbasis, rdm1->data(), nact, 0.0, buf2.get(), nbasis);
  dgemm_("T", "N", nbasis, nact, nbasis, prefactor, coeff, nbasis, buf2.get(), nbasis, 0.0, amat->element_ptr(0,nclosed), nbasis);

  // Half transformed DF vector
#if 0
  shared_ptr<const DFHalfDist> half = fci_->jop()->mo2e_1ext();
#else
  shared_ptr<const DFHalfDist> half = ref_->geom()->df()->compute_half_transform(acoeff, nact);
#endif
  shared_ptr<const DFFullDist> full = half->compute_second_transform(acoeff, nact)->apply_JJ();
  shared_ptr<const DFFullDist> fulld = full->apply_2rdm(rdm2->data());
  shared_ptr<const Matrix> jd = half->form_2index(fulld, 1.0);
  dgemm_("T", "N", nbasis, nact, nbasis, prefactor, coeff, nbasis, jd->data(), nbasis, 1.0, amat->element_ptr(0,nclosed), nbasis);

  return amat;
}
