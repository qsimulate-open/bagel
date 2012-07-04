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

#define CPHF_MAX_ITER 100
#define CPHF_THRESH 1.0e-8

using namespace std;

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix1e, Dvec> > grad, const vector<double>& eig, const shared_ptr<const DF_Half> h,
                   const shared_ptr<const DF_Half> h2, const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: solver_(new Linear<PairFile<Matrix1e, Dvec> >(CPHF_MAX_ITER, grad)), grad_(grad), eig_(eig), half_(h),
  halfjj_(h2), ref_(r), geom_(r->geom()), fci_(f) {

}


shared_ptr<PairFile<Matrix1e, Dvec> > CPCASSCF::solve() const {

  const size_t naux = geom_->naux();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = geom_->nbasis() - nocca;

  const size_t nbasis = geom_->nbasis();

  const double* const ocoeff = ref_->coeff()->data();
  const double* const vcoeff = ocoeff + nocca*nbasis;

  shared_ptr<Matrix1e> t0(new Matrix1e(geom_));
  shared_ptr<Dvec> t1(new Dvec(fci_->det(), ref_->nstate()));
  shared_ptr<PairFile<Matrix1e, Dvec> > t(new PairFile<Matrix1e, Dvec>(t0, t1));

  shared_ptr<Matrix1e> d0(new Matrix1e(geom_));
  for (int i = 0; i != nocca; ++i)
    for (int a = nocca; a != nvirt+nocca; ++a)
      d0->element(a,i) = eig_[a]-eig_[i];

  shared_ptr<const Civec> d1_each = fci_->denom();
  shared_ptr<const Dvec> d1(new Dvec(d1_each, ref_->nstate()));
  shared_ptr<Dvec> z1 = d1->clone();

  unique_ptr<double[]> jri(new double[nbasis*nocca]);
  unique_ptr<double[]> jai(new double[nvirt*nocca]);
  unique_ptr<double[]> kia(new double[nvirt*nocca]);

  cout << "  === CPCASSCF iteration ===" << endl << endl;


  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {

    shared_ptr<Matrix1e> amat = compute_amat(z1);
    amat->print();

#if 0
    solver_->orthog(t);

    shared_ptr<Matrix1e> sigma(new Matrix1e(geom_));
    // one electron part
    for (int i = 0; i != nocca; ++i)
      for (int a = nocca; a != nocca+nvirt; ++a)
        sigma->element(a,i) = (eig_[a]-eig_[i]) * t->element(a,i);

    // J part
    shared_ptr<Matrix1e> pbmao(new Matrix1e(geom_));
    {
      unique_ptr<double[]> pms(new double[nbasis*nocca]);
      dgemm_("T", "T", nocca, nbasis, nvirt, 1.0, t->element_ptr(nocca, 0), nbasis, vcoeff, nbasis, 0.0, pms.get(), nocca);
      dgemm_("N", "N", nbasis, nbasis, nocca, 1.0, ocoeff, nbasis, pms.get(), nocca, 0.0, pbmao->data(), nbasis);
    }
    pbmao->symmetrize();
    {
      unique_ptr<double[]> jrs = geom_->df()->compute_Jop(pbmao->data());
      dgemm_("N", "N", nbasis, nocca, nbasis, 1.0, jrs.get(), nbasis, ocoeff, nbasis, 0.0, jri.get(), nbasis);
    }
    dgemm_("T", "N", nvirt, nocca, nbasis, 4.0, vcoeff, nbasis, jri.get(), nbasis, 0.0, jai.get(), nvirt); 

    // K part
    {
      unique_ptr<double[]> kir = halfjj_->compute_Kop_1occ(pbmao->data());
      dgemm_("N", "N", nocca, nvirt, nbasis, -2.0, kir.get(), nocca, vcoeff, nbasis, 0.0, kia.get(), nocca); 
    }

    // Note that the alignment is different
    for (int i = 0; i != nocca; ++i)
      for (int a = 0; a != nvirt; ++a)
        sigma->element(a+nocca,i) += jai[a+nvirt*i] + kia[i+nocca*a];

    t = solver_->compute_residual(t, sigma);

    // TODO to be controlled by the input
    if (t->norm() < CPHF_THRESH) break;

    for (int i = 0; i != nocca; ++i)
      for (int a = nocca; a != nvirt+nocca; ++a)
        t->element(a,i) /= (eig_[a]-eig_[i]);

    cout << setw(6) << iter << setw(20) << setprecision(10) << t->norm() << endl;
#endif

  }

#if 0
  cout << endl;
  t = solver_->civec();
  t->fill_upper();
#endif
  return t;

}


// computes A matrix (scaled by 2 here)
shared_ptr<Matrix1e> CPCASSCF::compute_amat(shared_ptr<const Dvec> zvec) const {
  shared_ptr<Matrix1e> amat(new Matrix1e(ref_->geom())); 

  const size_t nbasis = geom_->nbasis();
  const int nclosed = ref_->nclosed();
  const int nact = ref_->nact();

  const double* const coeff = ref_->coeff()->data(); 
  const double* const acoeff = coeff + nclosed*nbasis;

  // compute RDMs 
  shared_ptr<const RDM<1> > rdm1 = ref_->rdm1_av();
  shared_ptr<const RDM<2> > rdm2 = ref_->rdm2_av();

  // core Fock operator
  shared_ptr<const Matrix1e> core_fock = fci_->jop()->core_fock();
  unique_ptr<double[]> buf(new double[nbasis*nact]);
  unique_ptr<double[]> buf2(new double[nbasis*nact]);
  dgemm_("N", "N", nbasis, nact, nbasis, 1.0, core_fock->data(), nbasis, acoeff, nbasis, 0.0, buf.get(), nbasis); 
  dgemm_("N", "N", nbasis, nact, nact, 1.0, buf.get(), nbasis, rdm1->data(), nact, 0.0, buf2.get(), nbasis); 
  dgemm_("T", "N", nbasis, nact, nbasis, 2.0, coeff, nbasis, buf2.get(), nbasis, 0.0, amat->element_ptr(0,nclosed), nbasis); 

amat->print();

  // Half transformed DF vector
  shared_ptr<const DF_Half> half = fci_->jop()->mo2e_1ext();
  shared_ptr<const DF_Full> full = half->compute_second_transform(acoeff, nact)->apply_JJ();
  shared_ptr<const DF_Full> fulld = full->apply_2rdm(rdm2->data());
  unique_ptr<double[]> jd = half->form_2index(fulld);
  dgemm_("T", "N", nbasis, nact, nbasis, 2.0, coeff, nbasis, jd.get(), nbasis, 1.0, amat->element_ptr(0,nclosed), nbasis); 

  return amat;
}
