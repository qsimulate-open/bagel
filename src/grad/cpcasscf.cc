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

CPCASSCF::CPCASSCF(const shared_ptr<const PairFile<Matrix1e, Civec> > grad, const vector<double>& eig, const shared_ptr<DF_Half> h,
                   const shared_ptr<const Reference> r, const shared_ptr<const FCI> f)
: solver_(new Linear<PairFile<Matrix1e, Civec> >(CPHF_MAX_ITER, grad)), grad_(grad), eig_(eig), halfjj_(h), ref_(r), geom_(r->geom()), fci_(f) {

}


shared_ptr<PairFile<Matrix1e, Civec> > CPCASSCF::solve() const {

  const size_t naux = geom_->naux();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = geom_->nbasis() - nocca;

  const size_t nbasis = geom_->nbasis();

  const double* const ocoeff = ref_->coeff()->data();
  const double* const vcoeff = ocoeff + nocca*nbasis;

  int la, lb;
  tie(la, lb) = fci_->len_string();

  shared_ptr<Matrix1e> t0(new Matrix1e(geom_));
  shared_ptr<Civec> t1(new Civec(lb, la));
  shared_ptr<PairFile<Matrix1e, Civec> > t(new PairFile<Matrix1e, Civec>(t0, t1));

  shared_ptr<Matrix1e> d0(new Matrix1e(geom_));
  shared_ptr<const Civec> d1 = fci_->denom();
  for (int i = 0; i != nocca; ++i)
    for (int a = nocca; a != nvirt+nocca; ++a)
      d0->element(a,i) = 1.0 / (eig_[a]-eig_[i]);

  unique_ptr<double[]> jri(new double[nbasis*nocca]);
  unique_ptr<double[]> jai(new double[nvirt*nocca]);
  unique_ptr<double[]> kia(new double[nvirt*nocca]);

  cout << "  === CPCASSCF iteration ===" << endl << endl;

#if 0
  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
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

  }

  cout << endl;
  t = solver_->civec();
  t->fill_upper();
#endif
  return t;

}
