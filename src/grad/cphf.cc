//
// BAGEL - Parallel electron correlation program.
// Filename: cphf.cc
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
#include <src/grad/cphf.h>
#include <src/util/matrix.h>

#define CPHF_MAX_ITER 100
#define CPHF_THRESH 1.0e-8

using namespace std;
using namespace bagel;

CPHF::CPHF(const shared_ptr<const Matrix> grad, const vector<double>& eig, const shared_ptr<const DFHalfDist> h,
           const shared_ptr<const Reference> r)
: solver_(new LinearRM<Matrix>(CPHF_MAX_ITER, grad)), grad_(grad), eig_(eig), halfjj_(h), ref_(r), geom_(r->geom()) {

}


shared_ptr<Matrix> CPHF::solve() const {

  const size_t naux = geom_->naux();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = geom_->nbasis() - nocca;

  const size_t nbasis = geom_->nbasis();

  shared_ptr<const Matrix> ocoeff = ref_->coeff()->slice(0, nocca);
  shared_ptr<const Matrix> vcoeff = ref_->coeff()->slice(nocca, nbasis);

  shared_ptr<Matrix> t(new Matrix(nbasis, nbasis));
  for (int i = 0; i != nocca; ++i)
    for (int a = nocca; a != nvirt+nocca; ++a)
      t->element(a,i) = grad_->element(a,i) / (eig_[a]-eig_[i]);

  Matrix jri(nbasis, nocca);
  Matrix jai(nvirt, nocca);
  Matrix kia(nocca, nvirt);

  cout << "  === CPHF iteration ===" << endl << endl;

  // TODO Max iter to be controlled by the input
  for (int iter = 0; iter != CPHF_MAX_ITER; ++iter) {
    solver_->orthog(t);

    shared_ptr<Matrix> sigma(new Matrix(nbasis, nbasis));
    // one electron part
    for (int i = 0; i != nocca; ++i)
      for (int a = nocca; a != nocca+nvirt; ++a)
        (*sigma)(a,i) = (eig_[a]-eig_[i]) * t->element(a,i);

    // J part
    Matrix pbmao(nbasis, nbasis);
    {
      Matrix pms(nocca, nbasis);
      dgemm_("T", "T", nocca, nbasis, nvirt, 1.0, t->element_ptr(nocca, 0), nbasis, vcoeff->data(), nbasis, 0.0, pms.data(), nocca);
      pbmao = *ocoeff * pms;
    }
    pbmao.symmetrize();
    jri = *geom_->df()->compute_Jop(pbmao.data()) * *ocoeff;
    jai = (*vcoeff % jri) * 4.0;

    // K part
    {
      // halfjj is an half transformed DF integral with J^{-1}_{DE}, given by the constructor
      shared_ptr<const Matrix> kir = halfjj_->compute_Kop_1occ(pbmao.data(), -2.0);
      kia = *kir * *vcoeff;
    }
    for (int i = 0; i != nocca; ++i)
      for (int a = 0; a != nvirt; ++a)
        (*sigma)(a+nocca,i) += jai[a+nvirt*i] + kia[i+nocca*a];

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
  return t;

}
