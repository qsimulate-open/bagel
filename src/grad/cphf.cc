//
// BAGEL - Brilliantly Advanced General Electronic Structure Library
// Filename: cphf.cc
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

#include <src/grad/cphf.h>
#include <src/util/math/linearRM.h>

using namespace std;
using namespace bagel;

CPHF::CPHF(const shared_ptr<const Matrix> grad, const VectorB& eig, const shared_ptr<const DFHalfDist> h,
           const shared_ptr<const Reference> r)
: grad_(grad), eig_(eig), halfjj_(h), ref_(r), geom_(r->geom()) {

}


shared_ptr<Matrix> CPHF::solve(const double zthresh, const int zmaxiter) {

  LinearRM<Matrix> solver(zmaxiter, grad_);

  const size_t nmobasis = ref_->coeff()->mdim();
  const size_t nocca = ref_->nocc();
  const size_t nvirt = nmobasis - nocca;

  const MatView ocoeff = ref_->coeff()->slice(0, nocca);
  const MatView vcoeff = ref_->coeff()->slice(nocca, nmobasis);

  auto t = make_shared<Matrix>(nmobasis, nmobasis);
  for (int i = 0; i != nocca; ++i)
    for (int a = nocca; a != nvirt+nocca; ++a)
      t->element(a,i) = grad_->element(a,i) / (eig_(a)-eig_(i));
  t->scale(1.0/t->norm());

  cout << "  === Z-vector iteration ===" << endl << endl;

  Timer timer;
  for (int iter = 0; iter != zmaxiter; ++iter) {

    auto sigma = make_shared<Matrix>(nmobasis, nmobasis);
    // one electron part
    for (int i = 0; i != nocca; ++i)
      for (int a = nocca; a != nocca+nvirt; ++a)
        (*sigma)(a,i) = (eig_(a)-eig_(i)) * t->element(a,i);

    // J part
    shared_ptr<const Matrix> tvo = t->get_submatrix(nocca, 0, nvirt, nocca);
    auto pbmao = make_shared<Matrix>(ocoeff ^ (vcoeff * *tvo));
    pbmao->symmetrize();
    Matrix jri = *geom_->df()->compute_Jop(pbmao) * ocoeff;
    Matrix jai = (vcoeff % jri) * 4.0;

    // K part
    // halfjj is an half transformed DF integral with J^{-1}_{DE}, given by the constructor
    shared_ptr<const Matrix> kir = halfjj_->compute_Kop_1occ(pbmao, -2.0);
    Matrix kia = *kir * vcoeff;

    for (int i = 0; i != nocca; ++i)
      for (int a = 0; a != nvirt; ++a)
        (*sigma)(a+nocca,i) += jai(a,i) + kia(i,a);

    t = solver.compute_residual(t, sigma);

    cout << setw(7) << iter << " " << setw(20) << setprecision(14) << t->rms() << setw(15) << setprecision(2) << timer.tick() << endl;
    if (t->rms() < zthresh) break;

    for (int i = 0; i != nocca; ++i)
      for (int a = nocca; a != nvirt+nocca; ++a)
        t->element(a,i) /= (eig_(a)-eig_(i));
    t->scale(1.0/t->norm());
  }

  cout << endl;
  t = solver.civec();
  t->fill_upper();
  return t;

}
