//
// BAGEL - Parallel electron correlation program.
// Filename: zcasscf_debug.cc
// Copyright (C) 2014 Toru Shiozaki
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

#include <src/zcasscf/zcasscf.h>

using namespace std;
using namespace bagel;

static const int morbital = 2;
static const int norbital = 3;

void ZCASSCF::___debug___orbital_rotation(const bool kramers) {
  // currently transforming n and mth orbitals;
  const double angle = idata_->get<double>("debugrot", 0.0);
  if (angle == 0.0) return;

  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << "perturbing an orbital by : " << setprecision(5) << angle << endl;
  cout << "kramers adaptation       : " << (kramers ? "on" : "off") << endl;

  shared_ptr<ZRotFile> atmp = make_shared<ZRotFile>(nclosed_*2, nact_*2, nvirt_*2, /*superci*/false);

  atmp->ele_vc(morbital, norbital) = angle;
  if (kramers)
    atmp->ele_vc(nvirt_+morbital, nclosed_+norbital);

  shared_ptr<ZMatrix> amattmp = atmp->unpack<ZMatrix>();

  // multiply -1 from the formula. multiply -i to make amat hermite (will be compensated)
  *amattmp *= -1.0 * complex<double>(0.0, -1.0);

  // restore the matrix from RotFile
  unique_ptr<double[]> eig(new double[amattmp->ndim()]);
  amattmp->diagonalize(eig.get());

  auto amattmp_sav = amattmp->copy();
  for (int i = 0; i != amattmp->ndim(); ++i) {
    complex<double> ex = exp(complex<double>(0.0, eig[i]));
    for_each(amattmp->element_ptr(0,i), amattmp->element_ptr(0,i+1), [&ex](complex<double>& a) { a *= ex; });
  }

  auto expatmp = make_shared<ZMatrix>(*amattmp ^ *amattmp_sav);

  if (kramers)
    kramers_adapt(expatmp);

  shared_ptr<const ZMatrix> pco = make_shared<const ZMatrix>(*coeff_ * *expatmp);
  auto tt = make_shared<ZMatrix>(*coeff_ - *pco);

  cout << "norm delta(coeff) = " << setprecision(8) << tt->norm() << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl;

  coeff_ = pco;
}


void ZCASSCF::___debug___print_gradient(shared_ptr<const ZRotFile> grad) const {
  cout << ">>>>>>>>>>>> debug >>>>>>>>>>>>" << endl;
  cout << setprecision(10) << grad->ele_vc(morbital, norbital) << endl;
  cout << "<<<<<<<<<<<< debug <<<<<<<<<<<<" << endl;
}
