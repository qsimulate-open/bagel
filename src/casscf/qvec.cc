//
// Newint - Parallel electron correlation program.
// Filename: qvec.cc
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


#include <src/casscf/qvec.h>

using namespace std;

Qvec::Qvec(const int n, const int m, shared_ptr<const DensityFit> df, shared_ptr<Coeff> coeff, const size_t nclosed, shared_ptr<FCI> fci)
 : QFile(n,m) {

  const int nbasis = df->nbasis0();
  assert(df->nbasis0() == df->nbasis1());

  shared_ptr<DF_Half> half = fci->jop()->mo2e_1ext();

  shared_ptr<DF_Full> full = half->compute_second_transform(coeff->data()+nclosed*nbasis, m)->apply_J()->apply_J(); 

  shared_ptr<DF_Full> prdm = full->apply_2rdm(fci->rdm2_av()->data());

  unique_ptr<double[]> tmp(new double[nbasis*m]);
  half->form_2index(tmp, prdm);
  dgemm_("T", "N", n, m, nbasis, 1.0, coeff->data(), nbasis, tmp.get(), nbasis, 0.0, data(), n);

}
