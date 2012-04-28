//
// Newint - Parallel electron correlation program.
// Filename: mp2.cc
// Copyright (C) 2012 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki.toru@gmail.com>
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


// implements the MP2-F12 theory

#include <src/mp2/mp2.h>
#include <iostream>
#include <iomanip>
#include <src/util/f77.h>
#include <src/smith/prim_op.h>

using namespace std;

MP2::MP2(const multimap<string, string> input, const shared_ptr<Geometry> g, shared_ptr<Reference> r) : geom_(g), ref_(r) {
  cout << endl << "  === DF-MP2 calculation ===" << endl << endl;

  if (!geom_->df()) throw logic_error("MP2 is only implemented in DF");

}

void MP2::compute() {
  // TODO this factor of 2 is very much error-prone..
  const size_t nocc = geom_->nocc() / 2;
  const size_t nvirt = geom_->nbasis() - nocc;
  assert(geom_->nbasis() == ref_->coeff()->mdim());
  const size_t nbasis = geom_->nbasis();

  // first compute half transformed integrals
  shared_ptr<DF_Half> half = geom_->df()->compute_half_transform(ref_->coeff()->data(), nocc);  
  // second transform for virtual index
  // this is now (naux, nocc, nvirt)
  shared_ptr<DF_Full> full = half->compute_second_transform(ref_->coeff()->data()+nocc*nbasis, nvirt)->apply_J();

  cout << "  * 3-index integral transformation done" << endl << endl;

  // assemble
  unique_ptr<double[]> buf(new double[nocc*nvirt*nocc]);
  vector<double> eig = ref_->eig();

  // TODO in priciple this should run over occupied (for optimal implementations)...
  double sum = 0.0;
  for (size_t i = 0; i != nvirt; ++i) {
    // nocc * nvirt * nocc
    unique_ptr<double[]> data = full->form_4index(full, i); 
    copy(data.get(), data.get()+nocc*nvirt*nocc, buf.get());

    // using SMITH's symmetrizer (src/smith/prim_op.h)
    SMITH::sort_indices<2,1,0,2,1,-1,1>(data, buf, nocc, nvirt, nocc);
    double* tdata = buf.get();
    for (size_t j = 0; j != nocc; ++j) {
      for (size_t k = 0; k != nvirt; ++k) {
        for (size_t l = 0; l != nocc; ++l, ++tdata) {
          *tdata /= -eig[i+nocc]+eig[j]-eig[k+nocc]+eig[l]; 
        }
      }
    }
    sum += ddot_(nocc*nvirt*nocc, data, 1, buf, 1);
  }
  cout << "  MP2 correlation energy: " << fixed << setw(15) << setprecision(10) << sum << endl;
}
