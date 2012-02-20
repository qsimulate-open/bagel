//
// Newint - Parallel electron correlation program.
// Filename: scf_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/scf_base.h>
#include <src/rysint/eribatch.h>
#include <src/util/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <boost/lexical_cast.hpp>

using namespace std;


SCF_base::SCF_base(multimap<string, string>& idat, const shared_ptr<Geometry> geom)
 : idata_(idat), geom_(geom), overlap_(new Overlap(geom)), hcore_(new Hcore(geom)) {

  unique_ptr<double[]> eig(new double[geom_->nbasis()]);
  eig_ = move(eig);
  hcore_->symmetrize();

  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  max_iter_ = read_input<int>(idata_, "maxiter_scf", max_iter_);
  diis_start_ = read_input<int>(idata_, "diis_start", 1);
  thresh_overlap_ = read_input<double>(idata_, "thresh_overlap", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh_scf", thresh_scf_);
  string dd = read_input<string>(idata_, "diis", "gradient");
  if (dd == "gradient") {
    density_change_ = false;
  } else if (dd == "density") {
    density_change_ = true;
  } else {
    throw runtime_error("unrecongnized option for DIIS error vectors");
  }

  { shared_ptr<TildeX>    tmp(new TildeX(overlap_, thresh_overlap_));    tildex_ = tmp; }

  init_schwarz();
}


void SCF_base::init_schwarz() {
  const vector<shared_ptr<Atom> > atoms = geom_->atoms(); 
  vector<shared_ptr<Shell> > basis; 
  for (auto aiter = atoms.begin(); aiter != atoms.end(); ++aiter) {
    const vector<shared_ptr<Shell> > tmp = (*aiter)->shells();
    basis.insert(basis.end(), tmp.begin(), tmp.end());  
  }
  const int size = basis.size();

  schwarz_.resize(size * size);
  for (int i0 = 0; i0 != size; ++i0) {
    const shared_ptr<Shell> b0 = basis[i0];
    for (int i1 = i0; i1 != size; ++i1) {
      const shared_ptr<Shell> b1 = basis[i1];

      vector<shared_ptr<Shell> > input;
      input.push_back(b1);
      input.push_back(b0);
      input.push_back(b1);
      input.push_back(b0);
      ERIBatch eribatch(input, 1.0);
      eribatch.compute();
      const double* eridata = eribatch.data();
      const int datasize = eribatch.data_size();
      double cmax = 0.0;
      for (int xi = 0; xi != datasize; ++xi, ++eridata) {
        const double absed = fabs(*eridata);
        if (absed > cmax) cmax = absed; 
      }
      schwarz_[i0 * size + i1] = cmax;
      schwarz_[i1 * size + i0] = cmax;
    }
  }
}
