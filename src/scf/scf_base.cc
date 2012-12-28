//
// BAGEL - Parallel electron correlation program.
// Filename: scf_base.cc
// Copyright (C) 2009 Toru Shiozaki
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


#include <src/scf/scf_base.h>
#include <src/rysint/eribatch.h>
#include <src/util/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;


SCF_base::SCF_base(const multimap<string, string>& idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : idata_(idat), geom_(geom), overlap_(new Overlap(geom)), hcore_(new Hcore(geom)) {

  unique_ptr<double[]> eig(new double[geom_->nbasis()]);
  eig_ = move(eig);
  hcore_->fill_upper();

  max_iter_ = read_input<int>(idata_, "maxiter", 100);
  max_iter_ = read_input<int>(idata_, "maxiter_scf", max_iter_);
  diis_start_ = read_input<int>(idata_, "diis_start", 1);
  thresh_overlap_ = read_input<double>(idata_, "thresh_overlap", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh", 1.0e-8);
  thresh_scf_ = read_input<double>(idata_, "thresh_scf", thresh_scf_);
  string dd = read_input<string>(idata_, "diis", "gradient");

  // so far assuming that this is RHF
  int nact = read_input<int>(idata_, "nact", 0);
  nocc_ = read_input<int>(idata_, "nocc", (geom_->nele()+nact)/2);
  noccB_ = nocc_ - nact;

  if (nocc_+noccB_ != geom_->nele()) throw runtime_error("nocc and nact are not consistently specified");

  tildex_ = shared_ptr<TildeX>(new TildeX(overlap_, thresh_overlap_));

  init_schwarz();

  // if ref is passed to this
  if (re != nullptr) {
    SCF_base::coeff_ = re->coeff();
  }
}


void SCF_base::init_schwarz() {
  schwarz_ = geom_->schwarz();
}

