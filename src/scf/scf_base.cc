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
#include <src/util/timer.h>
#include <src/math/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;


SCF_base::SCF_base(const shared_ptr<const PTree> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re, const bool need_schwarz)
 : idata_(idat), geom_(geom) {

  Timer scfb;
  overlap_ = make_shared<const Overlap>(geom);
  scfb.tick_print("Overlap matrix");
  hcore_ = make_shared<const Hcore>(geom);
  scfb.tick_print("Hcore matrix");

  eig_= unique_ptr<double[]>(new double[geom_->nbasis()]);

  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata_->get<int>("diis_start", 1);
  diis_size_ = idata_->get<int>("diis_size", 5);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh_scf", thresh_scf_);
  string dd = idata_->get<string>("diis", "gradient");

  const int ncharge = idata_->get<int>("charge", 0);
  const int nact    = idata_->get<int>("nact", (geom_->nele()-ncharge)%2);
  nocc_ = idata_->get<int>("nocc", (geom_->nele()-ncharge+nact)/2);
  noccB_ = nocc_ - nact;

  if (nocc_+noccB_ != geom_->nele()-ncharge) throw runtime_error("nocc and nact are not consistently specified");

  tildex_ = make_shared<TildeX>(overlap_, thresh_overlap_);

  scfb.tick_print("Overlap orthog");

  if (need_schwarz) {
    init_schwarz();
    scfb.tick_print("Schwarz matrix");
  }

  // if ref is passed to this
  if (re != nullptr) {
    coeff_ = re->coeff();
  }
  cout << endl;
}


void SCF_base::init_schwarz() {
  schwarz_ = geom_->schwarz();
}

