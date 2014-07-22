//
// BAGEL - Parallel electron correlation program.
// Filename: scf_base_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <RyanDReynolds@u.northwestern.edu>
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


#include <src/london/scf_base_london.h>
#include <src/london/reference_london.h>
#include <src/util/timer.h>
#include <src/math/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(SCF_base_London)

SCF_base_London::SCF_base_London(const shared_ptr<const PTree> idat, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re, const bool need_schwarz)
 : Method(idat, geom, re), eig_(cgeom_->nbasis()) {

  // if this is called by Opt
  do_grad_ = idata_->get<bool>("gradient", false);
  // enable restart capability
  restart_ = idata_->get<bool>("restart", false);

  Timer scfb;
  overlap_ = make_shared<const ZOverlap>(geom);
  scfb.tick_print("Overlap matrix");
  hcore_ = make_shared<const ZHcore>(geom);
  scfb.tick_print("Hcore matrix");

  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata_->get<int>("diis_start", 1);
  diis_size_ = idata_->get<int>("diis_size", 5);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh_scf", thresh_scf_);
  string dd = idata_->get<string>("diis", "gradient");

  multipole_print_ = idata_->get<int>("multipole", 1);

  const int ncharge = idata_->get<int>("charge", 0);
  const int nact    = idata_->get<int>("nact", (geom_->nele()-ncharge)%2);
  nocc_ = idata_->get<int>("nocc", (geom_->nele()-ncharge+nact)/2);
  noccB_ = nocc_ - nact;

  if (nocc_+noccB_ != geom_->nele()-ncharge) throw runtime_error("nocc and nact are not consistently specified");

  tildex_ = overlap_->tildex(thresh_overlap_);

  scfb.tick_print("Overlap orthog");

  if (need_schwarz) {
    init_schwarz();
    scfb.tick_print("Schwarz matrix");
  }

  // if ref is passed to this
  if (re != nullptr) {
    auto cref = dynamic_pointer_cast<const Reference_London>(re);
    assert(cref);
    coeff_ = cref->zcoeff();
  }

  cout << endl;
}


void SCF_base_London::init_schwarz() {
  schwarz_ = geom_->schwarz();
}

