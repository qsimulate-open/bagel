//
// BAGEL - Parallel electron correlation program.
// Filename: pscf_base.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Hai-Anh Le <anh@u.northwestern.edu>
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


#include <src/periodic/pscf_base.h>
#include <src/util/timer.h>
#include <src/util/math/diis.h>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(PSCF_base)

PSCF_base::PSCF_base(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
 : Method(idata, geom, re), dodf_(idata->get<bool>("df", true)) {

  Timer pscf;

  if (geom->dofmm()) {
    const int lmax   = idata->get<int>("l_max", 10);
    const int ws     = idata->get<int>("ws", 2);
    const double beta   = idata->get<double>("beta", 1.0);
    const int height = idata->get<int>("height", 21);
    const int k_param = idata->get<int>("k_parameter", 9);
    lattice_ = make_shared<const Lattice>(geom, k_param, ws, true, make_tuple(height, lmax, idata->get<bool>("contract", true), dodf_, idata->get<double>("thresh_fmm", PRIM_SCREEN_THRESH)));
    const bool doewald    = idata->get<bool>("ewald", false);
    form_pfmm(make_tuple(lmax, ws, beta, doewald, idata->get<int>("extent", 10)));
  } else {
    lattice_ = make_shared<const Lattice>(geom, idata->get<int>("extent", 0));
  }

  eig_.resize(lattice_->num_lattice_kvectors());
  for (auto& eigblock : eig_) eigblock = make_shared<VectorB>(geom->nbasis());

  restart_ = idata_->get<bool>("restart", false);
  overlap_ = make_shared<const POverlap>(lattice_);
  koverlap_ = overlap_->ft(lattice_->lattice_vectors(), lattice_->lattice_kvectors());
  pscf.tick_print("Periodic overlap matrix");
  hcore_ = make_shared<const PHcore>(lattice_);
  pscf.tick_print("Periodic hcore matrix");

  max_iter_ = idata_->get<int>("maxiter", 100);
  max_iter_ = idata_->get<int>("maxiter_scf", max_iter_);
  diis_start_ = idata_->get<int>("diis_start", 1);
  diis_size_ = idata_->get<int>("diis_size", 5);
  thresh_overlap_ = idata_->get<double>("thresh_overlap", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh", 1.0e-8);
  thresh_scf_ = idata_->get<double>("thresh_scf", thresh_scf_);
  string dd = idata_->get<string>("diis", "gradient");

  const int ncharge = idata_->get<int>("charge", 0);
  const int nact    = idata_->get<int>("nact", (geom_->nele() - ncharge) % 2);
  nocc_ = idata_->get<int>("nocc", (geom_->nele() - ncharge + nact) / 2);
  noccB_ = nocc_ - nact;

  if (nocc_ + noccB_ != geom_->nele() - ncharge)
    throw runtime_error("*** nocc and nact for the unit cell are not consistently specified!");

  ktildex_ = koverlap_->tildex(thresh_overlap_);
  pscf.tick_print("Periodic overlap orthogonalization");

  cout << endl;
}


void PSCF_base::form_pfmm(tuple<int, int, double, bool, int> fmmp) {

  // rectangular cells for now
  cout << "  PFMM option is specified: simulation cell will be constructed." << endl;
  Timer time;
  const int ndim = lattice_->ndim();
  bool is_rec = true;
  for (int i = 0; i != ndim; ++i)
    for (int j = 0; j != ndim; ++j) {
      if (i != j) {
        double dp = 0.0;
        for (int k = 0; k != 3; ++k)
          dp += lattice_->primitive_cell()->primitive_vectors(i)[k] * lattice_->primitive_cell()->primitive_vectors(j)[k];
        if (dp > numerical_zero__) {
          is_rec = false;
          break;
        }
      }
    }

  if (is_rec) {
    cout << "  Unit cell is rectangular, simulation cell is the same as primitive cell." << endl;
  } else {
    cout << "  Unit cell is non-rectangular, simulation cell is the smallest cubic cell that encloses the unit cell." << endl;
    throw runtime_error("  ***  Non-cubic cell under contruction... Oops sorry!");
  }

  fmm_ = make_shared<const PFMM>(lattice_, fmmp, dodf_);
  cout << "        elapsed time:  " << setw(10) << setprecision(2) << time.tick() << " sec." << endl << endl;
}
