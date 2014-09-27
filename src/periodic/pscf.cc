//
// BAGEL - Parallel electron correlation program.
// Filename: pscf.cc
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


#include <iomanip>
#include <algorithm>
#include <src/util/timer.h>
#include <src/periodic/pscf.h>

using namespace std;
using namespace bagel;

PSCF::PSCF(const shared_ptr<const PTree> idata, const shared_ptr<const Geometry> geom, const shared_ptr<const Reference> re)
  : PSCF_base(idata, geom, re), dodf_(idata->get<bool>("df",true)) {
  cout << "  *** Periodic Hartree--Fock ***" << endl << endl;
  if (dodf_)
    throw runtime_error("Periodic code does not work with density fitting yet!");

  lattice_->print_primitive_vectors();
  lattice_->print_lattice_coordinates();
}

void PSCF::compute() {

  Timer pscftime;

  cout << indent << "=== Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->primitive_cell()->nuclear_repulsion() << endl << endl;

  cout << indent << "=== Lattice Nuclear Repulsion ===" << endl << indent << endl;
  cout << indent << fixed << setprecision(10) << setw(15) << lattice_->nuclear_repulsion() << endl << endl;

  energy_ = 0.0;
  const double error = 0.0;
  const int iter = 1;

  cout << indent << "=== Periodic RHF iteration (" + geom_->basisfile() + ") ===" << endl << indent << endl;
  cout << indent << setw(5) << iter << setw(20) << fixed << setprecision(8) << energy_ << "   "
                                    << setw(17) << error << setw(15) << setprecision(2) << pscftime.tick() << endl;

}
