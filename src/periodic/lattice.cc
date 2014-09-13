//
// BAGEL - Parallel electron correlation program.
// Filename: lattice.cc
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


#include <src/periodic/lattice.h>

using namespace std;
using namespace bagel;

Lattice::Lattice(const shared_ptr<const Geometry> g) : unit_cell_(g) {

  ndim_ = g->lattice_vectors().size();
}


void Lattice::print_lattice_vectors() const {

  const string indent = "    ";
  cout << indent << setw(4) << " Basic lattice vectors:" << endl;
  for (int i = 0; i != ndim_; ++i)
    cout << indent << indent << setprecision(6) << "(" << setw(10) << unit_cell_->lattice_vectors(i)[0] << ", "
                                                       << setw(10) << unit_cell_->lattice_vectors(i)[1] << ", "
                                                       << setw(10) << unit_cell_->lattice_vectors(i)[2] << ") " << endl;
}
