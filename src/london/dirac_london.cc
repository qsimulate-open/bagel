//
// BAGEL - Parallel electron correlation program.
// Filename: dirac_london.cc
// Copyright (C) 2014 Toru Shiozaki
//
// Author: Ryan D. Reynolds <rreynoldschem@u.northwestern.edu>
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

#include <src/london/dirac_london.h>

using namespace bagel;
using namespace std;

BOOST_CLASS_EXPORT_IMPLEMENT(Dirac_London)

void Dirac_London::compute() {
  cout << "  Dirac-Fock method with London orbitals to be implemented soon!" << endl;
}


shared_ptr<const Reference> Dirac_London::conv_to_ref() const {
//  auto out = make_shared<Reference>(geom_, coeff(), nocc(), 0, coeff_->mdim()-nocc(), energy());
//  out->set_eig(eig_);
//  return out;
  return nullptr;
}
