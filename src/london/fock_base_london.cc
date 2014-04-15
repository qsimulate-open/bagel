//
// BAGEL - Parallel electron correlation program.
// Filename: fock_base_london.cc
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


#include <src/london/fock_london.h>
#include <src/scf/symmat.h>

using namespace std;
using namespace bagel;

BOOST_CLASS_EXPORT_IMPLEMENT(Fock_base_London)

Fock_base_London::Fock_base_London(const shared_ptr<const Geometry_London> geom, const shared_ptr<const ZMatrix> previous, const std::shared_ptr<const ZMatrix> den, const vector<double>& schwarz)
 : ZMatrix1e(geom), cgeom_(geom), previous_(previous), density_(den), schwarz_(schwarz) {

  schwarz_thresh_ = geom->schwarz_thresh();

  init(geom); // zero here
}


void Fock_base_London::fock_one_electron_part() {

  //const int nbasis = ndim_;
  assert(ndim_ == mdim_);

  const int nirrep = cgeom_->nirrep();
  if (nirrep != 1) throw std::runtime_error("London-based methods currently cannot make use of symmetry.");
  #if 0
    ZMatrix intermediate(nbasis, nbasis);
    for (int i = 1; i != nirrep; ++i) {
      SymMat symm(cgeom_, i);
      ZMatrix tmp = symm % (*this) * symm;
      intermediate += tmp;
    }
    *this += intermediate;
    *this *= 1.0/nirrep;
  #endif

  *this += *previous_;

  fill_upper_conjg();
}


void Fock_base_London::computebatch(const array<shared_ptr<const Shell>,2>& input, const int offsetb0, const int offsetb1, shared_ptr<const Molecule>) {

  // input = [b1, b0]
  assert(input.size() == 2);
  const int dimb0 = input[1]->nbasis();
  const int dimb1 = input[0]->nbasis();

  for (int i = offsetb0; i != dimb0 + offsetb0; ++i) {
    for (int j = offsetb1; j != dimb1 + offsetb1; ++j) {
      data_[i*ndim_+j] = 0.0;
    }
  }
}


