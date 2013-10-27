//
// BAGEL - Parallel electron correlation program.
// Filename: zqvec.cc
// Copyright (C) 2013 Toru Shiozaki
//
// Author: Toru Shiozaki <shiozaki@northwestern.edu>
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

#include <src/zcasscf/zqvec.h>

using namespace std;
using namespace bagel;
    

ZQvec::ZQvec(const int n, const int m, shared_ptr<const Geometry> geom, shared_ptr<const ZMatrix> coeff, shared_ptr<const ZHarrison> fci, const bool gaunt, const bool breit)
 : ZMatrix(n,m) {

  assert(gaunt || !breit);
  if (gaunt) throw logic_error("Gaunt not implemented yet in ZQvec");

  array<shared_ptr<const ZMatrix>,2> kcoeff = fci->kramers_coeff();
  array<array<shared_ptr<const Matrix>,4>,2> rocoeff;
  array<array<shared_ptr<const Matrix>,4>,2> iocoeff;
  assert(geom->nbasis()*4 == kcoeff[0]->ndim());

  const int nocc = kcoeff[0]->mdim();

  for (int k = 0; k != 2; ++k) {
    for (int i = 0; i != 4; ++i) {
      shared_ptr<const ZMatrix> oc = kcoeff[k]->get_submatrix(i*geom->nbasis(), 0, geom->nbasis(), nocc);
      rocoeff[k][i] = oc->get_real_part();
      iocoeff[k][i] = oc->get_imag_part();
    }
  }

  array<list<shared_ptr<RelDFHalf>>,2> half_coulomb = fci->jop()->half_complex_coulomb();
//array<list<shared_ptr<RelDFHalf>>,2> half_gaunt   = fci->jop()->half_complex_gaunt();

  unordered_map<bitset<2>, shared_ptr<const RelDFFull>> full = RelMOFile::compute_full(rocoeff, iocoeff, half_coulomb, false);  

}
