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
#include <src/smith/prim_op.h>

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

  // in principle this is redundant, but cheap..
  unordered_map<bitset<2>, shared_ptr<const RelDFFull>> full = RelMOFile::compute_full(rocoeff, iocoeff, half_coulomb, /*apply_J*/true);  
  assert(full.size() == 4);
  unordered_map<bitset<2>, shared_ptr<RelDFFull>> full_d;
  for (auto& i : full)
    full_d.insert(make_pair(i.first, i.second->clone()));

  for (auto& t : full_d) {
    for (auto& s : full) {
      bitset<4> b;
      b[0] = s.first[0];
      b[1] = t.first[0];
      b[2] = s.first[1];
      b[3] = t.first[1];
      // t^+ s^+ t s
      shared_ptr<const ZMatrix> rdmbuf = fci->jop()->mo2e(b);
      // t^+ t s^+ s
      shared_ptr<ZMatrix> rdm = rdmbuf->clone();
      assert(rdm->ndim() == rdm->mdim() && rdm->ndim() == nocc*nocc);
      SMITH::sort_indices<0,2,1,3,0,1,1,1>(rdmbuf->data(), rdm->data(), nocc, nocc, nocc, nocc); 
// TODO implement
//    *t.second += *s->apply_2rdm(rdm);
    }
  }

}
